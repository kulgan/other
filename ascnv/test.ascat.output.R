library(data.table)

setwd("~/SCRATCH/snp6cnv/test")

t.call = fread("./2f53bbf7-3b48-4092-8ac4-c22fb2bed7f0/birdseed.calls.txt")
names(t.call) = c("probeid", "t.call")

t.conf = fread("./2f53bbf7-3b48-4092-8ac4-c22fb2bed7f0/birdseed.confidences.txt")
names(t.conf)[2] = "t.conf"

t.lrrbaf = fread("./2f53bbf7-3b48-4092-8ac4-c22fb2bed7f0/2f53bbf7-3b48-4092-8ac4-c22fb2bed7f0.lrr_baf.txt")
names(t.lrrbaf) = c("probeid", "chr", "pos", "t.logr", "t.baf")

n.call = fread("./c23d61b6-4048-4aed-9647-fb920211890b/birdseed.calls.txt")
names(n.call) = c("probeid", "n.call")

n.conf = fread("./c23d61b6-4048-4aed-9647-fb920211890b/birdseed.confidences.txt")
names(n.conf)[2] = "n.conf"

n.lrrbaf = fread("./c23d61b6-4048-4aed-9647-fb920211890b/c23d61b6-4048-4aed-9647-fb920211890b.lrr_baf.txt")
names(n.lrrbaf) = c("probeid", "chr", "pos", "n.logr", "n.baf")

seg = fread("TARGET-ALL-P2.0161a290-9c04-59dc-a9c0-e027319d335d.ascat2.allelic_specific.seg.txt")
names(seg)[2:4] = c("chr", "start", "end")
setkey(seg, chr, start, end)

snp = fread("zcat /home/ubuntu/SCRATCH/snp6cnv/snp6.na35.remap.hg38.txt.gz")

# merge genotype, confidence and BAF together
geno = n.call
geno$t.call = t.call$t.call
geno$n.conf = n.conf$n.conf
geno$t.conf = t.conf$t.conf
geno = merge(geno, n.lrrbaf, by="probeid")
geno = merge(geno, t.lrrbaf, by="probeid")

# merge genotypes with snp location data
geno = merge(snp[, c("probeid", "chr", "pos"), with=F], geno, by="probeid")
geno$start = geno$pos
geno$end = geno$pos
setkey(geno, chr, start, end)

# add segmentation data
ov = foverlaps(geno, seg, type = "any", mult = "all", nomatch=0L, which = T)
geno$cn.b = seg$Minor_Copy_Number[ov$yid]
geno$cn.a = seg$Major_Copy_Number[ov$yid]
geno$cn = seg$Total_Copy_Number[ov$yid]
geno$seg.index = ov$yid

geno = geno[, -c("pos.x", "pos.y", "start", "end", "chr.x", "chr.y"), with=F]

saveRDS(geno, "combined.ascat.data.RDS")




library(data.table)
library(dplyr)
library(reshape2)

# read ASCAT metadata
meta = fread("TCGA.TARGET.SNP6.ASCAT.11582.tsv")
meta = meta[ascat_fail==F]
cols = c("probeid", meta$tumor_aliquot)
chrs = paste0("chr", c(1:22, "X", "Y"))

# Read all segments
seg = fread("zcat TCGA.TARGET.SNP6.ASCAT.11582.segment.txt.gz")
# Rename
setnames(seg, c("tumor_aliquot", "chr", "start", "end", "na", "nb", "cn"))
# Calculate minor allele frequencies
# Note: from ASCAT output, nMinor sometimes larger than nMajor in chrY.
seg = seg %>% 
  mutate(maf = ifelse(cn == 0, -1, nb / cn), 
         chr = factor(chr, levels = chrs)) %>% 
  select(-na, -nb) %>% 
  setDT() %>% 
  setkey(start, end)

# Read probeset position file
probe = fread("zcat snp6.na35.remap.hg38.txt.gz")
probe = probe %>% 
  mutate(start = pos, 
         end = pos) %>% 
  select(probeid, chr, start, end) %>% 
  setDT() %>%
  setkey(start, end)

chunk.size = 25000

for(i in 1:length(chrs)) {
  # Subset by chromosome
  chrom = chrs[i]
  print(chrom)
  chr.seg = seg[chr == chrom]
  chr.probe = probe[chr == chrom]

  # Number of Chunks
  num.chunks = ceiling(nrow(chr.probe) / chunk.size)

  for(chunk in 1:num.chunks) {
    # Subset chunks 
    my.probe = chr.probe[((i-1)*chunk.size + 1) : min(i*chunk.size, nrow(my.probe))]
    my.seg = chr.seg[start <= max(my.probe$end) & end >= min(my.probe$start)]

    # Overlap 
    ov = foverlaps(my.probe, my.seg, mult="all", type="any", nomatch=0L) %>%
      select(tumor_aliquot, probeid, maf, cn) 

    cn = acast(ov, probeid~tumor_aliquot, value.var="cn")
    cn = cn[order(match(rownames(cn), my.probe$probeid)), order(match(colnames(cn), meta$tumor_aliquot))]
    cn.file = paste0(chrom, ".chunk", chunk, ".ascat_copy_number.by_probe.tsv")
    write.table(cn, cn.file, col.names=T, row.names=T, sep="\t", quote=F)

    cn = acast(ov, probeid~tumor_aliquot, value.var="maf")
    cn = cn[order(match(rownames(cn), my.probe$probeid)), order(match(colnames(cn), meta$tumor_aliquot))]
    cn.file = paste0(chrom, ".chunk", chunk, ".ascat_minor_allele_frequency.by_probe.tsv")
    write.table(cn, cn.file, col.names=T, row.names=T, sep="\t", quote=F)

  }
}
