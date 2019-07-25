# This script accepts ASCAT segmentation file (single aggregated file) and spread 
# - total copy number values
# - minor_allele_frequency values (nMinor / total, if total is 0, set to -1)
# on probe level, and then transform into matrix format 
# - row as probeid
# - column as tumor_aliquot
# Please note that beause this is run by chunk to reduce memory burden, there is no 
# guarantee all tumor_aliquot were finally included (final result could have less columns 
# simply because there is no data). It's unlikely to happen to autosomes, but does happen in chrY

library(data.table)
library(dplyr)
library(reshape2)

# read ASCAT metadata
meta <- fread("TCGA.TARGET.SNP6.ASCAT.11582.tsv")
meta <- meta[ascat_fail==F]
cols <- c("probeid", meta$tumor_aliquot)
chrs <- paste0("chr", c(1:22, "X", "Y"))

# Read all segments
seg <- fread("zcat TCGA.TARGET.SNP6.ASCAT.11582.segment.txt.gz")
# Rename
setnames(seg, c("tumor_aliquot", "chr", "start", "end", "na", "nb", "cn"))

# Some adjacent ASCAT segment share boundary point, do the following to support matrix transformation
# - remove that point from a early segment "end" if late segment is a 1-base segment
# - otherwise, remove that point from later segment "start"
# identify late segment index that share "start" with a early segment "end"
key.end <- with(seg, paste(tumor_aliquot, chr, end))
key.end <- c("", key.end[-length(key.end)])
key.start <- with(seg, paste(tumor_aliquot, chr, start))
w <- which(key.end == key.start)
# find those later segment with start = end
w2 <- which(with(seg[w], start == end))
# change values
seg$start[w[-w2]] <- seg$start[w[-w2]] + 1
seg$end[w[w2]-1] <- seg$end[w[w2]-1] - 1

# Calculate minor allele frequencies
# Note: from ASCAT output, nMinor sometimes larger than nMajor in chrY.
seg <- seg %>% 
  mutate(maf = ifelse(cn == 0, -1, nb / cn), 
         chr = factor(chr, levels = chrs)) %>% 
  select(-na, -nb) %>% 
  setDT() %>% 
  setkey(start, end)

# Read probeset position file
probe <- fread("zcat snp6.na35.remap.hg38.txt.gz")
probe <- probe %>% 
  mutate(start = pos, 
         end = pos) %>% 
  select(probeid, chr, start, end) %>% 
  setDT() %>%
  setkey(start, end)

# define chunk size
chunk.size <- 25000

for(i in 1:length(chrs)) {
  # Subset by chromosome
  chrom <- chrs[i]
  print(chrom)
  chr.seg <- seg[chr == chrom]
  chr.probe <- probe[chr == chrom]

  # Number of Chunks
  num.chunks <- ceiling(nrow(chr.probe) / chunk.size)

  for(chunk in 1:num.chunks) {
    # Subset chunks 
    my.probe <- chr.probe[((i-1)*chunk.size + 1) : min(i*chunk.size, nrow(chr.probe))]
    my.seg <- chr.seg[start <= max(my.probe$end) & end >= min(my.probe$start)]

    # Overlap segment with probes
    ov <- foverlaps(my.seg, my.probe, mult="all", type="any", nomatch=0L) %>%
      select(tumor_aliquot, probeid, maf, cn) 

    # transform total copy number values into matrix and write to disk
    cn <- acast(ov, probeid~tumor_aliquot, value.var="cn")
    cn <- cn[order(match(rownames(cn), my.probe$probeid)), order(match(colnames(cn), meta$tumor_aliquot))]
    cn.file <- paste0(chrom, ".chunk", chunk, ".ascat_copy_number.by_probe.tsv")
    write.table(cn, cn.file, col.names=T, row.names=T, sep="\t", quote=F)

    # transform MAF values into matrix and write to disk
    cn <- acast(ov, probeid~tumor_aliquot, value.var="maf")
    cn <- cn[order(match(rownames(cn), my.probe$probeid)), order(match(colnames(cn), meta$tumor_aliquot))]
    cn.file <- paste0(chrom, ".chunk", chunk, ".ascat_minor_allele_frequency.by_probe.tsv")
    write.table(cn, cn.file, col.names=T, row.names=T, sep="\t", quote=F)

  }
}
