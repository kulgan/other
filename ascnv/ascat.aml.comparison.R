# This script collect SKCM SNP6 ASCAT output, and trying to detect probe positions that has significant differences
# in copy number gains or copy number losses, between solid primary tumor vs metastasis tumoe, using fisher exact tests, 
# with qvalue corrections
# - conclusion 1: fisher exact test p-values seem to be deflated 
# - conclusion 2: no probe is statistically significant. No permutation tests were run because preliminary results not promissing

# > summary(probe.summary$qvalue.loss)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8065  0.8636  0.9949  0.9463  1.0000  1.0000 
# > summary(probe.summary$qvalue.gain)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5753  0.5753  0.5753  0.6148  0.6414  0.7703 


library(data.table)
library(dplyr)

# read sample metadata
meta <- fread("GDC.SNP6.pairs.10878.02192024.tsv")
my.meta <- meta %>% 
  filter(project %in% c("TCGA-LAML", "TARGET-AML")) %>% 
  rename(type = project, aliquot = tumor_aliquot) %>%
  select(aliquot, type) %>%
  setDT()
chrs <- paste0("chr", 1:22)

# subsetting segment
seg <- fread("zcat TCGA.TARGET.SNP6.ASCAT.11582.1_off_modified.segment.txt.gz")
names(seg) <- c("aliquot", "chr", "start", "end", "na", "nb", "cn")
my.seg <- seg %>%
  filter(aliquot %in% my.meta$aliquot & chr %in% chrs) %>% 
  merge(my.meta, by = "aliquot") %>% 
  setDT() %>%
  setkey(chr, start, end)

# read probe location information
probe <- fread("zcat snp6.na35.remap.hg38.txt.gz")
probe <- probe %>% 
  filter(chr %in% chrs) %>%
  mutate(start = pos, 
         end = pos) %>% 
  select(probeid, chr, start, end) %>% 
  setDT() %>%
  setkey(chr, start, end)

# find overlapped probes with segments
ov <- foverlaps(my.seg, probe, mult="all", type="any", nomatch=0L) 
names(ov)[6:7] <- c("seg_start", "seg_end")

# collect counts for gain and loss
probe.summary <- ov %>% 
    group_by(probeid, chr, start) %>%
    summarise(kid_del = sum(cn == 0 & type == "TARGET-AML"),                      # both lost
              kid_loss_1 = sum(cn == 1 & type == "TARGET-AML"),                   # one lost, the other unchanged
              kid_gain_1 = sum(nb == 1 & na > 1 & type == "TARGET-AML"),          # one gain, the other unchanged
              kid_copy_neutral = sum(cn == 2 & type == "TARGET-AML"),             # copy number neutral (but could be 0/2)
              kid_copy_unchanged = sum(na == 1 & nb == 1 & type == "TARGET-AML"), # both unchanged (subset of copy number neutral)
              kid_gain = sum(cn > 2 & type == "TARGET-AML"),                      # copy number gain
              kid_loh = sum(nb == 0 & na > 0 & type == "TARGET-AML"),             # loh
              kid_total = sum(type == "TARGET-AML"),
              adult_del = sum(cn == 0 & type == "TCGA-LAML"),                      # both lost
              adult_loss_1 = sum(cn == 1 & type == "TCGA-LAML"),                   # one lost, the other unchanged
              adult_gain_1 = sum(nb == 1 & na > 1 & type == "TCGA-LAML"),          # one gain, the other unchanged
              adult_copy_neutral = sum(cn == 2 & type == "TCGA-LAML"),             # copy number neutral (but could be 0/2)
              adult_copy_unchanged = sum(na == 1 & nb == 1 & type == "TCGA-LAML"), # both unchanged (subset of copy number neutral)
              adult_gain = sum(cn > 2 & type == "TCGA-LAML"),                      # copy number gain
              adult_loh = sum(nb == 0 & na > 0 & type == "TCGA-LAML"),             # loh
              adult_total = sum(type == "TCGA-LAML")) %>% 
    rename(pos = start) %>%
    setDT()

methods = c("del", "loss_1", "gain_1", "copy_neutral", "copy_unchanged", "gain", "loh")

for(method in methods) {
  print(method)
  kid_count = paste("kid", method, sep="_")
  adult_count = paste("adult", method, sep="_")
  p.col = paste("p", method, sep="_")
  q.col = paste("q", method, sep="_")
  lfdr.col = paste("lfdr", method, sep="_")
  probe.summary[[p.col]] = apply(probe.summary[, c(kid_count, adult_count, "kid_total", "adult_total"), with=F], 1, function(x) {
                  ct = c(x[1], x[2], x[3]-x[1], x[4] - x[2]); 
                  dim(ct) = c(2, 2); 
                  fisher.test(ct, conf.int = F, alternative = "two.sided")$p.value})
  probe.summary[[p.col]][which(probe.summary[[p.col]] > 1)] = 1
  qobj = qvalue(probe.summary[[p.col]])
  probe.summary[[q.col]] = qobj$qvalues
  probe.summary[[lfdr.col]] = qobj$lfdr
}
fwrite(probe.summary, "probe.summary.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# > summary(probe.summary$q_loss_1)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.003837 1.000000 1.000000 0.999600 1.000000 1.000000
# > summary(probe.summary$q_loh)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.007988 1.000000 1.000000 0.980600 1.000000 1.000000

loh = probe.summary[q_loh <= 0.1]                                                                                                                                                        
fwrite(loh, "potential.loh.probe.summary.tsv", col.names=T, row.names=F, sep="\t", quote=F)
loh.summary = data.table(chr=c("chr2", "chr3", "chr3", "chr15"), 
                         start=c(144422326, 15740317, 152195459, 85404341), 
                         end = c(144558269, 15827904, 152356070, 85519400), 
                         nprobe = c(53, 116, 15, 5), 
                         range = c(135944, 87588, 160612, 115060), 
                         p.min = c(9.252e-07, 3.806e-07, 1.400e-06, 1.307e-05), 
                         q.min = c(1.307e-05, 0.007988, 0.01257, 0.09867))
loh.summary$previous.probe.pos = NA
loh.summary$next.probe.pos = NA
for(i in 1:nrow(loh.summary)) {
  pos = probe[chr == loh.summary$chr[1]]$start
  loh.summary$previous.probe.pos[i] = max(pos[which(pos < loh.summary$start[i])])
  loh.summary$next.probe.pos[i] = min(pos[which(pos > loh.summary$end[i])])
}
loh.summary$previous.gap = loh.summary$start - loh.summary$previous.probe.pos
loh.summary$next.gap = loh.summary$next.probe.pos - loh.summary$end
fwrite(loh.summary, "potential.loh.region.summary.tsv", col.names=T, row.names=F, sep="\t", quote=F)

gencode = fread("gencode.v22.genes.txt")
names(gencode)[1] = "chr"
setkey(loh.summary, chr, start, end)
loh.summary$region.index = 1:nrow(loh.summary)
ov = foverlaps(gencode, loh.summary, mult="all", nomatch=0L)
fwrite(ov, "potential.loh.gene.tsv", col.names=T, row.names=F, sep="\t", quote=F)








