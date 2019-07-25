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
meta <- fread("SKCM.collection_type.txt")
chrs <- paste0("chr", 1:22)

# subsetting segment
seg <- fread("zcat TCGA.TARGET.SNP6.ASCAT.11582.1_off_modified.segment.txt.gz")
names(seg) <- c("aliquot", "chr", "start", "end", "na", "nb", "cn")
seg <- seg %>%
  filter(aliquot %in% meta$tumor_aliquot & chr %in% chrs) %>% 
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
ov <- foverlaps(seg, probe, mult="all", type="any", nomatch=0L) 
names(ov)[6:7] <- c("seg_start", "seg_end")
# this should be move before overlapping
m <- match(ov$aliquot, meta$tumor_aliquot)
ov$collection_type = meta$collection_type[m]

# collect counts for gain and loss
cn.summary <- ov %>% 
    group_by(probeid) %>%
    summarise(primary_loss = sum(cn < 2 & collection_type == "primary_only"), 
              primary_neutral = sum(cn == 2 & collection_type == "primary_only"), 
              primary_gain = sum(cn > 2 & collection_type == "primary_only"), 
              metastasis_loss = sum(cn < 2 & collection_type == "metastatic_only"), 
              metastasis_neutral = sum(cn == 2 & collection_type == "metastatic_only"), 
              metastasis_gain = sum(cn > 2 & collection_type == "metastatic_only")) %>% 
    mutate(primary_count = primary_loss + primary_neutral + primary_gain, 
           metastatis_count = metastasis_loss + metastasis_neutral + metastasis_gain) %>% 
    setDT()

# Fisher Exact tests on gain and losses
cn.summary$p_value_loss <- apply(cn.summary[, 2:7, with=F], 1, function(x) {
                            ct = c(x[1], x[2]+x[3], x[4], x[5]+x[6]); 
                            dim(ct) = c(2,2);
                            fisher.test(ct, conf.int = F, alternative = "two.sided")$p.value})

cn.summary$p_value_gain <- apply(cn.summary[, 2:7, with=F], 1, function(x) {
                            ct = c(x[1] + x[2], x[3], x[4] + x[5], x[6]); 
                            dim(ct) = c(2,2);
                            fisher.test(ct, conf.int = F, alternative = "two.sided")$p.value})

# re-annotate each probe with position information
probe.summary <- probe %>% 
  mutate(pos = start, 
         chr = factor(chr, levels = chrs)) %>%
  select(probeid, chr, pos) %>% 
  arrange(chr, pos) %>% 
  merge(cn.summary, by = "probeid") %>%
  setDT()

# calculate q values
library(qvalue)
qobj.gain <- qvalue(probe.summary$p_value_gain)
qobj.loss <- qvalue(probe.summary$p_value_loss)
probe.summary$qvalue.gain <- qobj.gain$qvalues
probe.summary$lfdr.gain <- qobj.gain$lfdr
probe.summary$qvalue.loss <- qobj.loss$qvalues
probe.summary$lfdr.loss <- qobj.gain$lfdr
fwrite(probe.summary, "probe.summary.SKCM.FisherExact.tsv", col.names=T, row.names=F, sep="\t", quote=F)

