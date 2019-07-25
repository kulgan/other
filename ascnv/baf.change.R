library(data.table)
library(dplyr)
library(reshape2)

setwd("~/Downloads/baf")


pair.file = "./TCGA.TARGET.SNP6.ASCAT.11582.tsv"
baf.file = "./test1k"

pair = fread(pair.file)
# select unique tumor for each patient
# The following order is not for general use. It's just designed to break known ties
keys = c("1-10", "2-10", "6-10", "5-10", "4-10", "9-10", "1-11", "2-11", "6-11", "3-14", "4-14", "9-14", "1-12", "1-14", "3-10", "3-11")
pair = pair %>%
  filter(selected == T) %>%
  mutate(key = paste(tumor_sample_code, normal_sample_code, sep="-"), 
         key = factor(key, levels = keys)) %>% 
  arrange(patient, key, tumor_barcode) %>% 
  select(-key) %>% 
  setDT()
pair = pair[!duplicated(patient)]
# find column index of pairs
pair = pair %>%
   mutate(normal_index = match(normal_cel, colnames(baf)), 
          tumor_index = match(tumor_cel, colnames(baf))) %>% 
   setDT()
fwrite(pair, "GDC.SNP6.pairs.10879.02192019.tsv", col.names=T, row.names=F, sep="\t", quote=F)


data = fread("test1k", h=F)
setnames(data, names(header))
baf = as.matrix(data[, 2:ncol(data), with=F])
rownames(baf) = data$probeid

tumor.baf = baf[, pair$tumor_index]
normal.baf = baf[, pair$normal_index]
baf.diff = tumor.baf - normal.baf

normal.type = matrix("other", nrow(normal.baf), ncol(normal.baf))
normal.type[which(normal.baf > 0.9)] = "BB"
normal.type[which(normal.baf < 0.1)] = "AA"
normal.type[which(normal.baf > 0.4 & normal.baf < 0.6)] = "AB"



baf = as.matrix(data[, 2:ncol(data), with=F])
rownames(baf) = data$probeid


# looking for duplicated aliquot from the same sample
cel.file <- "./GDC.CEL.file.20190216.tsv"
cel <- fread(cel.file)
cel2 <- cel %>%
  filter(dnu == F) %>% 
  mutate(analyte = ifelse(grepl("^TCGA", barcode), substr(barcode, 1, 20), barcode)) 
cel3 <- merge(cel2, cel2, by=c("program", "project", "patient", "sample", "sample_code", "analyte")) %>%
  filter(aliquot.x != aliquot.y) %>%
  select(program, project, patient, sample, sample_code, analyte, aliquot.x, barcode.x, cel.x, aliquot.y, barcode.y, cel.y) %>%
  distinct(sample, .keep_all= TRUE) %>%
  mutate(sample_type = ifelse(sample_code < 10, "tumor", "normal")) %>%
  setDT()



baf1 = baf[, match(cel3$cel.x, colnames(baf))]
colnames(baf1) = cel3$analyte
baf1 = melt(baf1)
names(baf1) = c("probeid", "analyte", "baf1")

baf2 = baf[, match(cel3$cel.y, colnames(baf))]
colnames(baf2) = cel3$analyte
baf2 = melt(baf2)
names(baf2) = c("probeid", "analyte", "baf2")

baf1$baf2 = baf2$baf
baf2 = baf1
baf1$diff = baf1$baf2 - baf1$baf1
baf2$diff = baf2$baf1 - baf2$baf2
baf1$BAF = baf1$baf1
baf2$BAF = baf2$baf2
d = rbind(baf1, baf2) %>% 
  select(probeid, analyte, BAF, diff) 
d$tissue = "normal"
d$tissue[which(d$analyte %in% cel3[sample_type=="tumor"]$analyte)] = "tumor"
d = data.table(d)

l = lm(diff ~ BAF + tissue, data = baf.combined)

baf.diff = baf2 - baf1
baf.mean = (baf2 + baf1)/2
probe.summary = data.table(t(apply(baf.diff, 1, summary)))
probe.summary$sd = apply(baf.diff, 1, sd)



PairedTPermTest <- function(x, y, nperm) {
  n = length(x)
  p = t.test(x, y, paired=T)$p.value
  ps 
}




pair.file = "./GDC.SNP6.pairs.11678.02082019.tsv"


pair.file = "./TCGA.TARGET.SNP6.ASCAT.11582.tsv"

header = fread("tcga_snp6_baf_header", h=T)


baf = as.matrix(data[, 2:ncol(data), with=F])
rownames(baf) = data$probeid

pair = fread(pair.file)
tumor.baf = baf[, match(pair$tumor_cel, colnames(baf))]
normal.baf = baf[, match(pair$normal_cel, colnames(baf))]
baf.diff = tumor.baf - normal.baf

normal.type = matrix("other", nrow(normal.baf), ncol(normal.baf))
normal.type[which(normal.baf > 0.9)] = "BB"
normal.type[which(normal.baf < 0.1)] = "AA"
normal.type[which(normal.baf > 0.4 & normal.baf < 0.6)] = "AB"

> table(normal.type)/length(normal.type)
normal.type
       AA        AB        BB     other 
0.3352574 0.2246357 0.3118961 0.1282108 

p = apply(baf.diff, 1, function(x) t.test(x, alternative = "two.sided")$p.value)

baf.diff2 = baf.diff
baf.diff2[which(normal.type != "AB")] = NA
p = apply(baf.diff2, 1, function(x) t.test(x, alternative = "two.sided")$p.value)
p2 = sort(p)

cel = fread("../GDC.CEL.file.tsv")
dnu.list = fread("../SNP6.DNU.txt", h=F)$V1

cel = cel %>% 
  mutate(portion = ifelse(program == "TCGA", substr(barcode, 1, 20), barcode), 
         dnu = barcode %in% dnu.list, 
         tissue = ifelse(sample_code >= 10 & sample_code <= 20, "normal", "tumor")) %>%
  setDT() 

pair = fread("../GDC.SNP6.pairs.tsv")
snp = fread("zcat /home/ubuntu/SCRATCH/snp6cnv/snp6.na35.remap.hg38.txt.gz")

d = fread("xaa", h=T)


##############################
# Portion Duplicates
##############################

# re-order of cel row to make it the same as sample in d 
sample.list = gsub(".B Allele Freq", "", names(d))
m = match(sample.list, cel$filename)
# cel = cel[m]
cel$index = 1:nrow(cel)

# generate list of portion duplicates for future control analysis
cel2 = cel[dnu == F][order(barcode)]
temp = cel2[which(duplicated(cel2$portion))]
dup = data.table(portion = temp$portion, barcode2 = temp$barcode, index2 = temp$index)
temp = cel2[which(!duplicated(cel2$portion))]
m = match(dup$portion, temp$portion)
dup$barcode1 = temp$barcode[m]
dup$index1 = temp$index[m]
dup$program = temp$program[m]
dup$sample_code = temp$sample_code[m]

##############################
# Good pair
##############################
# generate list of good tumor-only paired samples for analysis
temp = pair %>%
  mutate(dnu = tumor_barcode %in% dnu.list | normal_barcode %in% dnu.list) 


FindReplacement = function(x, cel){
  w = which(grepl(substr(x, 1, 16), cel$barcode) & cel$dnu == F)
  if(length(w) == 0) {
    w = which(grepl(substr(x, 1, 15), cel$barcode) & cel$dnu == F)
    if(length(w) == 0) {
      tissue_type = cel[barcode == x]$tissue[1]
      w = cel[tissue ==  tissue_type & patient == substr(x, 1, 12) & dnu == F]$index
    }
  }
# return(length(w))  
  return(paste(cel$barcode[w], collapse=" "))
}


pair1 = pair[!(normal_barcode %in% dnu.list | tumor_barcode %in% dnu.list)]
pair2 = pair[normal_barcode %in% dnu.list | tumor_barcode %in% dnu.list]

pair3 = pair2
pair3$normal_dnu = pair3$normal_barcode %in% dnu.list
pair3$tumor_dnu = pair3$tumor_barcode %in% dnu.list

l = c(pair3$normal_barcode, pair3$tumor_barcode)
l = l[which(l %in% dnu.list)]
d = data.table(origin = l)
d$n = sapply(d$origin, function(x) FindReplacement(x, cel))
d$new = sapply(d$origin, function(x) FindReplacement(x, cel))
