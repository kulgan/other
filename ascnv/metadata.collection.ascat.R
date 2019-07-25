library(data.table)
library(dplyr)

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

fwrite(pair1, "pair.keep.11648.02082019.txt", col.names=T, row.names=F, sep="\t", quote=F)
fwrite(pair2, "pair.dnu.162.02082019.txt", col.names=T, row.names=F, sep="\t", quote=F)
fwrite(temp, "pair.dnu_replaced.30.02082019.txt", col.names=T, row.names=F, sep="\t", quote=F)
fwrite(d, "GDC.SNP6.pairs.11678.02082019.tsv", col.names=T, row.names=F, sep="\t", quote=F)












