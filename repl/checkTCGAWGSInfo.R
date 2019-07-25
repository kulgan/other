library(data.table)
library(dplyr)

setwd("~/Downloads/target/")

rg = fread("tcga.wgs.rg.tsv", na.strings = "")
names(rg)[1] = "rg.index"

sr = fread("tcga.wgs.sr.tsv", na.strings = "")
names(sr)[1] = "sr.index"

ar = fread("tcga.wgs.ar.tsv", na.strings = "")
names(ar)[1] = "ar.index"


##########################################
# AR Check
##########################################
# test if only one alignement wf per AR
sum(ar$ar.num.align_wf != 1) == 0 

##########################################
# SR Check
##########################################
# test if all aligned SR are SURs 
table(sr[!is.na(align_wf)]$sr.data_format)
# FASTQ 
#  1453

# Look for alignment wf without AR
sr[!align_wf %in% ar$align_wf & !is.na(align_wf)]

# Look for SR not linked to RG
sum(! sr$sr %in% rg$sr)
# 0 

# Look for duplicated alignment wf on SR
sum(duplicated(sr[!is.na(align_wf)]$align_wf))
# 0

##########################################
# RG Check
##########################################
# Look for RG name == 0
sum(rg$rg_name==0)
# 1336
# read_length could also be 0 

# combine SR with AR
m = match(sr$align_wf, ar$align_wf)
sr$ar = ar$ar[m]
sr$has_ar = !is.na(sr$ar)

# combine RG with SR
m = match(rg$sr, sr$sr)
rg$sr.batch_id = sr$sr.batch_id[m]
rg$sr.file_name = sr$sr.file_name[m]
rg$has_ar = sr$has_ar[m]

# 
tar.rg = rg[sr.format == "tar"]$rg
fq.rg = rg[sr.format == "fq"]$rg
bam.rg = rg[sr.format == "bam"]$rg
rg = rg %>% 
  mutate(has_tar = rg %in% tar.rg, 
         has_fq = rg %in% fq.rg, 
         has_bam = rg %in% bam.rg) %>% 
  setDT()


aliquot = rg %>% 
  group_by(aliquot) %>%
  summarise(has_tar = "tar" %in% sr.format, 
            has_fq = "fq" %in% sr.format, 
            has_bam = "bam" %in% sr.format) %>% 
  setDT()
> xtabs(~ has_fq + has_bam + has_tar, aliquot)
, , has_tar = FALSE
       has_bam
has_fq  FALSE TRUE
  FALSE     0    2
  TRUE    187   84

, , has_tar = TRUE
       has_bam
has_fq  FALSE TRUE
  FALSE     0    1
  TRUE      0  354
# suggest to remove all BAMs

# assume all bams were removed
rg = rg[sr.format != "bam"]

aliquot = rg %>% 
  group_by(aliquot) %>%
  summarise(has_tar = "tar" %in% sr.format, 
            has_fq = "fq" %in% sr.format) %>% 
  setDT()
xtabs(~ has_fq + has_tar, aliquot)
       has_tar
has_fq  FALSE TRUE
  FALSE     0    1
  TRUE    271  354
# which aliquot is that? 
aliquot[has_tar & !has_fq]$aliquot
# "b37e7042-0d58-5fa7-9d08-24f8b8c62e46"
rg[aliquot=="b37e7042-0d58-5fa7-9d08-24f8b8c62e46"]$sr
# "69b872cd-8934-4cdc-8bee-33f13fb76c29" this one should be split

xtabs(~ sr.format + sr.batch_id, data=rg)
         sr.batch_id
sr.format    1    4    5
      fq     0 1170  166
      tar  356    0    0

# look into fastq sizes when both fq and tar exist
size = data.table(aliquot = aliquot[has_tar & has_fq]$aliquot, 
                  tar.size = 0, 
                  fq.size = 0)
for(i in 1:nrow(size)) {
  a = size$aliquot[i]
  size$tar.size[i] = sr[sr == unique(rg[aliquot == a & sr.format == "tar"]$sr)]$sr.file_size
  size$fq.size[i] = sum(sr[sr %in% rg[aliquot == a & sr.format == "fq"]$sr]$sr.file_size)
}
> summary(size$fq.size / size$tar.size)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  1.000   1.000   1.000   1.001   1.000   1.197 


library(stringr)

# only need fastqs
rg = rg[sr.format == "fq"]
rg = rg %>% 
  mutate(direction = str_sub(gsub(".fastq.gz", "", sr.file_name), -1, -1), 
         corename = gsub("^(.+)(_(r||R|X|x)(1|2))\\.fastq.gz$", "\\1", sr.file_name)) %>% 
  setDT() 

sum(duplicated(rg$sr.file_name))
# 0

s = rg %>% 
  group_by(aliquot, corename) %>%
  summarise(num_fq = length(sr)) %>% 
  setDT()

table(s$num_fq)
#  2 
# 668 


w = which(duplicated(s$aliquot))
a = s$aliquot[w]



