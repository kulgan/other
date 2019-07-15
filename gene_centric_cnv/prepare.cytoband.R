# The script read cytoBand.hg38.txt.gz file downloaded from UCSC genome browser
# Output 2 files
# 1. chromosome.info.human.hg38.txt that contains information of each major chromosome with columns chr/ centromere/ chr_len 
# 2. cytoband.info.human.hg38.txt that contains information of each cytoband with columns chr/ start/ end/ cytoband/ gstain/ arm 

library(data.table)
library(dplyr)

options(scipen = 999, digits=4)

# input cytoband file from UCSC
cytoband.file <- "gunzip -c cytoBand.hg38.txt.gz"

# define major contigs and arms
contigs <- c(1:22, "X", "Y")
arms <- c(t(outer(contigs, c("p", "q"), FUN = paste0)))

# read cytoband file and set colnames
cyto <- fread(cytoband.file, h=F )
setnames(cyto, c("chr", "start", "end", "cytoband", "gstain"))

# filter for major contigs, remove 'chr' in contig names, add arm information and sort
# One additional base is added to all start positions to convert bed coordinates to ranges
cyto <- cyto %>% 
  mutate(chr = gsub('chr', '', chr)) %>%
  filter(chr %in% contigs) %>% 
  mutate(chr = factor(chr, levels = contigs), 
      arm = paste0(chr, substr(cytoband, 1, 1)), 
      start = start + 1) %>% 
  arrange(chr, start)

# extract centromere information by using the end location of acen band on p arms
centromere.info <- cyto %>% 
  filter(gstain == 'acen' & substr(cytoband, 1, 1) == "p") %>% 
  mutate(centromere = end) 

# extract chromosome length, and merge with centromere location
chr.info <- cyto %>% 
  group_by(chr) %>%
  summarise(chr_len = max(end)) %>%
  mutate(centromere = centromere.info$centromere) 

# get information of p arms
p.arm.info <- chr.info %>% 
  mutate(start = 1, 
         end = centromere, 
         arm = paste0(chr, "p")) 

# get information of q arms
q.arm.info <- chr.info %>% 
  mutate(start = centromere + 1, 
         end = chr_len, 
         arm = paste0(chr, "q"))

# combine information of p arms and q arms, and order by arm
arm.info <- rbind(p.arm.info, q.arm.info) %>%
    mutate(arm = factor(arm, levels = arms)) %>%
    arrange(arm) %>%
    select(chr, arm, start, end)

# output to disk
write.table(cyto, 'cytoband.info.human.hg38.txt', col.names=T, row.names=F, sep='\t', quote=F)
write.table(chr.info, 'chromosome.info.human.hg38.txt', col.names=T, row.names=F, sep='\t', quote=F)
write.table(arm.info, 'arm.info.human.hg38.txt', col.names=T, row.names=F, sep='\t', quote=F)

