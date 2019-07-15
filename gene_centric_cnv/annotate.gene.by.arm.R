# The script read cytoBand.hg38.txt.gz file downloaded from UCSC genome browser
# Output 2 files
# 1. chromosome.info.human.hg38.txt that contains information of each major chromosome with columns chr/ centromere/ chr_len 
# 2. cytoband.info.human.hg38.txt that contains information of each cytoband with columns chr/ start/ end/ cytoband/ gstain/ arm 

library(data.table)
library(dplyr)

file <- 'gunzip -c cytoBand.hg38.txt.gz'
contigs <- c(1:22, 'X', 'Y')

# read file  
cyto <- fread(file, h=F )
setnames(cyto, c('chr', 'start', 'end', 'cytoband', 'gstain'))

# filter for major contigs, remove 'chr' in contig names, add arm information and sort
# 1 additional base is added to all start positions to convert bed coordinates to segment coordinates
cyto <- cyto %>% 
  mutate(chr = gsub('chr', '', chr)) %>%
  filter(chr %in% contigs) %>% 
  mutate(chr = factor(chr, levels = contigs), 
      arm = paste0(chr, substr(cytoband, 1, 1)), 
      start = start + 1) %>% 
  arrange(chr, start)


# extract centromere information by using the end location of acen band on p arms
centromere = cyto %>% 
  filter(gstain == 'acen' & substr(cytoband, 1, 1) == 'p') %>% 
  select(chr, end) %>% 
  rename(centromere = end) 

# extract chromosome length, and merge with centromere location
chr.info <- cyto %>% 
  group_by(chr) %>%
  summarise(chr_len = max(end)) %>% 
  inner_join(centromere, ., by = 'chr') %>% 
  arrange(chr)

# output to disk
write.table(chr.info, 'chromosome.info.human.hg38.txt', col.names=T, row.names=F, sep='\t', quote=F)
write.table(cyto, 'cytoband.info.human.hg38.txt', col.names=T, row.names=F, sep='\t', quote=F)
