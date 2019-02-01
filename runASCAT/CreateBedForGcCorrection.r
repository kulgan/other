library(data.table)
library(dplyr)

setwd("~/SCRATCH/snp6cnv/gc")

snp <- fread("SNP6.hg38.remap.position.txt")
setnames(snp, c("probeid", "chr", "pos", "type"))
nsnp <- nrow(snp)

width <- c(12,25,50,100,250,500,1000,2500,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000)
window <- windows * 2
window[1] <- 25
nwindow <- length(windows)

sizes <- fread("hg38.chrom.sizes", h=F)
setnames(sizes, c("chr", "len"))

snp <- merge(snp[, 1:3, with=F], sizes)
snp <- snp[rep(1:nsnp, each = nwindow)]
snp$width <- rep(width, nsnp)
snp$window <- rep(window, nsnp)

snp <- snp %>% 
  mutate(start = pmax(pos - width, 0), 
         end = pmin(pos + width, len)) %>%
  select(chr, start, end, probeid, pos, window) %>%
  setDT()

snp$start <- bit64::as.integer64(snp$start)
snp$end <- bit64::as.integer64(snp$end)
snp$pos <- bit64::as.integer64(snp$pos)
snp$window <- bit64::as.integer64(snp$window)

fwrite(snp, "snp.window.bed", col.names=F, row.names=F, sep="\t", quote=F)
