# GetCommonSegment accept input of a segmentation data frame, and output all non-overlapping sub-segments
# At default setting, the function will use "chr" column as the key, and use "start" amd "end" column as inclusive segment boundaries
# These can be changed by passing parameter of keys, start and end. In particular, a typical use is to use an array of c("sample", "chr") as keys
# Parameters
#  - dt: a data frame contains keys, start and end
#  - keys: aggregation key column(s). Default: chr
#  - start: start coordinate
#  - end: end coordinate
#  - reduce: reduce output segment numbers to those have been used; or otherwise, all possible segments. Default: TRUE
#  - min: min coordinate. Default: 1
#  - max: max coordinate. Dafault: 1e9. Please pick a number that is larger than the longest contig length. 1e9 is good for human genome
# Outout: a data.table object that have columns of keys, start, end (all have been setkey'ed) 
GetCommonSegment <- function(dt, keys = "chr", start = "start", end = "end", reduce = T, min = 1, max = 1e9) {
  # clean up data, and rename columns
  library(data.table)

  dt <- dt %>% 
    select(keys, start, end) %>% 
    rename(start = start, 
           end = end) %>%
    distinct() %>%
    setDT()

  dt.keys <- dt %>% 
    select(keys) %>% 
    distinct()
  
  # find all boundaries: orginal start and end; start - 1 as new end; end + 1 and new start; min as new start and max as new end
  bound1 <- melt(dt, id.var=keys, measure.var=c("start", "end"), variable.name="type", value.name="position") %>%
    distinct()    
  bound2 <- bound1 %>% 
    mutate(position = ifelse(type == "start", position - 1, position + 1), 
           type = ifelse(type == "start", "end", "start")) %>% 
    filter(position >= min & position <= max)  
  bound3 <- cbind(dt.keys, type = "start", position = min)
  bound4 <- cbind(dt.keys, type = "end", position = max)
  bound <- rbind(bound1, bound2, bound3, bound4) %>% 
    distinct() %>% 
    mutate(type = factor(type, levels = c("start", "end"))) %>%
    arrange_(.dots=c(keys, "position", "type")) %>%
    setDT()

  # convert long format into segment
  w <- seq(1, nrow(bound), 2)
  seg <- cbind(bound[w, ..keys], start = bound$position[w], end = bound$position[-w]) %>%
    setDT() %>%
    setkeyv(c(keys, "start", "end"))

  # validate segment and remove those are not used
  if(reduce) {
    setkeyv(dt, c(keys, "start", "end"))
    ov <- foverlaps(seg, dt, mult="all", nomatch=0L, which=T)
    seg <- seg[unique(ov$xid)] %>%
      setkeyv(c(keys, "start", "end"))
  }

  # return seg
  return(seg)
}


library(data.table)
library(dplyr)

# s1 is the ascat output based on SNP6; 
# s2 is the SNP6 based segmentation; 
# s3 is the WGS based segmentation

# read and reformat
meta <- fread("GDC.SNP6.pairs.11678.02082019.tsv")
s1 <- fread("gunzip -c TCGA.TARGET.SNP6.ASCAT.11582.1_off_modified.segment.txt.gz")
names(s1) <- c("aliquot_id", "chr", "start", "end", "major", "minor", "cn")
s1$aliquot_barcode <- meta$tumor_barcode[match(s1$aliquot_id, meta$tumor_aliquot)]
s1$portion <- substr(s1$aliquot_barcode, 1, 20)

s2 <- fread("~/Downloads/wgs/ucec.seg133.v2.gistic.seg.txt")
names(s2) <- c("aliquot_barcode", "chr", "start", "end", "nprobe", "segmean")
s2$portion <- substr(s2$aliquot_barcode, 1, 20)
s2$chr = paste0("chr", s2$chr)

s3 <- fread("~/Downloads/wgs/ucec.wgs133.cnvkit.seg.txt")
names(s3)[1] = "chr"
names(s3)[5] = "segmean"

# subset s1, s2 and s3 to have common aliquot, and only keep autosomes
common.portion <- intersect(intersect(s1$portion, s2$portion), s3$portion)
s1 <- s1[portion %in% common.portion & chr %in% paste0("chr", 1:22)]
s2 <- s2[portion %in% common.portion & chr %in% paste0("chr", 1:22)]
s3 <- s3[portion %in% common.portion & chr %in% paste0("chr", 1:22)]

# compare DNACopy vs ASCAT to figure out cutoffs
t1 <- s1 %>%
  select(-c(aliquot_id, major, minor, aliquot_barcode)) %>%
  mutate(ascat_index = row_number()) %>%
  setDT() %>%
  setkey(portion, chr, start, end)
t2 <- s2 %>%
  select(-aliquot_barcode) %>% 
  mutate(dnacopy_index = row_number()) %>%
  setDT() %>%
  setkey(portion, chr, start, end)

ov <- foverlaps(t1, t2, type="any", mult="all", nomatch=0L) 

# exclude multiple mappings
dup.dnacopy_index = ov$dnacopy_index[which(duplicated(ov$dnacopy_index))]
ov2 <- ov[!dnacopy_index %in% dup.dnacopy_index]
ov2$cna <- "del"
ov2$cna[ov2$cn>=1] <- "loss"
ov2$cna[ov2$cn>=2] <- "neutral"
ov2$cna[ov2$cn>=3] <- "gain"
ov2$cna[ov2$cn>=4] <- "amp"
ov2$cna = factor(ov2$cna, levels = c("del", "loss", "neutral", "gain", "amp"))


ov <- foverlaps(s1, s2, type="any", mult="all", nomatch=0L) 



barcode <- intersect(meta$tumor_barcode, unique(s2$aliquot))
meta <- meta[tumor_barcode %in% barcode]
s1 <- s1[aliquot %in% meta$tumor_aliquot]
s1$aliquot <- meta$tumor_barcode[match(s1$aliquot, meta$tumor_aliquot)]
s1$portion <- substr(s1$aliquot, 1, 16)


s2 <- s2[aliquot %in% barcode]
s3 <- s3[aliquot %in% barcode]





# formating segments
f1 <- s1 %>% 
  mutate(portion = substr(aliquot, 1, 20)) %>%
  rename(segmean.v1 = segmean) %>% 
  select(-c(aliquot, nprobe)) %>%
  setDT() %>%
  setkey(portion, chr, start, end)

f2 <- s2 %>% 
  mutate(portion = substr(aliquot, 1, 20)) %>%
  rename(segmean.v2 = segmean) %>%
  select(-c(aliquot, nprobe)) %>%  
  setDT() %>%
  setkey(portion, chr, start, end)

f3 <- s3 %>% 
  mutate(portion = substr(aliquot, 1, 20)) %>%
  rename(segmean.v3 = segmean) %>%
  select(-c(aliquot, nprobe)) %>%
  setDT() %>%
  setkey(portion, chr, start, end)


dt <- rbind(f1[, -"segmean.v1"], f2[, -"segmean.v2"], f3[, -"segmean.v3"])
keys <- c("portion", "chr")
seg <- GetCommonSegment(dt, keys)


# find boundaries
ov1 <- foverlaps(f1, seg, mult="all", nomatch=NA)
ov2 <- foverlaps(f2, seg, mult="all", nomatch=NA)
ov3 <- foverlaps(f3, seg, mult="all", nomatch=NA)
ov <- merge(merge(ov1, ov2, by=c("portion", "chr", "start", "end"), all=T), ov3, by=c("portion", "chr", "start", "end"), all=T) %>%
  select(portion, chr, start, end, segmean.v1, segmean.v2, segmean.v3) %>%
  setDT()
fwrite(ov, "ucec.segment.block.tsv", col.names=T, row.names=F, sep="\t", quote=F)

library(wCorr)
ov$width <- ov$end - ov$start + 1


ov.cor <- ov[!is.na(segmean.v1) & !is.na(segmean.v2) & !is.na(segmean.v3)] %>% 
  group_by(portion) %>% 
  summarise(pearson.v12 = weightedCorr(x=segmean.v1, y=segmean.v2, weights = width, method="Pearson"), 
            pearson.v13 = weightedCorr(x=segmean.v1, y=segmean.v3, weights = width, method="Pearson"), 
            pearson.v23 = weightedCorr(x=segmean.v2, y=segmean.v3, weights = width, method="Pearson"), 
            spearman.v12 = weightedCorr(x=segmean.v1, y=segmean.v2, weights = width, method="Spearman"), 
            spearman.v13 = weightedCorr(x=segmean.v1, y=segmean.v3, weights = width, method="Spearman"), 
            spearman.v23 = weightedCorr(x=segmean.v2, y=segmean.v3, weights = width, method="Spearman")) %>% 
  setDT()

ov12.cor <- ov[!is.na(segmean.v1) & !is.na(segmean.v2)] %>% 
  group_by(portion) %>% 
  summarise(pair.pearson.v12 = weightedCorr(x=segmean.v1, y=segmean.v2, weights = width, method="Pearson"), 
            pair.spearman.v12 = weightedCorr(x=segmean.v1, y=segmean.v2, weights = width, method="Spearman")) %>%
  setDT()

ov13.cor <-  ov[!is.na(segmean.v1) & !is.na(segmean.v3)] %>% 
  group_by(portion) %>% 
  summarise(pair.pearson.v13 = weightedCorr(x=segmean.v1, y=segmean.v3, weights = width, method="Pearson"), 
            pair.spearman.v13 = weightedCorr(x=segmean.v1, y=segmean.v3, weights = width, method="Spearman")) %>%
  setDT()

ov23.cor <-  ov[!is.na(segmean.v2) & !is.na(segmean.v3)] %>% 
  group_by(portion) %>% 
  summarise(pair.pearson.v23 = weightedCorr(x=segmean.v2, y=segmean.v3, weights = width, method="Pearson"), 
            pair.spearman.v23 = weightedCorr(x=segmean.v2, y=segmean.v3, weights = width, method="Spearman")) %>%
  setDT()

corr <-  merge(ov.cor, merge(merge(ov12.cor, ov13.cor, by="portion"), ov23.cor, by="portion"), by="portion")
fwrite(corr, "ucec.segment.correlations.tsv", col.names=T, row.names=F, sep="\t", quote=F)

ov.complete <- ov[!is.na(segmean.v1) & !is.na(segmean.v2) & !is.na(segmean.v3)]
with(ov.complete, weightedCorr(x=segmean.v1, y=segmean.v2, weights = width, method="Pearson")) 
0.9991632
with(ov.complete, weightedCorr(x=segmean.v1, y=segmean.v3, weights = width, method="Pearson")) 
0.8476169
with(ov.complete, weightedCorr(x=segmean.v2, y=segmean.v3, weights = width, method="Pearson")) 
0.8479212
with(ov.complete, weightedCorr(x=segmean.v1, y=segmean.v2, weights = width, method="Spearman")) 
0.9987158
with(ov.complete, weightedCorr(x=segmean.v1, y=segmean.v3, weights = width, method="Spearman")) 
0.5849407
with(ov.complete, weightedCorr(x=segmean.v2, y=segmean.v3, weights = width, method="Spearman")) 
0.5853945

with(ov[!is.na(segmean.v1) & is.na(segmean.v2) & !is.na(segmean.v3)], weightedCorr(x=segmean.v1, y=segmean.v3, weights = width, method="Pearson")) 
0.6539076
with(ov[is.na(segmean.v1) & !is.na(segmean.v2) & !is.na(segmean.v3)], weightedCorr(x=segmean.v2, y=segmean.v3, weights = width, method="Pearson")) 
0.7540452
with(ov[!is.na(segmean.v1) & is.na(segmean.v2) & !is.na(segmean.v3)], weightedCorr(x=segmean.v1, y=segmean.v3, weights = width, method="Spearman")) 
0.7481753
with(ov[is.na(segmean.v1) & !is.na(segmean.v2) & !is.na(segmean.v3)], weightedCorr(x=segmean.v2, y=segmean.v3, weights = width, method="Spearman")) 
0.8956037


library(matrixStats)
wmean1 <- weightedMean(ov.complete$segmean.v1, ov.complete$width)
wmean2 <- weightedMean(ov.complete$segmean.v2, ov.complete$width)
wmean3 <- weightedMean(ov.complete$segmean.v3, ov.complete$width)
wsd1 <- weightedSd(ov.complete$segmean.v1, ov.complete$width)
wsd2 <- weightedSd(ov.complete$segmean.v2, ov.complete$width)
wsd3 <- weightedSd(ov.complete$segmean.v3, ov.complete$width)
ovz <- ov.complete %>% 
  mutate(z1 = (segmean.v1 - wmean1) / wsd1, 
         z2 = (segmean.v2 - wmean2) / wsd2, 
         z3 = (segmean.v3 - wmean3) / wsd3) %>% 
  setDT()

> with(ovz, sum(abs(z1-z2)*width)/sum(width))
[1] 0.002912635
> with(ovz, sum(abs(z1-z3)*width)/sum(width))
[1] 0.3106342
> with(ovz, sum(abs(z2-z3)*width)/sum(width))
[1] 0.3103565

> sum(ov[!is.na(segmean.v1)]$width)/133
[1] 2737152494
> sum(ov[!is.na(segmean.v2)]$width)/133
[1] 2737593380
> sum(ov[!is.na(segmean.v3)]$width)/133
[1] 2758395772

sum(ov[!is.na(segmean.v1) & is.na(segmean.v2) & !is.na(segmean.v3)]$width) / 133
350946.8
sum(ov[is.na(segmean.v1) & !is.na(segmean.v2) & !is.na(segmean.v3)]$width) / 133
788369.3
sum(ov[!is.na(segmean.v1) & is.na(segmean.v2)]$width) / 133
358807.9
sum(ov[is.na(segmean.v1) & !is.na(segmean.v2)]$width) / 133
799694.2


sum(ov$width) / 133
2803029176  87.34121
sum(ov[!is.na(segmean.v1) & !is.na(segmean.v2) & !is.na(segmean.v3)]$width) / 133
2692179468  83.88718
sum(ov[!is.na(segmean.v1) & !is.na(segmean.v2) & is.na(segmean.v3)]$width) / 133
44614218  1.39016
sum(ov[!is.na(segmean.v1) & is.na(segmean.v2) & !is.na(segmean.v3)]$width) / 133
350946.8  0.01093535
sum(ov[is.na(segmean.v1) & !is.na(segmean.v2) & !is.na(segmean.v3)]$width) / 133
788369.3  0.02456525
sum(ov[!is.na(segmean.v1) & is.na(segmean.v2) & is.na(segmean.v3)]$width) / 133
7861.083  0.000244948
sum(ov[is.na(segmean.v1) & !is.na(segmean.v2) & is.na(segmean.v3)]$width) / 133
11324.9 0.0003528791
sum(ov[is.na(segmean.v1) & is.na(segmean.v2) & !is.na(segmean.v3)]$width) / 133
65076988  2.027771

total = 3209286105

wmean1 <- weightedMean(ov$segmean.v1, ov$width, na.rm=T)
wmean2 <- weightedMean(ov$segmean.v2, ov$width, na.rm=T)
wmean3 <- weightedMean(ov$segmean.v3, ov$width, na.rm=T)
wsd1 <- weightedSd(ov$segmean.v1, ov$width, na.rm=T)
wsd2 <- weightedSd(ov$segmean.v2, ov$width, na.rm=T)
wsd3 <- weightedSd(ov$segmean.v3, ov$width, na.rm=T)
ovz <- ov %>% 
  mutate(z1 = (segmean.v1 - wmean1) / wsd1, 
         z2 = (segmean.v2 - wmean2) / wsd2, 
         z3 = (segmean.v3 - wmean3) / wsd3) %>% 
  setDT()

## work on GISTIC peaks
peak1 <- fread("~/SCRATCH/wgs/v1_output/regions_track.conf_99.bed")
names(peak1) <- c("chr", "start", "end", "type")
peak1$type <- ifelse(grepl("Any-AP", peak1$type), "amp", "del")
setkey(peak1, chr, start, end)

peak2 = fread("~/SCRATCH/wgs/v2_output/regions_track.conf_99.bed")
names(peak2) = c("chr", "start", "end", "type")
peak2$type <- ifelse(grepl("Any-AP", peak2$type), "amp", "del")
setkey(peak2, chr, start, end)

peak3 = fread("~/SCRATCH/wgs/v3_output/regions_track.conf_99.bed")
names(peak3) = c("chr", "start", "end", "wgs.type")
peak3$wgs.type <- ifelse(grepl("Any-AP", peak3$wgs.type), "amp", "del")
peak3$width <- peak3$end - peak3$start + 1
setkey(peak3, chr, start, end)

seg <- GetCommonSegment(rbind(peak1[, -"type1"], peak2[, -"type2"], peak3[, -"type3"]))

ov1 <- foverlaps(peak1, seg, mult="all", nomatch=NA)
ov2 <- foverlaps(peak2, seg, mult="all", nomatch=NA)
ov3 <- foverlaps(peak3, seg, mult="all", nomatch=NA)
ov <- merge(merge(ov1, ov2, by=c("chr", "start", "end"), all=T), ov3, by=c("chr", "start", "end"), all=T) %>%
  select(chr, start, end, type1, type2, type3) %>%
  mutate(width = end - start + 1) %>%
  setDT()

ov.complete <- ov[!is.na(type1) & !is.na(type2) & !is.na(type3)]
ov.del <- ov[type3 == "del"]
ov.amp <- ov[type3 == "amp"]


amp <- peak3[wgs.type=="amp"]
amp$index = 1:nrow(amp)

amp1 <- foverlaps(amp, peak1[type=="amp"], nomatch=NA)

del1 <- unlist(fread("~/SCRATCH/wgs/v1_output/del_genes.conf_99.txt", h=F)[1])
del1 <- del1[-which(is.na(del1) | del1 == "cytoband")]
del2 <- unlist(fread("~/SCRATCH/wgs/v2_output/del_genes.conf_99.txt", h=F)[1])
del2 <- del2[-which(is.na(del2) | del2 == "cytoband")]
del3 <- unlist(fread("~/SCRATCH/wgs/v3_output/del_genes.conf_99.txt", h=F)[1])
del3 <- del3[-which(is.na(del3) | del3 == "cytoband")]

> length(del1)
[1] 41
> length(del2)
[1] 37
> length(del3)
[1] 73
> length(intersect(del1,del3))
[1] 7
> length(intersect(del2,del3))
[1] 7
> length(intersect(del1,del2))
[1] 36


amp1 <- unlist(fread("~/SCRATCH/wgs/v1_output/amp_genes.conf_99.txt", h=F)[1])
amp1 <- amp1[-which(is.na(amp1) | amp1 == "cytoband")]
amp2 <- unlist(fread("~/SCRATCH/wgs/v2_output/amp_genes.conf_99.txt", h=F)[1])
amp2 <- amp2[-which(is.na(amp2) | amp2 == "cytoband")]
amp3 <- unlist(fread("~/SCRATCH/wgs/v3_output/amp_genes.conf_99.txt", h=F)[1])
amp3 <- amp3[-which(is.na(amp3) | amp3 == "cytoband")]

length(amp1)
length(amp2)
length(amp3)
length(intersect(amp1,amp3))
length(intersect(amp2,amp3))
length(intersect(amp1,amp2))


fwrite(ov, "ucec.segment.block.tsv", col.names=T, row.names=F, sep="\t", quote=F)







# find boundaries
ov1 <- foverlaps(f1, seg, mult="all", nomatch=NA)
ov2 <- foverlaps(f2, seg, mult="all", nomatch=NA)
ov3 <- foverlaps(f3, seg, mult="all", nomatch=NA)
ov <- merge(merge(ov1, ov2, by=c("portion", "chr", "start", "end"), all=T), ov3, by=c("portion", "chr", "start", "end"), all=T) %>%
  select(portion, chr, start, end, segmean.v1, segmean.v2, segmean.v3) %>%
  setDT()



w_pearson(d[, 1:2], d$width)
> w_pearson(d[, 2:4], weight=d$width, use="complete.obs")
           segmean.v1 segmean.v2 segmean.v3
segmean.v1  1.0000000  0.9676108  0.8491583
segmean.v2  0.9676108  1.0000000  0.8381360
segmean.v3  0.8491583  0.8381360  1.0000000


