library(data.table)
library(dplyr)
library(ggplot2)

# get centromere positions
cyto = fread("~/SCRATCH/gistic/cyto.tsv")
centromere = cyto %>% 
	filter(chrn %in% 1:23) %>%
	mutate(arm = gsub("([0-9X]+[pq]).+", "\\1", name)) %>% 
	filter(grepl("q", arm)) %>% 
	group_by(chrn) %>% 
	summarize(centromere = min(start)) 
write.table(centromere, "centromere.hg38.txt", col.names=T, row.names=F, sep="\t", quote=F)

# extract sample data
sample = fread("segment.meta.txt")
sample = sample %>% 
	mutate(sample_type = substr(aliquot_barcode, 14, 15)) %>% 
	select(aliquot_id, aliquot_barcode, sample_type)

# read segs
segs = fread("TCGA.BRCA.nocnv_grch38.seg.v2.txt")

# add chrn
segs$chrn = as.integer(ifelse(segs$Chromosome=="X", 23, segs$Chromosome))
segs = merge(segs, centromere, by="chrn")

# add sample_type information
segs = merge(segs, sample, by.x = "Sample", by.y = "aliquot_id")

# add arm level, chromosome level and sample level segmeans
segs = segs %>% 
	mutate(	pq = ifelse(Start > centromere, "q", "p"), 
			arm = paste0(chrn, pq), 
			sample.arm = paste(Sample, arm, sep="_"),
			sample.chr = paste(Sample, chrn, sep="_")) %>% 
	group_by(sample.arm) %>%
	mutate(arm_seg_sum = sum(Segment_Mean * Num_Probes), arm_probe_num = sum(Num_Probes)) %>%
	mutate(arm_segmean = arm_seg_sum/arm_probe_num, arm_segmean_exclude = (arm_seg_sum - (Segment_Mean * Num_Probes)) / (arm_probe_num - Num_Probes)) %>% 
	group_by(sample.chr) %>%
	mutate(chr_seg_sum = sum(Segment_Mean * Num_Probes), chr_probe_num = sum(Num_Probes)) %>%
	mutate(chr_segmean = chr_seg_sum/chr_probe_num, chr_segmean_exclude = (chr_seg_sum - arm_seg_sum) / (chr_probe_num - arm_probe_num)) %>% 
	group_by(Sample) %>%
	mutate(sample_seg_sum = sum(Segment_Mean * Num_Probes), sample_probe_num = sum(Num_Probes)) %>%
	mutate(sample_segmean = sample_seg_sum/sample_probe_num, sample_segmean_exclude = (sample_seg_sum - chr_seg_sum) / (sample_probe_num - chr_probe_num))  

# Check importance of arm, chr and sample level segmean. 
summary(lm(Segment_Mean ~ arm_segmean, data=segs))$adj.r.squared
> 0.1025107
summary(lm(Segment_Mean ~ arm_segmean_exclude, data=segs))$adj.r.squared
> 0.07997801

summary(lm(Segment_Mean ~ chr_segmean, data=segs))$adj.r.squared
> 0.06827442
summary(lm(Segment_Mean ~ chr_segmean_exclude, data=segs))$adj.r.squared
> 0.001423884

summary(lm(Segment_Mean ~ sample_segmean, data=segs))$adj.r.squared
> 0.06526786
> 0.002850883


# add left segment
left  = segs %>% select(Sample, Start, End, Num_Probes, Segment_Mean, sample.arm) %>% mutate(density = (End - Start + 1) / Num_Probes)
left = rbind(left[nrow(left), ], left[-nrow(left), ])
segs$left.segmean = left$Segment_Mean
segs$left.dist = segs$Start - left$End
segs$left.density = left$density
w = which(left$sample.arm != segs$sample.arm)
segs$left.segmean[w] = NA
segs$left.dist[w] = NA
segs$left.density[w] = NA
left = rbind(left[nrow(left), ], left[-nrow(left), ])
segs$left2.segmean = left$Segment_Mean
segs$left2.dist = segs$Start - left$End
segs$left2.density = left$density
w = which(left$sample.arm != segs$sample.arm)
segs$left2.segmean[w] = NA
segs$left2.dist[w] = NA
segs$left2.density[w] = NA


# add right segment
right  = segs %>% select(Sample, Start, End, Num_Probes, Segment_Mean, sample.arm) %>% mutate(density = (End - Start + 1) / Num_Probes)
right = rbind(right[-1, ], right[1, ])
segs$right.segmean = right$Segment_Mean
segs$right.dist = right$Start - segs$End
segs$right.density = right$density
w = which(right$sample.arm != segs$sample.arm)
segs$right.segmean[w] = NA
segs$right.dist[w] = NA
segs$right.density[w] = NA
right = rbind(right[-1, ], right[1, ])
segs$right2.segmean = right$Segment_Mean
segs$right2.dist = right$Start - segs$End
segs$right2.density = right$density
w = which(right$sample.arm != segs$sample.arm)
segs$right2.segmean[w] = NA
segs$right2.dist[w] = NA
segs$right2.density[w] = NA

# add closest segment
segs$closest.segmean = segs$right.segmean
segs$closest.dist = segs$right.dist
w = with(segs, which(right.dist > left.dist | is.na(right.dist)))
segs$closest.segmean[w] = segs$left.segmean[w]
segs$closest.dist[w] = segs$left.dist[w]

segs$closest2.segmean = segs$right2.segmean
segs$closest2.dist = segs$right2.dist
w = with(segs, which(right2.dist > left2.dist | is.na(right2.dist)))
segs$closest2.segmean[w] = segs$left2.segmean[w]
segs$closest2.dist[w] = segs$left2.dist[w]







segs = data.table(segs)

summary(lm(Segment_Mean ~ right.segmean, data=segs))$adj.r.squared
0.05209542
summary(lm(Segment_Mean ~ right2.segmean, data=segs))$adj.r.squared
0.3961172


summary(lm(Segment_Mean ~ left.segmean, data=segs))$adj.r.squared
0.05209542
summary(lm(Segment_Mean ~ left2.segmean, data=segs))$adj.r.squared
0.396117


summary(lm(Segment_Mean ~ left.segmean + right.segmean, data=segs))$adj.r.squared
0.06764362
summary(lm(Segment_Mean ~ left2.segmean + right2.segmean, data=segs))$adj.r.squared
0.5323152


summary(lm(Segment_Mean ~ closest.segmean, data=segs))$adj.r.squared
0.05563075
summary(lm(Segment_Mean ~ closest2.segmean, data=segs))$adj.r.squared
0.3964362

> summary(segs$left.dist)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      0     437    1463    8243    3573 4139869   41560 
> summary(segs$left2.dist)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
        3     68516    629235   3752954   3492558 136095756     74338 
> summary(segs$closest.dist)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      0     219     637    3225    1838 4139869    8782 
> summary(segs$closest2.dist)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
        3     24103    286398   2930694   2053633 136095756     22107 

segs$factor2 = cut(segs$closest2.dist, breaks = c(0, 100, 200, 500, 1000, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9))


> summary(lm(Segment_Mean ~ 0 + closest2.segmean, data=segs))$adj.r.squared
[1] 0.3962876

> segs[!is.na(closest2.segmean)] %>% group_by(factor2) %>% summarise(summary(lm(Segment_Mean ~ 0 + closest2.segmean))$adj.r.squared)
# A tibble: 10 x 2
   factor2       `summary(lm(Segment_Mean ~ 0 + closest2.segmean))$adj.r.squared`
   <fct>                                                                    <dbl>
 1 (0,100]                                                                  0.979
 2 (100,200]                                                                0.993
 3 (200,500]                                                                0.738
 4 (500,1e+03]                                                              0.867
 5 (1e+03,1e+04]                                                            0.369
 6 (1e+04,1e+05]                                                            0.480
 7 (1e+05,1e+06]                                                            0.509
 8 (1e+06,1e+07]                                                            0.315
 9 (1e+07,1e+08]                                                            0.137
10 (1e+08,1e+09]                                                            0.200

> summary(lm(Segment_Mean ~ 0 + closest2.segmean + arm_segmean_exclude, data=segs))$adj.r.squared
[1] 0.407111

> segs[!is.na(closest2.segmean)] %>% group_by(factor2) %>% summarise(summary(lm(Segment_Mean ~ 0 + closest2.segmean + arm_segmean_exclude))$adj.r.squared)
# A tibble: 10 x 2
   factor2       `summary(lm(Segment_Mean ~ 0 + closest2.segmean + arm_segmean_exclude))$adj.r.squared`
   <fct>                                                                                          <dbl>
 1 (0,100]                                                                                        0.979
 2 (100,200]                                                                                      0.993
 3 (200,500]                                                                                      0.751
 4 (500,1e+03]                                                                                    0.869
 5 (1e+03,1e+04]                                                                                  0.419
 6 (1e+04,1e+05]                                                                                  0.487
 7 (1e+05,1e+06]                                                                                  0.513
 8 (1e+06,1e+07]                                                                                  0.327
 9 (1e+07,1e+08]                                                                                  0.171
10 (1e+08,1e+09]                                                                                  0.223


segs = segs %>% mutate(linear.segmean = (right2.dist * left2.segmean + left2.dist * right2.segmean) / (right2.dist + left2.dist))
w = which(is.na(segs$linear.segmean))
segs$linear.segmean[w] = segs$closest2.segmean[w]
segs = data.table(segs)

> summary(lm(Segment_Mean ~ 0 + linear.segmean, data=segs))$adj.r.squared
[1] 0.432708

> segs[!is.na(linear.segmean)] %>% group_by(factor2) %>% summarise(r.adj = summary(lm(Segment_Mean ~ 0 + linear.segmean))$adj.r.squared, beta = summary(lm(Segment_Mean ~ 0 + linear.segmean))$coefficients[1,1]) 
# A tibble: 10 x 2
   factor2       `summary(lm(Segment_Mean ~ 0 + linear.segmean))$adj.r.squared`
   <fct>                                                                  <dbl>
 1 (0,100]                                                                0.979
 2 (100,200]                                                              0.993
 3 (200,500]                                                              0.749
 4 (500,1e+03]                                                            0.873
 5 (1e+03,1e+04]                                                          0.409
 6 (1e+04,1e+05]                                                          0.515
 7 (1e+05,1e+06]                                                          0.553
 8 (1e+06,1e+07]                                                          0.352
 9 (1e+07,1e+08]                                                          0.145
10 (1e+08,1e+09]                                                          0.200









