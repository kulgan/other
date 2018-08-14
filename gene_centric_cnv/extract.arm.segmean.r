library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)

setwd("~/Downloads/segs")
theme_set(theme_gray(base_size = 18))

# read centromere loc
centromere = fread("centromere.hg38.txt")

# read All TCGA segments
segs = fread("gunzip -c All.TCGA.nocnv_grch38.seg.v2.txt.gz")
segs = segs %>% mutate(chrn = as.integer(ifelse(Chromosome == "X", 23, Chromosome)))
segs = data.table(merge(segs, centromere, by = "chrn"))
segs = segs[order(Sample, chrn, Start)]

# read All TCGA GISTIC2 sample Cutoffs
cutoff = fread("tcga.gistic2.sample.cutoff.txt")
p.low.boxplot = ggplot(cutoff, aes(Project, Low, fill = Project)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + guides(fill = F)
p.high.boxplot = ggplot(cutoff, aes(Project, High, fill = Project)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + guides(fill = F)
grid.arrange(p.high.boxplot, p.low.boxplot, ncol=1)


# get arm-level segmean
arm.summary = data.table(segs %>% 
	mutate(arm = ifelse(Start > centromere, "q", "p")) %>% 
	group_by(Sample, Chromosome, arm) %>%
	summarise(arm.segmean = sum(Num_Probes * Segment_Mean) / sum(Num_Probes)))

# get sample-level segmean
sample.summary = data.table(segs %>% 
	group_by(Sample) %>%
	summarise(sample.segmean = sum(Num_Probes * Segment_Mean) / sum(Num_Probes)))

# get min and max of arm-level segmean per sample
sample.extreme = data.table(arm.summary %>%
	group_by(Sample) %>%
	summarise(low = min(arm.segmean), high = max(arm.segmean)))

# get gap distance and plot basic stats of segmentation data
temp1 = segs[-nrow(segs), ]
names(temp1) = paste0("left.", names(temp1))
temp2 = segs[-1, ]
names(temp2) = paste0("right.", names(temp2))
temp = data.table(cbind(temp1, temp2))
temp$Gap_Length = temp$right.Start - temp$left.End

p.gaplen.dist = ggplot(temp[left.Sample == right.Sample & left.Chromosome == right.Chromosome], aes(Gap_Length)) + geom_density(fill="gray")
p.segmean.dist = ggplot(segs, aes(Segment_Mean)) + geom_density(fill = "gray")
p.seglen.dist = ggplot(segs, aes(End - Start + 1)) + geom_density(fill = "gray") + xlab("Segment_Length") 
p.nprobe.dist = ggplot(segs, aes(Num_Probes)) + geom_density(fill = "gray") 
grid.arrange(p.segmean.dist, p.nprobe.dist, p.seglen.dist + scale_x_log10(), p.gaplen.dist + scale_x_log10(), ncol=2)








qt = c(0, 0.001, 0.01, seq(0.1, 0.9, 0.1), 0.99, 0.999, 1)

quantile(arm.summary$arm.segmean, qt)
          0%         0.1%           1%          10%          20%          30%
-2.135803153 -0.953841225 -0.674253946 -0.257472095 -0.103923316 -0.034647253
         40%          50%          60%          70%          80%          90%
-0.003665163  0.006924811  0.025966201  0.060178677  0.120615347  0.251448819
         99%        99.9%         100%
 0.661432769  1.073425451  2.032341005



sample.summary = merge(sample.summary, sample.extreme)

quantile(sample.summary$sample.segmean, qt)
           0%          0.1%            1%           10%           20%
-0.1813374129 -0.0661761146 -0.0075550582  0.0007780342  0.0011865124
          30%           40%           50%           60%           70%
 0.0014688068  0.0017562429  0.0020948814  0.0026625850  0.0045979959
          80%           90%           99%         99.9%          100%
 0.0104102833  0.0167474709  0.0237832464  0.0285977329  0.0398945292

quantile(sample.summary$low, qt)
         0%        0.1%          1%         10%         20%         30%
-2.13580315 -1.51265753 -1.12835426 -0.78015072 -0.66811709 -0.57802570
        40%         50%         60%         70%         80%         90%
-0.49950000 -0.42727468 -0.35390000 -0.27580000 -0.17700559 -0.05146432
        99%       99.9%        100%
-0.00130000  0.01245413  0.01600000

quantile(sample.summary$high, qt)
          0%         0.1%           1%          10%          20%          30%
-0.071500000 -0.006847520  0.005409957  0.044526583  0.192500000  0.315315120
         40%          50%          60%          70%          80%          90%
 0.397917045  0.468825353  0.538676565  0.617727551  0.718534529  0.865796472
         99%        99.9%         100%
 1.280286274  1.578685138  2.032341005

quantile(segs$Num_Probes, qt)
       0%      0.1%        1%       10%       20%       30%       40%       50%
     1.00      2.00      2.00      5.00     20.00     81.00    273.00    753.00
      60%       70%       80%       90%       99%     99.9%      100%
  1628.00   3122.00   5982.00  13660.00  64341.45 107215.00 132317.00

quantile(segs$End - segs$Start, qt)
       0%      0.1%        1%       10%       20%       30%       40%       50%
        0         7        72      4896     28132    137926    494658   1369694
      60%       70%       80%       90%       99%     99.9%      100%
  2988375   5719921  10886098  24954042 124926658 199632188 244349219


