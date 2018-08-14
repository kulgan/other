library(data.table)
library(dplyr)

library(ggplot2)
library(gridExtra)


gain.threshold = 0.2
loss.threshold = -0.2
amp.threshold = 0.8
deletion.threshold = -0.6


setwd("~/Downloads/segs")
theme_set(theme_gray(base_size = 18))


# read All TCGA segments
segs = fread("gunzip -c All.TCGA.nocnv_grch38.seg.v2.txt.gz")
samples = unique(segs$Sample)
nsample = length(samples)

# sort segs by absolute value of Segmean
segs = segs[order(-abs(segs$Segment_Mean))]

# read gene locations
gene = fread("gencode.v22.genes.txt")
ngene = nrow(gene)
gene[, chr := gsub("chr", "", seqname)]
gene = gene[, c("chr", "start", "end", "gene_id")]
setkey(gene, chr, start, end)

# count how many genes overlaps with multiple segments and how many with none
for(i in 1:nsample) {
	print(i)
	sample = samples[i]
	# get overlaps between segment and genes
	ov = foverlaps(segs[Sample == sample], gene, by.x = c("Chromosome", "Start", "End"), type="any", mult="all", nomatch=0L)
	# find extreme value of segmean, as Segment_Mean in segs were pre-sorted by absolute values
	ov.extreme = ov[!duplicated(ov$gene_id)]
	# find those genes that overlaps with more than one segment with segments mean exceeding both middle gain/loss threshold
	ov.multi = ov %>% 
				filter(gene_id %in% unique(ov[duplicated(ov$gene_id)]$gene_id)) %>% 
				group_by(gene_id) %>% summarise(count = length(gene_id), loss = min(Segment_Mean) <= loss.threshold, gain = max(Segment_Mean >= gain.threshold)) %>% 
				filter(count > 1 & loss == T & gain == T)
	# output result into data matrix
	extreme[, i] = ov.extreme[match(rownames(extreme), gene_id)]$Segment_Mean
	multi[which(rownames(multi) %in% ov.multi$gene_id), i] = T
}



# find overlaps between segs and gene loci. If multiple overlaps exist for 
# one gene, only keep the segs with the largest absolute segmean value

# matrix extreme keep values of the extreme segmean per sample per gene
extreme = matrix(NA, ngene, nsample)
rownames(extreme) = gene$gene_id
colnames(extreme) = samples
# matrix multi has value T only one the indicated gene in that sample overlaps with more than one segments exceeding both middle gain/loss threshold
multi = matrix(F, ngene, nsample)
rownames(multi) = gene$gene_id
colnames(multi) = samples

for(i in 1:nsample) {
	print(i)
	sample = samples[i]
	# get overlaps between segment and genes
	ov = foverlaps(segs[Sample == sample], gene, by.x = c("Chromosome", "Start", "End"), type="any", mult="all", nomatch=0L)
	# find extreme value of segmean, as Segment_Mean in segs were pre-sorted by absolute values
	ov.extreme = ov[!duplicated(ov$gene_id)]
	# find those genes that overlaps with more than one segment with segments mean exceeding both middle gain/loss threshold
	ov.multi = ov %>% 
				filter(gene_id %in% unique(ov[duplicated(ov$gene_id)]$gene_id)) %>% 
				group_by(gene_id) %>% summarise(count = length(gene_id), loss = min(Segment_Mean) <= loss.threshold, gain = max(Segment_Mean >= gain.threshold)) %>% 
				filter(count > 1 & loss == T & gain == T)
	# output result into data matrix
	extreme[, i] = ov.extreme[match(rownames(extreme), gene_id)]$Segment_Mean
	multi[which(rownames(multi) %in% ov.multi$gene_id), i] = T
}
#
# 0.018% has multi == T
(-100,-0.6] (-0.6,-0.2]  (-0.2,0.2]   (0.2,0.8]   (0.8,100]
0.022715697 0.106428718 0.743164519 0.120043752 0.007647314

(-100,-0.6] (-0.6,-0.2]  (-0.2,0.2]   (0.2,0.7]   (0.7,100]
 0.02271570  0.10642872  0.74316452  0.11593651  0.01175455

(-100,-0.6] (-0.6,-0.2]  (-0.2,0.2]   (0.2,0.6]   (0.6,100]
 0.02271570  0.10642872  0.74316452  0.10897146  0.01871961

(-100,-0.6] (-0.6,-0.3]  (-0.3,0.3]   (0.3,0.6]   (0.6,100]
 0.02271570  0.06796771  0.82848431  0.06211268  0.01871961


# filter result above to remove chrY amd chrM
gene = gene[chr %in% c(1:22, "X")]
ngene = nrow(gene)
w = which(rownames(extreme) %in% gene$gene_id)
extreme = extreme[w, ]
multi = multi[w, ]

# add centromere and determine arms
centromere = fread("centromere.hg38.txt")
centromere$chrn[which(centromere$chrn == 23)] = "X"
gene = data.table(merge(gene, centromere, by.x = "chr", by.y = "chrn") %>% 
		mutate(arm = ifelse(start > centromere, "q", "p")))



gene.sortbyend = data.table(merge(gene, centromere, by.x = "chr", by.y = "chrn") %>% 
		mutate(arm = ifelse(start > centromere, "q", "p")) %>% 
		arrange(chr, end))

for(i in 1:nsample) {
	print(i)
	d = merge(gene, data.table(gene_id = rownames(extreme), value = extreme[, i]))
	for(j in 1:ngene) {
		if(is.na(d$value[j])) {
			my.gene = d[j]
			arm.extreme = d[chr == d$chr[j] & arm == d$arm[j] & !is.na(value)]
			right.gene = arm.extreme[start > my.gene$end][order(start)][1]
			left.gene = arm.extreme[end < my.gene$start][order(-end)][1]
			if(is.na(left.gene$gene_id)) {
				d$value = right.gene$value
			} elseif (is.na(right.gene$gene_id)) {
				d$value = left.gene$value
			} else {
				d$value = ((my.gene$start - left.gene$end) * right.gene$value +
							(right.gene$start - my.gene$end) * left.gene$value ) /
							(right.gene$start - left.gene$end)
			}
		}
	}

}

temp = foverlaps(segs, gene, by.x = c("Chromosome", "Start", "End"), type="any")

my.segs = segs[Sample == "154e9edd-b354-4e62-be0f-8b26ebddffbc"]
temp = foverlaps(my.segs, gene, by.x = c("Chromosome", "Start", "End"), type="any", which=T, mult="all", nomatch=0L)








# read centromere loc





g
data.table(merge(gene, centromere, by.x = "chr", by.y = "chrn"))
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


