#!/usr/bin/env Rscript

# SEG2GENECNV
# author: Zhenyu Zhang
# date: 08/15/2018
# require: R > 3.20 
#	   	data.table
#		futile.logger
#		optparse
#		tools


library(data.table)
library(dplyr)
library(futile.logger)

# library(ggplot2)
# library(gridExtra)

setwd("~/Downloads/segs")
# theme_set(theme_gray(base_size = 18))

# read centromere loc

min.arm.evidence = 50
cutoffs = c(-0.6, -0.2, 0.2, 0.6)
min.seg.evidence = 4
min.bed.segment.len = 50
neighbor.distance = 1e6

# ReadGenomicTSV read a TSV file, and return a data.table of file content
# - if file has a header, it must include a column named Chromosome
# - if file does not have a header, default colnames of Chromosome/Start/End will be added to the first 3 columns
# the function will also remove "chr" string from the contig names, and chang 23 to X
ReadGenomicTSV = function(file, header = T) {
	# permissible values
	contigs = c(1:23, "X")
	chr.contigs = paste0("chr", contigs)

	# read file
	tsv = fread(file, h = header)

	# add default header for the first 3 columns if header == F
	if (!header) {
		colnames(tsv)[1:3] = c("Chromosome", "Start", "End")
	}

	# filter out contigs not in chr1-22/X
	tsv = tsv[Chromosome %in% c(contigs, chr.contigs)]

	# reformat contig names without "chr"
	if(sum(tsv$Chromosome %in% contigs) == 0) {
		tsv$Chromosome = gsub("chr", "", tsv$Chromosome)
	}

	# update 23 to X
	tsv$Chromosome[which(tsv$Chromosome == "23")] = "X"

	# return resulting data.table
	return(tsv)
}



# GetArm will return arm values of segment in GenomicTSV data.table based on the centromere position provided
GetArm = function(genomictsv, centromere) {
	# get arm information of each segment and return
	arm.info = genomictsv %>% 
		merge(., centromere, by = "Chromosome") %>% 
		mutate(arm = ifelse(Start > Centromere, paste0(Chromosome,"q"), paste0(Chromosome, "p"))) 
	return(arm.info$arm)	
}


# SplitSegByCentromere will split segment by centromere, and adjust Num_Probes and Star or End values
SplitSegByCentromere = function(segs, centromere) {
	# annotate each segment by centromere position, and whether centromere is within the segment
	segs = 	data.table(segs %>% 
				merge(., centromere, by = "Chromosome") %>% 
				mutate(	Cross_Centromere = Start <= Centromere & End > Centromere))

	# if cross-centromere segments exist, they are split into two segments
	if(sum(segs$Cross_Centromere) > 0) {
		# extract segs that cross centromere
		cross.segs = segs[Cross_Centromere == T] %>%
			rowwise() %>% 
			mutate(Num_P_Probes = min(Num_Probes -1, ceiling((Centromere - Start) / (End - Start) * Num_Probes)))

		# split to get p-arm segments 
		p.segs = cross.segs %>% 
			mutate(Num_Probes = Num_P_Probes, End = Centromere) %>%
			select(-Num_P_Probes) 

		# split to get q-arm segments
		q.segs = cross.segs %>% 
			mutate(Num_Probes = Num_Probes - Num_P_Probes, Start = Centromere + 1) %>%
			select(-Num_P_Probes) 	

		# merge p-arm and q-arm segments back 
		segs = rbind(segs[Cross_Centromere == F], p.segs, q.segs)	
	}

	# drop unused columns and return
	return(segs[, -c("Centromere", "Cross_Centromere")])
} 

# GetArmSegmean returns a data.table of chromosome arms and their weighted-average segmeans
GetArmSegmean = function(segs, min.arm.evidence, centromere, bed = NULL, min.bed.segment.len = 50) {
	# define default contigs
	contigs = c(1:22, "X")

	# get all possible arms
	arms = c(outer(contigs, c("p", "q"), FUN = paste0))

	# calcualte arm-level segment means
	# if number of probes on an arm is less than min.arm.evidence, set segmean = 0
	arm.summary = segs %>% 
		mutate(arm = GetArm(., centromere)) %>% 
		group_by(arm) %>%
		summarise( 	num_probes = sum(Num_Probes), 
					segmean = sum(Num_Probes * Segment_Mean) / sum(Num_Probes)) %>% 
		mutate(segmean = ifelse(num_probes < min.arm.evidence, 0, segmean))

	# get arm.summary of arms that do not have any probes
	arm.summary = rbind(arm.summary, data.table(arm = arms, num_probes = 0, segmean = 0)[!arm %in% arm.summary$arm])

	# if number of bed segment on an arm is less than min.arm.evidence, set segmean = 0
	# here, we assume bed file is in good shape (not cover centromeres)
	if(!is.null(bed)) {
		# annotate bed by arm
		bed.arm.summary = bed %>% 
			mutate(arm = GetArm(., centromere)) %>% 
			filter(End - Start >= min.bed.segment.len) %>%
			group_by(arm) %>%
			summarise(nsegment = length(arm))

		# set nsegment < min.arm.evidence arms to have segmean = 0	
		arm.summary[arm %in% bed.arm.summary$arm[which(bed.arm.summary$nsegment < min.arm.evidence)]]$segmean = 0
	}

	return(arm.summary[, -"num_probes"])
}

# GetCutoffs return a new cutoffs array that has modified values
# the 1st element of cutoffs: low_cutoff = min(segmean_of_each_arm, default_low_cutoff)
# the 4th element of cutoffs: high_cutoff = max(segmean_of_each_arm, default_high_cutoff) 
GetCutoffs = function(cutoffs, arm.segmean) {
	cutoffs[1] = min(arm.segmean$segmean, cutoffs[1])
	cutoffs[4] = max(arm.segmean$segmean, cutoffs[4])
	return(cutoffs)
}

# GetThresholedValue return thresholded CN values by cutoffs
# please notice that these bins are closed on the left and open on the right
GetThresholedValue = function(segmean, cutoffs) {
	# assign -2/-1/0/1/2 to each segmean
	cutoffs = c(-Inf, cutoffs, +Inf)
	thresholded = cut(segmean, cutoffs, right = F, labels = -2:2)
	return(as.integer(as.character(thresholded)))
}

# AddFilter add a new filter to existing filter string
# if existing filter is "PASS", replaced with new filter
# else, append new filter to the end, and separate with semicolon
AddFilter = function(old, new) {
	return(ifelse(old == "PASS", new, paste(old, new, sep=";")))
}

# GetGemeSegmean find all potential segmeans for each gene:
# 1. overlapped segment segmeans
# 2. arm-level segmean if no overlapped segments or the only overlapped segment is mall (less than min.seg.evidence probe support)
# 3. segmean of left and right neighbor segments (within neighbor.distance) if no overlapped segments
# It returns a data.frame of 3 columns Gene, potential_segmean and source.  
# "source" indicates the source of the segmean: overlap, arm, left_neighbor, right_neighbor. 
GetGeneSegmean = function(segs, arm.segmean, gene, neighbor.distance = 1e6, min.seg.evidence = 4) {
	# sort segs, define small_segment, and setkey for later foverlaps 
	segs = segs %>% 
		mutate(is_small_segment = Num_Probes < min.seg.evidence) %>% 
		rename(seg_segmean = Segment_Mean) %>%
		select(Chromosome, Start, End, seg_segmean, is_small_segment) %>%
		arrange(Chromosome, Start) %>% 
		setDT() %>% 
		setkey(., Chromosome, Start, End)


	# get arm segmean of each gene	
	gene.info = gene %>% 
		mutate(arm = GetArm(., centromere)) %>% 
		inner_join(., arm.segmean, by="arm") %>%
		rename(arm_segmean = segmean) %>% 
		select(Chromosome, Start, End, Gene, arm_segmean) %>% 
 		setDT()

	# overlap between gene and segment
	ov = foverlaps(gene.info, segs, type = "any", mult = "all", nomatch = 0) 

	# get seg_segmean as potential_segmean for genes overlapped with segments
	gene.segmean = ov %>% 
		mutate(	potential_segmean = seg_segmean, 
				source = "overlap") %>%
		select(Gene, potential_segmean, source)

	# get arm_segmean as potential_segmean for genes not overlaps with any segments
	gene.info.in_gap = gene.info[!Gene %in% ov$Gene] 
	gene.segmean0 = gene.info.in_gap %>%
		mutate(	potential_segmean = arm_segmean, 
				source = "arm") %>% 
		select(Gene, potential_segmean, source)

	# get left neighbor segmean as potential_segmean for genes not overlaps with any segments
	gene.segmean0.left = gene.info.in_gap %>% 
		mutate(Start = Start - neighbor.distance) %>% 
		setDT() %>%
		foverlaps(., segs, type = "any", mult = "last", nomatch = 0) %>% 
		mutate(	potential_segmean = seg_segmean, 
				source = "left_neighbor") %>% 
		select(Gene, potential_segmean, source)

	# get right neighbor segmean as potential_segmean for genes not overlaps with any segments		
	gene.segmean0.right = gene.info.in_gap %>% 
		mutate(End = End + neighbor.distance) %>% 
		setDT() %>%
		foverlaps(., segs, type = "any", mult = "first", nomatch = 0) %>% 
		mutate(	potential_segmean = seg_segmean, 
				source = "right_neighbor") %>% 
		select(Gene, potential_segmean, source)	

	# get arm_segmean as potential_segmean for genes that only overlaps with a small segment
	count = table(ov$Gene)
	gene.segmean1 = ov %>% 
		filter(Gene %in% names(count)[count == 1] & is_small_segment) %>% 
		mutate(	potential_segmean = arm_segmean, 
				source = "arm") %>% 
		select(Gene, potential_segmean, source)

	# merge all potential segmeans. (this is the most important intermediate dataset) 
	gene.segmean = rbind(gene.segmean, gene.segmean0, gene.segmean0.left, gene.segmean0.right, gene.segmean1)	

	
}


# SelectOneValue accept an array of numeric values, and return only one of them based on method
# In the default settings, it returns the item with the largest absolute value, and favor negative number when tie
# The function is not used for performance reason 
SelectOneValue = function(arr, method = "extreme", favor.positive = F) {
	# sort arr first to favor selection of positive or negative values with the same absolute value
	arr = sort(arr, decreasing = favor.positive)
	return(arr[which.max(abs(arr))])
}
 
# AggregateSegmean select only one segmean if multiple exist for each gene
# The default method is "extreme", and favor negative values
AggregateSegmean = function(gene.segmean = gene.segmean, cutoffs = cutoffs) {

	setDT(gene.segmean)
	imputed.genes = setdiff(gene.segmean[source == "arm"]$Gene, gene.segmean[source == "overlap"]$Gene)

	# aggregate, get min and max, find extreme value, and thresholded to integer scores by cutoffs
	score = gene.segmean %>%
		group_by(Gene) %>% 
		summarise(	min_segmean = min(potential_segmean), 
					max_segmean = max(potential_segmean)) %>% 
		mutate(	min_score = GetThresholedValue(min_segmean, cutoffs), 
				max_score = GetThresholedValue(max_segmean, cutoffs), 
				similar = abs(min_score - max_score) <= 1, 
				one_segmean = ifelse(max_segmean > 0 & abs(max_segmean) > abs(min_segmean), max_segmean, min_segmean), 
				one_score = GetThresholedValue(one_segmean, cutoffs), 
				filter = "PASS", 
				filter= ifelse(Gene %in% imputed.genes, AddFilter(filter, "impute"), filter), 
				filter = ifelse(similar, filter, AddFilter(filter, "nondet"))) %>%
		select(Gene, one_segmean, one_score, filter)

	score$filter = "PASS"

		
		
	return(score)
}


	# threshold arm segmean 
	# arm.segmean$arm_score = GetThresholedValue(arm.segmean$segmean, cutoffs) 
	# threshold segment segmean

		mutate(	seg_score = GetThresholedValue(Segment_Mean, cutoffs),





check start < End
check at least 2 probes


# read All TCGA segments
segs = ReadGenomicTSV("DUNGS_p_TCGA_b84_115_SNP_N_GenomeWideSNP_6_A07_771588.hg19.seg.txt")
centromere = ReadGenomicTSV("centromere.position.hg38.txt")
gene = ReadGenomicTSV("gencode.v22.gene.short.txt")
bed = ReadGenomicTSV("../vcrome_2_1_hg19_capture_targets.hg19.bed", header = F)



# split segments by centromere 
split.segs = SplitSegByCentromere(segs, centromere)

# we assume bed file and gene file do not have regions that contains centromere 

# calcualte arm-level segmean value
arm.segmean = GetArmSegmean(split.segs, min.arm.evidence, centromere, bed, min.bed.segment.len) 

# modify cutoffs by max and min of arm.segmean
cutoffs = GetCutoffs(cutoffs, arm.segmean)

# get all potential gene-level segmean values
gene.segmean = GetGeneSegmean(segs = split.segs, arm.segmean = arm.segmean, gene = gene, neighbor.distance = neighbor.distance, min.seg.evidence = min.seg.evidence)

# get aggregated gene-level CNV scores
gene.score = AggregateSegmean(gene.segmean = gene.segmean, cutoffs = cutoffs)




gene.segmean = AddBedFilter(gene.segmean, bed)












segs = merge(segs, centromere, by = "Chromosome") %>% 
		mutate(arm = ifelse(Start > Centromere, paste0(Chromosome,"q"), paste0(Chromosome, "p"))) 



cutoffs[1] = min(cutoffs[1], arm.summary$segmean) 
cutoffs[4] = max(cutoffs[4], arm.summary$segmean) 


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











read.meta = function(file) {
	# read SNP6 metadata information in rda/rdata/gz/txt/tsv/csv format
	# get file extension
	flog.info("Loading metadata: %s", file)

	# read and format
	options(warn = -1)
	# gunzip to open gz data. zcat has compatibility issue in OSX
	meta = fread(paste("gunzip -c", file))
	options(warn = 0)
	meta$strand = factor(meta$strand)
	meta$type = factor(meta$type)
	meta$chr = factor(meta$chr, levels = c(1:22, "X", "Y"))
	# return
	return(meta)
}

get.segment = function(data, mode = "allcnv", sampleid = "Sample", seed=12345678) {
        # load data that at least contains segmean, chr, pos, return copy number segments
        if (mode == "nocnv") {
        	# remove freqcnv probes and chromosome Y
        	data = data[freqcnv == FALSE & chr != "Y"]
        }
        CNA.object = CNA( genomdat = data$segmean, chrom = data$chr, maploc = as.numeric(data$pos), data.type = "logratio", sampleid = sampleid)
        CNA.smoothed = smooth.CNA(CNA.object)
        set.seed(seed)
        segs = segment(CNA.smoothed, nperm=10000, alpha=0.01)
        return(segs$output)
}

format.segment = function(segment, gender=NA) {
	# format segment output, modify column names, and remove chrY is gender is female	
	names(segment) = c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
	if (gender == "female") {
		segment = segment[segment$Chromosome != "Y", ]
	}
	return(segment)
}

read.tangentcn = function(file) {
	# read tangentcn data
	flog.info("Loading tangent normalized copy number input file: %s", file)
	options(warn = -1)
	data = fread(file, skip=1, col.names=c("probeid", "signal"))
	options(warn = 0)
	return(data)	
}

annotate.tangentcn = function(tangentcn, meta) {
	# merge tangentcn data with meta data, and add segmean
	flog.info("Annotating tangent data with probeset information from metadata")
	data = merge(meta, tangentcn, by="probeid")
	# remove probes with missing information
	data = data[(! is.na(data$strand)) & (! is.na(data$signal)), ]
	# calculate probeset segment mean 
	data$segmean = log2(data$signal) - 1
	return(data)
}

get.opt = function() {
	# define option_list
	option_list = list(
		make_option(c("-f", "--file"), type = "character", default = NULL, 
	              	help = "tangent copy number file name"),
		make_option(c("--out1"), type = "character", default = "allcnv.txt", 
	              	help = "output copy number segmentation file name [default = %default]"), 
		make_option(c("--out2"), type = "character", default = "nocnv.txt", 
	              	help = "output germline masked copy number segmentation file name [default = %default]"), 
	    make_option(c("-m", "--meta"), type = "character", default = "snp6.na35.liftoverhg38.txt.gz", 
	              	help = "gzipped snp6 probeset metadata name [default = %default]"), 
	    make_option(c("--gender"), type = "character", default = "unknown", 
	              	help = "gender of the patient [default = %default]\n\t\t\t\t      [allowed = male, female, unknown]"),
	    make_option(c("-s", "--sample"), type = "character", default = "Sample", 
	              	help = "Sample ID in output files [default = %default]"), 
	    make_option(c("-l", "--log"), type = "character", default = "stdout", 
	              	help = "log file [default = %default]")   
	)
	
	# parse 
	opt_parser = OptionParser(option_list = option_list);
	opt = parse_args(opt_parser);

	# if log output is not "stdout", write to log file
	if(opt$log != "stdout") {
		flog.appender(appender.file(opt$log))
	}

	# stop if input file is not provided
	if (is.null(opt$file)){
	  print_help(opt_parser)
	  flog.error("no input file provided")
	  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
	}

	# change gender to unknown if not male or female
	opt$gender = tolower(opt$gender)
	if (! opt$gender %in% c("male", "female")) {
		opt$gender = "unknown"
	}

	# display all input options
	flog.info("input: \t\t\t%s", opt$file)
	flog.info("gender: \t\t\t%s", opt$gender)
	flog.info("regular cnv output: \t\t%s", opt$out1)
	flog.info("germline masked cnv output: \t%s", opt$out2)
	flog.info("metadata used: \t\t%s", opt$meta)
	flog.info("sample id used: \t\t%s", opt$sample)

	# return
	return(opt)
}


# load libraries
require(data.table)
require(DNAcopy)
require(optparse)
require(futile.logger)
require(tools)

options(scipen=999)

# get options
opt = get.opt()

# if log output is not "stdout", write to log file
if(opt$log != "stdout") {
	flog.appender(appender.file(opt$log))
}

# check input existance
# this script assume input file is valid TCGA tangent copy number file and will not check file content
if (! file.exists(opt$file)) {
	flog.error("input file does not exist")
	stop("Please provide input file\n", call.=FALSE)
}

# check metadata existance and file md5
# this script assume the metadata has md5sum fce277e26e4d65d187e1ea9800628bb9
md5 = md5sum(opt$meta)
if (is.na(md5)) {
	flog.error("metadata does not exist")
	stop("Please provide metadata\n", call.=FALSE)
} else if (md5 != "fce277e26e4d65d187e1ea9800628bb9") {
	flog.warning("metadata md5sum does not match fce277e26e4d65d187e1ea9800628bb9")
}

# read SNP6 metadata
meta = read.meta(opt$meta)
# read TCGA level 2 tangent CNV file
tangentcn = read.tangentcn(opt$file)
# merge level 2 data with probeset metadata
data = annotate.tangentcn(tangentcn = tangentcn, meta = meta)

# exclude PAR if gender is male or unknown
if (opt$gender %in% c("male", "unknown")) {
	data = data[par == FALSE]
}

# calculate, format and write segment
flog.info("Calculating regular CNV segments")
mode = "allcnv"
segs = get.segment(data, mode = mode, sampleid = opt$sample)
formatted.segs = format.segment(segs, gender = opt$gender)
write.table(formatted.segs, opt$out1, col.names=T, row.names=F, sep="\t", quote=F)
flog.info("Regular CNV Segnment Completed")

# calculate, format and write nocnv segment
flog.info("Calculating masked CNV segments")
mode = "nocnv"
nocnv.segs = get.segment(data, mode = mode, sampleid = opt$sample)
formatted.nocnv.segs = format.segment(nocnv.segs, gender = opt$gender)
write.table(formatted.nocnv.segs, opt$out2, col.names=T, row.names=F, sep="\t", quote=F)
flog.info("Masked CNV Segnment Completed")

