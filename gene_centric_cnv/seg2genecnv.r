#!/usr/bin/env Rscript

# SEG2GENECNV
# author: Zhenyu Zhang
# date: 08/15/2018
# require: R > 3.20 
#	   	data.table
#		futile.logger
#		optparse
#		tools


# ReadGenomicTSV read a TSV file, and return a data.table of file content
# - if file has a header, it must include a column named Chromosome
# - if file does not have a header, default colnames of Chromosome/Start/End will be added to the first 3 columns
# the function will also remove "chr" string from the contig names, and change 23 to X, 24 to Y
ReadGenomicTSV = function(file, header = T, segment = T, keep.chry = F) {
	# permissible values
	contigs = c(1:23, "X", "Y")
	chr.contigs = paste0("chr", contigs)

	# read file
	tsv = fread(file, h = header)


	# add default header for the first 3 columns if header == F
	if (!header & segment) {
		colnames(tsv)[1:3] = c("Chromosome", "Start", "End")
	}

	# filter out contigs not in chr1-22/X
	tsv = tsv[Chromosome %in% c(contigs, chr.contigs)]

	# reformat contig names without "chr"
	if(sum(tsv$Chromosome %in% contigs) == 0) {
		tsv$Chromosome = gsub("chr", "", tsv$Chromosome)
	}

	# update 23 to X, and 24 to Y
	tsv$Chromosome[which(tsv$Chromosome == "23")] = "X"
	tsv$Chromosome[which(tsv$Chromosome == "24")] = "Y"

	# remove Y chromosome is needed
	if(! keep.chry) {
		tsv = tsv[Chromosome != "Y"]
	}

	# sort 
	tsv$Chromosome = factor(tsv$Chromosome, levels = contigs)
	if(segment){
		tsv = tsv[order(Chromosome, Start, End)]
	} else {
		tsv = tsv[order(Chromosome)]
	}

	# return resulting data.table
	return(tsv)
}



# GetArm will return arm values of segment in GenomicTSV data.table based on the centromere position provided
GetArm = function(genomictsv, centromere) {
	# get arm information of each segment and return
	arm.info = genomictsv %>% 
		left_join(., centromere, by = "Chromosome") %>% 
		mutate(arm = ifelse(Start > Centromere, paste0(Chromosome,"q"), paste0(Chromosome, "p"))) 
	return(arm.info$arm)	
}


# SplitSegByCentromere will split segment by centromere, and adjust Num_Probes and Star or End values
SplitSegByCentromere = function(segs, centromere) {
	# annotate each segment by centromere position, and whether centromere is within the segment
	segs = 	data.table(segs %>% 
				left_join(., centromere, by = "Chromosome") %>% 
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
		select(Chromosome, i.Start, i.End, Gene, potential_segmean, source) %>%
		rename(End = i.End, Start = i.Start)

	# get arm_segmean as potential_segmean for genes not overlaps with any segments
	gene.info.in_gap = gene.info[!Gene %in% ov$Gene] 
	gene.segmean0 = gene.info.in_gap %>%
		mutate(	potential_segmean = arm_segmean, 
				source = "arm") %>% 
		select(Chromosome, Start, End, Gene, potential_segmean, source)

	# get left neighbor segmean as potential_segmean for genes not overlaps with any segments
	gene.segmean0.left = gene.info.in_gap %>% 
		mutate(Start = Start - neighbor.distance) %>% 
		setDT() %>%
		foverlaps(., segs, type = "any", mult = "last", nomatch = 0) %>% 
		mutate(	potential_segmean = seg_segmean, 
				source = "left_neighbor") %>% 
		select(Chromosome, Start, End, Gene, potential_segmean, source)

	# get right neighbor segmean as potential_segmean for genes not overlaps with any segments		
	gene.segmean0.right = gene.info.in_gap %>% 
		mutate(End = End + neighbor.distance) %>% 
		setDT() %>%
		foverlaps(., segs, type = "any", mult = "first", nomatch = 0) %>% 
		mutate(	potential_segmean = seg_segmean, 
				source = "right_neighbor") %>% 
		select(Chromosome, Start, End, Gene, potential_segmean, source)	

	# get arm_segmean as potential_segmean for genes that only overlaps with a small segment
	count = table(ov$Gene)
	gene.segmean1 = ov %>% 
		filter(Gene %in% names(count)[count == 1] & is_small_segment) %>% 
		mutate(	potential_segmean = arm_segmean, 
				source = "arm") %>% 
		select(Chromosome, Start, End, Gene, potential_segmean, source)

	# merge all potential segmeans. (this is the most important intermediate dataset) 
	gene.segmean = rbind(gene.segmean, gene.segmean0, gene.segmean0.left, gene.segmean0.right, gene.segmean1)	

	return(gene.segmean)
	
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
# filter values are also added to indicate the reliability of the segmeans
AggregateSegmean = function(gene.segmean = gene.segmean, cutoffs = cutoffs) {

	setDT(gene.segmean)
	imputed.genes = setdiff(gene.segmean[source == "arm"]$Gene, gene.segmean[source == "overlap"]$Gene)

	# aggregate, get min and max, find extreme value, and thresholded to integer scores by cutoffs
	score = gene.segmean %>%
		group_by(Gene) %>% 
		summarise(	Chromosome = Chromosome[1], 
					Start = Start[1], 
					End = End[1],
					min_segmean = min(potential_segmean), 
					max_segmean = max(potential_segmean)) %>% 
		mutate(	min_score = GetThresholedValue(min_segmean, cutoffs), 
				max_score = GetThresholedValue(max_segmean, cutoffs), 
				similar = abs(min_score - max_score) <= 1, 
				one_segmean = ifelse(max_segmean > 0 & abs(max_segmean) > abs(min_segmean), max_segmean, min_segmean), 
				one_score = GetThresholedValue(one_segmean, cutoffs), 
				filter = "PASS", 
				filter = ifelse(Gene %in% imputed.genes, AddFilter(filter, "impute"), filter), 
				filter = ifelse(similar, filter, AddFilter(filter, "non_det"))) %>%
		select(Chromosome, Start, End, Gene, one_segmean, one_score, filter)

	return(data.table(score))
}

# add off-target filter to gene.score
FilterByBed = function(gene.score, bed) {
	# normalize bed coordinates with 1-indexed interval coordinates
	bed$Start = bed$Start + 1
	# overlap between gene and bed
	setkey(bed, Chromosome, Start, End)
	ov = foverlaps(gene.score, bed, type = "any", mult = "first", nomatch = NA, which = T) 
	
	# get genes without overlaps, and add filter "off_target"
	off.target.index = which(is.na(ov))
	gene.score$filter[off.target.index] = AddFilter(gene.score$filter[off.target.index], "off_target")

	return(gene.score)
}

# reformat and write  
WriteScores = function(gene.score, gene, filters, out.file) {
	# re-oder genes and change colnames
	gene.score = gene.score[match(gene$Gene, gene.score$Gene)] %>%
		mutate(Chromosome = paste0("chr", Chromosome)) %>%
		rename(	Segment_Mean = one_segmean, 
				Copy_Number_Score = one_score, 
				Filter = filter)

	# write date
	date = format(Sys.Date(), "%Y%m%d")
	cat("##fileDate=", date, "\n", sep ="", file = out.file)

	# write filters	
	for(i in 1:nrow(filters)) {
		cat("##Filter=<ID=", filters$ID[i], ",Description=\"", filters$Description[i], "\">", "\n", sep = "", file = out.file, append = T)
	}

	# write content
	suppressWarnings(write.table(gene.score, file = out.file, append = T, col.names = T, row.names = F, sep="\t", quote=F))
}


library(data.table)
library(dplyr)
library(futile.logger)

options(scipen = 999, digits=4)

setwd("~/Downloads/segs")
# theme_set(theme_gray(base_size = 18))


minArmEvidence = 50
minSegEvidence = 4
minBedSegLength = 50
neighborDistance = 1e6

cutoffs = c(-0.6, -0.2, 0.2, 0.6)

filters = data.table(ID = c("PASS", 
							"impute", 
							"non_det", 
							"off_target"), 
					Description = c("All filters passed", 
							"gene-level value is imputed because gene is not overlapped with any segments", 
							"gene-level value is not deterministic due to aggregation of segments with very different segment means", 
							"gene not overlapped with any targeted capture regions"))

bed.file = "../vcrome_2_1_hg19_capture_targets.hg19.bed"
out.file = "output.tsv"


# read All TCGA segments
segs = ReadGenomicTSV("DUNGS_p_TCGA_b84_115_SNP_N_GenomeWideSNP_6_A07_771588.hg19.seg.txt", header = T, segment = T, keep.chry = F)
# centromere location is extracted from UCSC genome browser http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
centromere = ReadGenomicTSV("centromere.position.hg38.txt", header = T, segment = F, keep.chry = F)
gene = ReadGenomicTSV("gencode.v22.gene.short.txt", header = T, segment = T, keep.chry = F)
bed = ReadGenomicTSV(bed.file, header = F, segment = T, keep.chry = F)


# split segments by centromere 
split.segs = SplitSegByCentromere(segs, centromere)

# we assume bed file and gene file do not have regions that contains centromere 

# calcualte arm-level segmean value
arm.segmean = GetArmSegmean(split.segs, minArmEvidence, centromere, bed, minBedSegLength) 

# modify cutoffs by max and min of arm.segmean
new.cutoffs = GetCutoffs(cutoffs, arm.segmean)

# get all potential gene-level segmean values
gene.segmean = GetGeneSegmean(segs = split.segs, arm.segmean = arm.segmean, gene = gene, neighbor.distance = neighborDistance, min.seg.evidence = minSegEvidence)

# get aggregated gene-level CNV scores
gene.score = AggregateSegmean(gene.segmean = gene.segmean, cutoffs = new.cutoffs)

# add off_target filter
if(!is.null(bed.file)){
	gene.score = FilterByBed(gene.score = gene.score, bed = bed)
}

# format and output
WriteScores(gene.score, gene, filters, out.file)














# sample flog code to be integrated later
if(F) {

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


flog.info("Regular CNV Segnment Completed")

# calculate, format and write nocnv segment
flog.info("Masked CNV Segnment Completed")

}