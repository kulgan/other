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
library(tidyr)
library(futile.logger)
library(matrixStats)

options(scipen = 999, digits=4)

setwd("~/Downloads/segs")
# theme_set(theme_gray(base_size = 18))

source("~/github/other/gene_centric_cnv/cnv.utils.r")

minArmEvidence <- 100
minSegEvidence <- 4
noise <- 0.1
cap <- 3
markerless <- F
markerDistance <- 2000


# minBedSegLength = 50
# neighborDistance = 1e6
allowedContigs <- c(1:22, "X")


filters <- data.table(ID = c("PASS", 
							"impute", 
							"non_det", 
							"off_target"), 
					Description = c("All filters passed", 
							"gene-level value is imputed because gene is not overlapped with any segments", 
							"gene-level value is not deterministic due to aggregation of segments with very different segment means", 
							"gene not overlapped with any targeted capture regions"))

# bed.file = "../vcrome_2_1_hg19_capture_targets.hg19.bed"
out.file <- "output.tsv"


# read segmentation file
raw.segs <- ReadAndFilterSegmentation(file = "TCGA.BRCA.10.samples", contigs = allowedContigs, h = F)

# in markerless mode, need to add support probes
if(markerless) {
  raw.segs$nprobes <- with(raw.segs, ceiling((end - start + 1) / markerDistance))
}

# read arm information
arm.info <- ReadAndFilterArmInfo(file = "arm.info.human.hg38.txt", contigs = allowedContigs, h = T)

# clean up segmentation
clean.segs <- CleanUpSegmentations(segs = raw.segs, arm.info = arm.info, min.seg.evidence = minSegEvidence)

# calculate copy number from segmean
clean.segs$copy_number <- SegmeanToCopyNumber(segmean = clean.segs$segmean, cap = cap)

# calculate arm and sample level copy numbers 
arm.copy_number <- GetArmCopyNumber(segs = clean.segs, min.arm.evidence = minArmEvidence) 
sample.copy_number <- GetSampleCopyNumber(arm.copy_number = arm.copy_number)

# Integrate noise levels
# calculate sample threshold 
sample.threshold <- GetSampleThreshold(sample.copy_number = sample.copy_number, noise = noise)



# centromere location is extracted from UCSC genome browser http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
# centromere = ReadAndFilterCentromereFile("centromere.position.hg38.txt", contigs = allowedContigs, h = T)
# gene = ReadAndFilterGeneFile("gencode.v22.gene.short.txt", contigs = allowedContigs, h = T)
# marker = ReadGenomicLocation("gunzip -c marker_file.tsv.gz", h=F)
# bed = ReadGenomicSegment(bed.file, header = F)
gene.info <- ReadAndFilterGene(file = "gencode.v22.gene.short.txt", contigs = allowedContigs, header = T)








event.model = CalcualteCopyNumberEventModel(segs)



]


# sample.cn.summary = GetSampleCopyNumberSummary(arm.cn.summary)

gene.seg.overlaps = GetGeneSegmentOverlaps(segs, arm.cn, genes, minSegEvidence)



GetThresholedValue(gene.segment.overlaps, sample.th)








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








brca = readRDS("BRCA.segs.rds")
brca = data.table(brca[, 2:7])
samples = unique(brca$Sample)
centromere = ReadGenomicTSV("~/SCRATCH/github/other/gene_centric_cnv/centromere.position.hg38.txt", header = T, segment = F, keep.chry = F)
gene = ReadGenomicTSV("~/SCRATCH/github/other/gene_centric_cnv/gencode.v22.gene.short.txt", header = T, segment = T, keep.chry = F)

out = matrix(0, nrow(gene), length(samples))
rownames(out) = gene$Gene
colnames(out) = samples

cutoffs = c(-0.6, -0.1, 0.1, 0.6)

for(i in 1: length(samples)) {
	print(i)
	sample = samples[i]
	segs = brca[Sample == sample]
	split.segs = SplitSegByCentromere(segs, centromere)
	arm.segmean = GetArmSegmean(split.segs, minArmEvidence, centromere, bed = NULL, minBedSegLength) 
	new.cutoffs = GetCutoffs(cutoffs, arm.segmean)
	gene.segmean = GetGeneSegmean(segs = split.segs, arm.segmean = arm.segmean, gene = gene, neighbor.distance = neighborDistance, min.seg.evidence = minSegEvidence)
	gene.score = AggregateSegmean(gene.segmean = gene.segmean, cutoffs = new.cutoffs)
	gene.score = gene.score[match(gene$Gene, gene.score$Gene)] 
	gene.score$one_score[which(gene.score$filter != "PASS")] = NA
	out[, i] = gene.score$one_score
}

saveRDS(out, file="6116.rds")




gistic = fread("~/SCRATCH/gistic/data/all_thresholded.by_genes/BRCA.all_thresholded.by_genes.txt")
m = match(rownames(d), gistic$"Gene Symbol")
gistic = gistic[m]
aliquot = intersect(colnames(gistic), colnames(new))
gistic = gistic[, aliquot, with=F]


d = readRDS("6116.rds")
d = d[, aliquot]
n = nrow(d) * ncol(d)
table(d, useNA="always") / n * 100
temp = table(as.matrix(d) - as.matrix(gistic), useNA="always") / n * 100
temp[names(temp) == "0"]
sum(temp[names(temp)  %in% c("0", "1", "-1")])
sum(temp[! names(temp)  %in% c("0", "1", "-1", NA)])



max ~ 1.85045
min ~ 1.30202

gistic


7227
        -2         -1          0          1          2       <NA>
 0.6461463 14.6607624 66.2230594 13.2782971  1.2424378  3.9492971
same 78.78448
similar 95.94168
diff 0.1090227


6226
        -2         -1          0          1          2       <NA>
0.01028773 0.14276693 0.66223059 0.12909240 0.01610227 0.03952007

          -4           -3           -2           -1            0            1
1.027248e-06 1.142058e-05 8.891734e-04 9.233743e-02 7.845477e-01 8.250499e-02
           2            3         <NA>
1.869742e-04 1.193420e-06 3.952007e-02
same 78.45477
similar 95.93901
diff 0.108983 


5225
       -2        -1         0         1         2      <NA>
1.649323 13.654042 66.223059 12.450781  2.066926  3.955868
same 77.80784
similar 95.93518
diff 0.1089502

4224
       -2        -1         0         1         2      <NA>
 2.335178 12.966083 66.223059 12.068763  2.447299  3.959616
same 77.06105
similar 95.93145
diff 0.1089351



out = d[0]
for (disease in diseases) {
	print(disease)
	f = paste("TCGA", disease, "nocnv_grch38.seg.v2.txt", sep=".")
	segs = ReadGenomicTSV(f, header = T, segment = T, keep.chry = F)
	split.segs = SplitSegByCentromere(segs, centromere)
	samples = unique(segs$Sample)
	d = data.table(project=disease, sample=samples, min=0, max=0)

	for(i in  1:length(samples)) {
		sample = samples[i]
		arm.segmean = GetArmSegmean(split.segs[Sample==sample], minArmEvidence, centromere, bed=NULL, minBedSegLength)
		d$min[i] = min(arm.segmean$segmean)
		d$max[i] = max(arm.segmean$segmean)	
	}
	out = rbind(out, d)
}


tcga.segs = fread("zcat All.TCGA.nocnv_grch38.seg.v2.txt.gz")
tcga.split.segs = SplitSegByCentromere(tcga.segs, centromere)
tcga.arm.segmean = out



# modify cutoffs by max and min of arm.segmean
new.cutoffs = GetCutoffs(cutoffs, arm.segmean)

# get all potential gene-level segmean values
arm.segmean = GetArmSegmean(split.segs, minArmEvidence, centromere, bed, minBedSegLength) 
gene.segmean = GetGeneSegmean(segs = split.segs, arm.segmean = arm.segmean, gene = gene, neighbor.distance = neighborDistance, min.seg.evidence = minSegEvidence)

spilt.segs = tcga.split.segs[Sample == sample]
run = function(split.segs, new.cutoffs = new.cutoffs, centromere, gene) {
	arm.segmean = GetArmSegmean(split.segs, min.arm.evidence = 4, centromere, bed = NULL, min.bed.segment.len = 50) 
	new.cutoffs = GetCutoffs(new.cutoffs, arm.segmean)
	gene.segmean = GetGeneSegmean(segs = split.segs, arm.segmean = arm.segmean, gene = gene, neighbor.distance = 1e6, min.seg.evidence = 50)	
	gene.score = AggregateSegmean(gene.segmean = gene.segmean, cutoffs = new.cutoffs)
	score = gene.score[match(gene$Gene, gene.score$Gene)]$one_score
	score[which(! gene.score$filter %in% c("impute", "PASS"))] = NA
	return(score)
}
run(split.segs, new.cutoffs, centromere, gene)


library(parallel)
num_parallel = 35


samples = colnames(gistic)

p1 = proc.time()
new.cutoffs = c(-0.6, -0.1, 0.1, 0.6)
cl = makeCluster(num_parallel, type="FORK")
out = parLapply(cl, samples, function(x) 
			run(tcga.split.segs[Sample == x], new.cutoffs, centromere, gene))	
stopCluster(cl)
p2 = proc.time()
p2 - p1
out = unlist(out)
dim(out) = c(59852, length(samples))
colnames(out) = samples
rownames(out) = gene$Gene
save(out, file="out.6116.rds")

temp = table(as.vector(gistic - out), useNA="always")




cl = makeCluster(num_parallel, type="FORK")
out = parLapply(cl, samples, function(x) 
			tcga.split.segs[Sample == x])	
stopCluster(cl)



fread(f[3])[1:3,1:3]

 run(tcga.split.segs[Sample == samples[2]], new.cutoffs, centromere, gene)[1:10]



# get aggregated gene-level CNV scores
gene.score = AggregateSegmean(gene.segmean = gene.segmean, cutoffs = new.cutoffs)

# add off_target filter
if(!is.null(bed.file)){
	gene.score = FilterByBed(gene.score = gene.score, bed = bed)
}

# format and output
WriteScores(gene.score, gene, filters, out.file)





# split segments by centromere 
split.segs = SplitSegByCentromere(segs, centromere)

# we assume bed file and gene file do not have regions that contains centromere 

# calcualte arm-level segmean value
arm.segmean = GetArmSegmean(split.segs, minArmEvidence, centromere, bed, minBedSegLength) 


> summary(lm(gistic_high ~ arm_max, data = a))
            Estimate Std. Error t value             Pr(>|t|)    
(Intercept) -0.03693    0.00753    -4.9           0.00000095 ***
arm_max      1.98210    0.01346   147.2 < 0.0000000000000002 ***

> summary(lm(gistic_high ~ arm_max + 0, data = a))
        Estimate Std. Error t value            Pr(>|t|)    
arm_max  1.92654    0.00728     265 <0.0000000000000002 ***


> summary(lm(gistic_low ~ arm_min, data = a))
(Intercept) -0.13593    0.00296   -45.9 <0.0000000000000002 ***
arm_min      1.06032    0.00579   183.3 <0.0000000000000002 ***

> summary(lm(gistic_low ~ arm_min + 0, data = a))
        Estimate Std. Error t value            Pr(>|t|)    
arm_min  1.28495    0.00336     383 <0.0000000000000002 ***

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