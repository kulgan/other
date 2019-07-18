#!/usr/bin/env

# Seg2Gene
# author: Zhenyu Zhang
# date: 07/15/2019
# require: R > 3.40 
#	   	data.table
#			futile.logger
#			optparse
#			matrixStats



GetOptions = function() {
	# define option_list
	option.list = list(
		make_option(c("-f", "--file"), type = "character", default = NULL, 
	              	help = "integer value copy number segmentation file"),
	  make_option(c("-g", "--gene"), type = "character", default = NULL, 
	              	help = "gene location file name"), 
		make_option(c("-o", "--out"), type = "character", default = NULL, 
	              	help = "output gene level copy number file name"), 
	  make_option(c("-l", "--log"), type = "character", default = "stdout", 
	              	help = "log file [default = %default]")   
	)
	
	# parse 
	opt.parser = OptionParser(option_list = option.list);
	opt = parse_args(opt.parser);

	# if log output is not "stdout", write to log file
	if(opt$log != "stdout") {
		flog.appender(appender.file(opt$log))
	}

	# stop if input segmentation file is not provided
	if(is.null(opt$file) | file.exists(opt$seg)){
	  print_help(opt_parser)
	  flog.error("Input segmentation file was not provided")
	  stop("Input segmentation file is required\n", call.=FALSE)
	}

	# stop if gene file is not provided
	if(is.null(opt$gene) | file.exists(opt$gene)){
	  print_help(opt_parser)
	  flog.error("Gene location file was not provided")
	  stop("Gene location is required\n", call.=FALSE)
	}

	# stop if exon output files is not provided & default file name already exists
	if(is.null(opt$out)) {
		opt$out = paste0("./", basename(as.character(opt$file)), ".gene_level_copy_number.tsv")
		if(file.exists(opt$out)) {
			print_help(opt_parser)
	  	flog.error("Output file name was not provided and default output %s already exist", opt$out)
	  	stop("Output file name is required\n", call.=FALSE)
		}
	}

	# display all input options
	flog.info("input segmentation file: \t%s", opt$file)
	flog.info("intput gene location file: \t%s", opt$gene)
	flog.info("output gene level cnv file: \t%s", opt$out)

	# return
	return(opt)
}

# ReadSegment reads and reformat input segmentation file
ReadSegment <- function(file, log) {

	flog.info("Reading copy number segmentation file")

	# define dfault column names nad chromosome names
	col.names <- c( "chromosome", "start", "end", "copy_number")

	# if log output is not "stdout", write to log file
	if(opt$log != "stdout") {
		flog.appender(appender.file(opt$log))
	}

# read data and change column names into lower case
	data <- fread(file)
	setnames(data, tolower(names(data)))
	flog.info("- %s segments read", nrow(data))

# check if required columns exist:
# if yes, only keep those required columns;
# if no, keep column 2 ~ 5, and set them as required columns
	if(sum(! col.names %in% names(data)) == 0) {
		data <- data[, col.names, with=F]
	} else {
		data <- data[, 2:5]
		setnames(data, col.names)
	}

	# format column names into chr1-22/X/Y, and drop others
	data$chromosome <- paste0("chr", gsub("^chr", "", data$chromosome))
	data$chromosome[which(data$chromosome == "chr23")] <- "chrX"
	data$chromosome[which(data$chromosome == "chr24")] <- "chrY"

	# return
	return(data)
}


# ReadSegment reads and reformat input gene location file
ReadGene <- function(file, log) {
	flog.info("Reading gene location file")

	# define expected column names nad chromosome names
	col.names <- c( "chromosome", "start", "end", "gene_id", "gene_name")

	# if log output is not "stdout", write to log file
	if(opt$log != "stdout") {
		flog.appender(appender.file(opt$log))
	}

# read data and change column names into lower case
	data <- fread(file)
	setnames(data, tolower(names(data)))
	flog.info("- %s genes read", nrow(data))

# check if required columns exist
# if yes, only keep those required columns
# if not, output error message
	if(sum(! col.names %in% names(data)) == 0) {
		data <- data[, col.names, with=F]
	} else {
		flog.error("The expected gene location file should have at least the following columns: %s", paste0(col.names, collapse=","))
		stop("Program aborted. Please provide a gene location file with expected format\n", call.=FALSE)
	}

	# format column names into chr1-22/X/Y, and drop others
	data$chromosome <- paste0("chr", gsub("^chr", "", data$chromosome))
	data$chromosome[which(data$chromosome == "chr23")] <- "chrX"
	data$chromosome[which(data$chromosome == "chr24")] <- "chrY"

	# return
	return(data)
}


GetGeneLevel <- function(seg, gene, contigs, ties, log) {
	flog.info("Generating gene-level copy number calls")

	# only keep segments and genes in defined "contigs"
	seg <- seg[chromosome %in% contigs]
	gene <- gene[chromosome %in% contigs]

	# find segmentation overlaps and calcualte overlapped length
	setkey(seg, chromosome, start, end)
	ovlp <- foverlaps(gene, seg, type = "any", nomatch=NA, mult = "all")
	ovlp$size <- with(ovlp, pmin(end, i.end) - pmax(start, i.start) + 1)

	# calcualte gene-level copy number based on weighted median of the overlapped regions
	cnv <- ovlp %>% 
		group_by(gene_id) %>% 
		summarise(gene_name = gene_name[1], 
						 	chromosome = chromosome[1], 
						 	start = i.start[1],
						 	end = i.end[1],  
						 	cnv = weightedMedian(x = copy_number, w = size, ties = ties), 
						 	min_copy_number = min(copy_number), 
							max_copy_number = max(copy_number)) %>%
		rename(copy_number = cnv) %>% 
		mutate(chromosome = factor(chromosome, levels = contigs)) %>%
		arrange(chromosome, start, end) %>%
		setDT()

	# information output
	flog.info("- %s genes are kept", nrow(cnv))
	flog.info("- %s genes have assgined copy number values", sum(!is.na(cnv$copy_number)))
	flog.info("  - %s genes have one uniform copy number acrossing the gene region", nrow(cnv[min_copy_number == max_copy_number]))
	flog.info("  - %s genes overlap with segments with different copy numbers, and a unique value is used by weighted median method", nrow(cnv[min_copy_number != max_copy_number]))
	flog.info("- %s genes do not overlap with any input copy number segments, thus do not have a copy number value", sum(is.na(cnv$copy_number)))

	# return
	return(cnv)
}


##############################
# Main Program
##############################

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(futile.logger)))
suppressWarnings(suppressMessages(library(matrixStats)))

options(scipen = 999, digits = 0)


# pre-defined variables
# contigs: which contigs will be kept from final output
# ties: ties method from matrixStats::weightedMedian to break ties when two possible median values exist. Acceptable values are "min", "max", "mean", "weighted"
contigs <- paste0("chr", c(1:22, "X", "Y"))
ties = "min"

# get options
opt <- GetOptions()

# if log output is not "stdout", write to log file
if(opt$log != "stdout") {
	flog.appender(appender.file(opt$log))
}

# read segmentation file
seg <- ReadSegment(file = opt$file, log = opt$log) 

# read gene file
gene <- ReadGene(file = opt$gene, log = opt$log)

# calculate gene-level copy number
cnv <- GetGeneLevel(seg = seg, gene = gene, contigs = contigs, ties = ties, log = opt$log) 

# write to disk
fwrite(cnv, opt$out, col.names=T, row.names=F, sep="\t", quote=F)


