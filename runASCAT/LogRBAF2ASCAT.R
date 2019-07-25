# Usage Rscript LogRBAF2ASCAT.R tumor_uuid normal_uuid tumor_aliquot_id project_id gender

# read {uuid}.lrr_baf.txt.gz file, reformat and write LogR and BAF file to disk
WriteLogRBAF <- function(uuid, snp, aliquot.id){
  # input
  lrrbaf.file <- paste0("zcat /home/ubuntu/SCRATCH/snp6cnv/output/", uuid, ".lrr_baf.txt.gz")

  # output  
  logr.file <- paste0(aliquot.id, ".", uuid, ".LogR.txt")
  baf.file <- paste0(aliquot.id, ".", uuid, ".BAF.txt")

  # read lrrbaf file and merge with snp file
  lrrbaf <- fread(lrrbaf.file)
  setnames(lrrbaf, c("probeid", "chr", "pos", "lrr", "baf"))
  lrrbaf <- merge(lrrbaf[, -c("chr", "pos"), with=F], snp, all=F)

  # Set CN probe BAF to be NA
  cn.index <- which(lrrbaf$type == "CN")
  lrrbaf$baf[cn.index] <- NA

  # Normalize LogR of CN and SNP probes separately
  cn.mean <- mean(lrrbaf$lrr[cn.index], na.rm = T)
  lrrbaf$lrr[cn.index] <- lrrbaf$lrr[cn.index] - cn.mean

  snp.mean <- mean(lrrbaf$lrr[-cn.index], na.rm = T)
  lrrbaf$lrr[-cn.index] <- lrrbaf$lrr[-cn.index] - snp.mean

  # order lrrbaf file by chr and position
  lrrbaf <- lrrbaf[order(chr, pos)]

  # reformat and write logR file to disk
  lrr <- lrrbaf[, c("probeid", "chr", "pos", "lrr")]
  setnames(lrr, c("", "chrs", "pos", aliquot.id))
  fwrite(lrr, logr.file, col.names = T, row.names = F, sep="\t", quote=F)

  # reformat and write BAF file to disk
  baf <- lrrbaf[, c("probeid", "chr", "pos", "baf")]
  setnames(baf, c("", "chrs", "pos", aliquot.id))
  fwrite(baf, baf.file, col.names = T, row.names = F, sep="\t", quote=F)
}


library(data.table)
library(ASCAT)

# get options
args <- commandArgs(trailingOnly=TRUE)
tumor.uuid <- args[1]
normal.uuid <- args[2]
aliquot.uuid <- args[3]
project.id <- args[4]
gender <- args[5]

# transform gender to ASCAT format
if(is.null(gender)){
  ascat.gender <- "XX"
} else if (tolower(gender) != "male") {
  ascat.gender <- "XX"
} else {
  ascat.gender <- "XY"
}

# Input file
snp.file <- "zcat /home/ubuntu/SCRATCH/snp6cnv/snp6.na35.remap.hg38.txt.gz"
ascat.gc.file <- "/home/ubuntu/SCRATCH/snp6cnv/GC_AffySNP6_GDC_GRCh38_20190128.txt"

####################################################
# Generate ASCAT Input LogR and BAF files
####################################################

# read snp file
snp <- fread(snp.file)
names(snp)[1:3] <- c("probeid", "chr", "pos")
snp$chr <- factor(snp$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))

# generate LogR and BAF file from PennCNV output
WriteLogRBAF(tumor.uuid, snp, aliquot.uuid)
WriteLogRBAF(normal.uuid, snp, aliquot.uuid)

####################################################
# Run ASCAT 
####################################################
# input (output from the previous WriteLogRBAF step)
tumor.logr.file <- paste0(aliquot.uuid, ".", tumor.uuid, ".LogR.txt")
normal.logr.file <- paste0(aliquot.uuid, ".", normal.uuid, ".LogR.txt")
tumor.baf.file <- paste0(aliquot.uuid, ".", tumor.uuid, ".BAF.txt")
normal.baf.file <- paste0(aliquot.uuid, ".", normal.uuid, ".BAF.txt")

# output 
prefix <- paste0(project.id, ".", aliquot.uuid)
segment.file <- paste0(prefix, ".ascat2.allelic_specific.seg.txt")
acf.file <- paste0(prefix, ".ascat2.acf.txt")
ploidy.file <- paste0(prefix, ".ascat2.ploidy.txt")
error.log.file <- paste0(prefix, ".log.txt")

chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
# load input
ascat.bc <- ascat.loadData(tumor.logr.file, tumor.baf.file, normal.logr.file, normal.baf.file, chrs = chrs, gender = ascat.gender, sexchromosomes = c("chrX", "chrY"))
# GC correction
ascat.bc <- ascat.GCcorrect(ascat.bc, ascat.gc.file)
# ascat.plotRawData(ascat.bc, img.prefix = prefix)
# ASPCF segmentation
ascat.bc <- ascat.aspcf(ascat.bc)
# ascat.plotSegmentedData(ascat.bc, img.prefix = prefix)
# run ASCAT
tryCatch({
  ascat.output <- ascat.runAscat(ascat.bc)
  # output segment 
  segment <- ascat.output$segments
  setnames(segment, c("GDC_Aliquot", "Chromosome", "Start", "End", "Major_Copy_Number", "Minor_Copy_Number"))
  segment$Total_Copy_Number <- segment$Major_Copy_Number + segment$Minor_Copy_Number
  write.table(segment, segment.file, sep = "\t", quote = F, row.names = F, col.names = T)
  # output aberrantcellfraction
  acf <- data.frame(ascat.output$aberrantcellfraction)
  write.table(acf, acf.file, sep = "\t", quote = F, row.names = T, col.names = F)
  # output ploidy
  ploidy <- data.frame(ascat.output$ploidy)
  write.table(ploidy, ploidy.file, sep = "\t", quote = F, row.names = T, col.names = F)
}, warning = function(war) {
  msg = "ASCAT could not find an optimal ploidy and cellularity value"
  if(grepl(msg, war$message)){
    cat(msg, file=error.log.file, sep="\n")
  } else {
    cat("other warnings", file=error.log.file, sep="\n")
  }
}, error = function(err) {
  cat("other errors", file=error.log.file, sep="\n")
}, finally = {
  # cleanup
  system(paste0("rm ", aliquot.uuid, ".*"))
})




