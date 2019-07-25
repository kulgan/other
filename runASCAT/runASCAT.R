# Usage Rscript runASCAT.R --vanilla tumor_uuid normal_uuid tumor_aliquot_id project_id gender

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
if(tolower(gender) == "male") {
  ascat.gender <- "XY"
} else {
  ascat.gender <- "XX"
}

# Input file
snp.file <- "zcat ~/SCRATCH/snp6.na35.remap.hg38.txt.gz"
ascat.gc.file <- "~/SCRATCH/GC_AffySNP6_GDC_GRCh38_20190128.txt"

####################################################
# Generate ASCAT Input LogR and BAF files
####################################################

# read snp file
snp <- fread(snp.file)
names(snp)[1:3] <- c("probeid", "chr", "pos")
snp$chr <- factor(snp$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))

####################################################
# Run ASCAT 
####################################################
# input (output from the previous WriteLogRBAF step)
tumor.logr.file <- paste0("../lrrbaf/", tumor.uuid, ".LogR.txt")
normal.logr.file <- paste0("../lrrbaf/", normal.uuid, ".LogR.txt")
tumor.baf.file <- paste0("../lrrbaf/", tumor.uuid, ".BAF.txt")
normal.baf.file <- paste0("../lrrbaf/", normal.uuid, ".BAF.txt")

# output 
prefix <- paste0(project.id, ".", aliquot.uuid)
segment.file <- paste0(prefix, ".ascat2.allelic_specific.seg.txt")
acf.file <- paste0(prefix, ".ascat2.acf.txt")
ploidy.file <- paste0(prefix, ".ascat2.ploidy.txt")

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

# clean other files
file.remove(paste0(tumor.uuid, ".LogR.PCFed.txt"))
file.remove(paste0(tumor.uuid, ".BAF.PCFed.txt"))
file.remove(paste0(tumor.uuid, ".sunrise.png"))
file.remove(paste0(tumor.uuid, ".rawprofile.png"))
file.remove(paste0(tumor.uuid, ".ASCATprofile.png"))




