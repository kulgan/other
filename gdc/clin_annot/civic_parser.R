
# example command to run this script RScript aggregator.R --gtf gencode.v22.annotation.gtf.gz
# setwd("/Users/zhenyu/github/other/gdc/clin_annot")

# GetOptions read parameters from input
GetOptions = function() {
  # define option_list
  option_list = list(
    make_option(c("-g", "--gtf"), type = "character", default = NULL, 
                  help = "path to the GTF file"),
    make_option(c("-c", "--checkURL"), type = "logical", default = F, 
                  help = "whether to check the existance of output URLs [default = %default]"), 
    make_option(c("-d", "--wkdir"), type = "character", default = "auto", 
                  help = "path to the working directory where outputs are saved (Attention: file overwritten risk) [default = %default]"), 
    make_option(c("-l", "--log"), type = "character", default = "log.txt", 
                  help = "log file (file name only, no path) [default = %default]")   
  )
  
  # parse 
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)

  # normalize path
  opt$gtf = normalizePath(opt$gtf)
  if(opt$wkdir != "auto") {
    opt$wkdir = normalizePath(opt$wkdir)
  }

  # stop if input gtf file is not provided
  if (is.null(opt$gtf)){
    print_help(opt_parser)
    flog.error("no input GTF file provided")
    stop("At least one argument must be supplied (input file)\n", call.=FALSE)
  }

  # display all input options
  flog.info("GTF File: \t\t\t%s", opt$gtf)
  flog.info("CheckURL: \t\t\t%s", opt$checkURL)
  flog.info("Working Directory: \t%s", opt$wkdir)
  flog.info("Log File: \t\t\t%s", opt$log)

  # return
  return(opt)
}

# CreateWkDir creates and sets the current working diretory to be ./run_YYYYMMDD/
# If it already exists, a underscore and a number is added as postfix of the directory name, starting from "_1"
CreateWkDir <- function(working.dir) {
  flog.info("Initialization ...")

  if (working.dir == "auto") {
    current.date <- format(Sys.time(), "%Y%m%d")
    working.dir <- paste0(getwd(), "/run_", current.date)

    # loop until the new directory name does not exist
    i <- 1
    new.dir <- working.dir
    while(dir.exists(new.dir)) {
      new.dir <- paste0(working.dir, "_", i) 
      i <- i + 1
    }
    working.dir <- new.dir
  }

  # create directory and setwd
  dir.create(working.dir)
  setwd(working.dir)
  flog.info("%s is created and set as output directory", working.dir)
}

# ExtractGtfAttr accepts an array of GTF attribute, and attribute feature name, and return attribute values
ExtractGtfAttr <- function(gtf.attr, feature) {
  pattern <- paste0(".*", feature, " \"((\\w|\\.|-)+)\".*")
  value <- gsub(pattern, "\\1", gtf.attr)

  # when no pattern exist, set value to NA
  value[which(value == gtf.attr)] = NA
  return(value)
}

# GetGeneModel accept input of a GTF file, and return a list of gene data.table and transcript data.table
GetGeneModel <- function(gtf.file) {
  # parse gencode gene and transcript
  flog.info("  - Reading GTF file")

  # special treatment to gz file to allow reading
  if(file_ext(gtf.file) == "gz") {
    if(Sys.info()[['sysname']] == "Darwin") {
      # if we have MAC, use gunzip 
      gtf.file = paste("gunzip -c", gtf.file)
    } else if (Sys.info()[['sysname']] == "Linux") {
      # if we have Linux, use zcat
      gtf.file = paste("zcat", gtf.file)
    } else {
      flog.error("Unsupported OS")
    }
  }

  gtf <- fread(gtf.file, h=F)
  colnames(gtf) <- c("seqname", "source", "feature", "start", "end",  "score", "strand", "frame", "attribute")

  # select genes and extract gene_name/ gene_id 
  flog.info("  - Extracting gene features")
  gene <- gtf %>% 
    filter(feature == "gene") %>%
    mutate(gdc.gene.id = ExtractGtfAttr(attribute, "gene_id"), 
           gdc.gene.name = ExtractGtfAttr(attribute, "gene_name")) %>% 
    setDT()
  flog.info("    %s genes detected", nrow(gene))       

  flog.info("  - Extracting transcript features")
  transcript <- gtf %>%
    filter(feature == "transcript") %>%
    mutate(gdc.gene.id = ExtractGtfAttr(attribute, "gene_id"), 
           gdc.gene.name = ExtractGtfAttr(attribute, "gene_name"),
           gdc.transcript.id = ExtractGtfAttr(attribute, "transcript_id"), 
           gdc.isoform = gsub("^(\\w+).\\d+$", "\\1", gdc.transcript.id)) %>% 
    select(gdc.gene.id, gdc.gene.name, gdc.transcript.id, gdc.isoform) %>%
    setDT()                 
  flog.info("    %s transcripts detected", nrow(transcript))

  # return list of gene and transcript data tables
  out <- list(gene, transcript)
  names(out) <- c("gene", "transcript")
  return(out)
}




# ProcessCivic downloads and reformats CIVIC JSON raw data from https://civicdb.org/api/variants?count=
# URLs are then created as https://civicdb.org/events/genes/{gene_id}/summary/variants/{id}/summary
ProcessCivic <- function(allowedTypes, gene.model, max.record = 99999) {  

  # Download raw JSON
  civic.variant.url <- paste0("https://civicdb.org/api/variants?count=", max.record)
  flog.info("  - Downloading CIVIC raw variant JSON file")    
  civic.raw.file <- "civic.raw.json"
  download.file(civic.variant.url, civic.raw.file)
  flog.info("    Raw variant file saved as %s", civic.raw.file)  

  # Extract variant records
  flog.info("  - Reading CIVIC variant sets") 
  civic.raw <- jsonlite::fromJSON(civic.raw.file)
  expected.count = civic.raw$`_meta`$total_count
  records = civic$records
  count = nrow(records)

  # Validation of counts
  if (expected.count >  max.record) {
    flog.warn("Number of variants in CIVIC is more than %s. Please change the default parameter max.record to be larger", max.record)
  }
  if (expected.count != count) {
    flog.warn("Number of variants downloaded (%s) does not match the expected count (%s) from field _meta`$total_count", count, expected.count)
  }

  # In the records, one variant could have more than one variant_types
  # the following to flatten variant_types and extract "name" and "display_name"
  flog.info("  - Flattening JSON variant_types fields")    
  max.variant.type = max(sapply(records$variant_types, function(x) nrow(x)))
  if(max.variant.type > 2) {
    flog.error("Detected at least one variant with more than 2 variant types. Script needs to be updated to parse more than two variant_types")
  }
  # Extract variant_types$name
  variant.type.names <- sapply(records$variant_types, function(x) x$name)
  records$variant.type.name <- unlist(lapply(variant.type.names, function(x) ifelse(is.null(x), NA, x[1])))
  records$variant.type.name2 <- unlist(lapply(variant.type.names, function(x) ifelse(is.null(x), NA, x[2])))
  # Extract variant_types$displayname
  variant.type.display.names <- sapply(records$variant_types, function(x) x$display_name)
  records$variant.type.display.name <- unlist(lapply(variant.type.display.names, function(x) ifelse(is.null(x), NA, x[1])))
  records$variant.type.display.name2 <- unlist(lapply(variant.type.display.names, function(x) ifelse(is.null(x), NA, x[2])))
  # Additional information
  t = table(sapply(variant.type.names, length))
  flog.info("    %s variants have no variant_type names", ifelse(length(t[names(t) == 0]) == 0, 0, t[names(t) == 0]))
  flog.info("    %s variants have one variant_type names", ifelse(length(t[names(t) == 1]) == 0, 0, t[names(t) == 1]))
  flog.info("    %s variants have two variant_type names", ifelse(length(t[names(t) == 2]) == 0, 0, t[names(t) == 2]))

  # Flatten coordinates and extract "representative_transcript"
  records = cbind(records, records$coordinates) %>% 
    select(-c(variant_types, coordinates)) %>% 
    setDT()

write.table(records, "all.civic.variant.records.txt", col.names=T, row.names=F, sep="\t", quote=F)



  civic$transcript <-civic$coordinates$representative_transcript
  civic$transcript2 <- civic$coordinates$representative_transcript2

  # Drop columns of variant_types, coordinates and description
  civic.complete <- subset(civic, select=-c(variant_types, coordinates)) %>%
    mutate(variant.type.name = ifelse(variant.type.name == "N/A", NA, variant.type.name)) 

  # add index for easy debug and save 
  civic.complete$index <- 1:nrow(civic.complete)
  civic.complete.file <- "civic.complete.tsv"
  fwrite(civic.complete, civic.complete.file, sep="\t", quote=F)
  flog.info("    %s variants read", nrow(civic.complete))
  flog.info("    Complete variant file saved as %s", civic.complete.file)

  # Fix errors in CIVIC raw data 
  flog.info("  - Transforming CIVIC variant sets")    
  civic <- civic.complete %>% 
    select(-description) %>%
    mutate(alteration = gsub(" *\\((\\w|\\.| |_|-|>|\\+|\\?)+\\)", "", name), # removing tailing bracket comments
           alteration = gsub("^p\\.", "", alteration), # remove p. from begining             
           alteration = gsub("^([ACDEFGHIKLMNPQRSTVWY]\\d+([ACDEFGHIKLMNPQRSTVWY])?fs)Ter(\\d+)$", "\\1*\\3", alteration), # change Ter into *
           alteration = gsub("^([ACDEFGHIKLMNPQRSTVWY]\\d+([ACDEFGHIKLMNPQRSTVWY])?)FS((\\*\\d+)?)$", "\\1fs\\3", alteration), # change FS into fs
           alteration = gsub("^(([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+)DUP$", "\\1dup", alteration), # change DUP into dup
           alteration = gsub("^(([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+)DEL$", "\\1del", alteration), # change DEL into del
           alteration = gsub("^(([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+)INS([ACDEFGHIKLMNPQRSTVWY*]+)$", "\\1ins\\3", alteration), # change INS into ins
           alteration = gsub("^(([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+)DELins([ACDEFGHIKLMNPQRSTVWY*]+)$", "\\1delins\\3", alteration), # change DELins into delins
           alteration = gsub("^([ACDEFGHIKLMNPQRSTVWY\\*]\\d+)X$", "\\1*", alteration), # change terminal X into *
           alteration = case_when(
             alteration == "Asn67fs"    ~ "N67fs", 
             alteration == "Glu34Lys"   ~ "E34K", 
             alteration == "EZH2 Y641F" ~ "Y641F", 
             alteration == "MLL-MLLT3"  ~ "KMT2A-MLLT3", 
             alteration == "BCR-ABL"    ~ "BCR-ABL1", 
             TRUE                       ~ alteration
            )
          ) 

  # Categorize variants by variant types
  civic <- civic %>%
    mutate(type = case_when(
      grepl("^[ACDEFGHIKLMNPQRSTVWY\\*]\\d+[ACDEFGHIKLMNPQRSTVWY\\*]$", alteration)                                  ~ "point_mutation", 
      grepl("^\\w+ (FUSIONS|Fusions|fusions)$", alteration)                                                          ~ "fusion_ns", 
      grepl("^(\\w|-)+ (FUSION|Fusion|fusion)$", alteration)                                                         ~ "fusion",
      grepl("^(\\w|-)+$", alteration) & ( grepl("fusion", variant.type.name) | grepl("fusion", variant.type.name2))  ~ "fusion", 
      grepl("^[ACDEFGHIKLMNPQRSTVWY]\\d+([ACDEFGHIKLMNPQRSTVWY])?fs(\\*\\d+)?$", alteration)                         ~ "frameshift",
      grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+dup$", alteration)                             ~ "duplication",
      grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+del$", alteration)                             ~ "small_deletion",
      grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+ins[ACDEFGHIKLMNPQRSTVWY*]+$", alteration)     ~ "insertion",
      grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+delins[ACDEFGHIKLMNPQRSTVWY*]+$", alteration)  ~ "delins",
      grepl("^[A-Z0-9]+-[A-Z0-9]+$", alteration)                                                                     ~ "fusion",                              
      alteration == "OVEREXPRESSION"                                                                                 ~ "overexpression",
      alteration == "UNDEREXPRESSION"                                                                                ~ "underexpression",
      alteration == "DELETION"                                                                                       ~ "large_deletion",
      alteration == "EXPRESSION"                                                                                     ~ "expression",
      alteration == "AMPLIFICATION"                                                                                  ~ "amplification",
      alteration == "MUTATION"                                                                                       ~ "mutation_ns",
      alteration == "LOSS-OF-FUNCTION"                                                                               ~ "loss-of-function",
      alteration == "TRUNCATING MUTATION"                                                                            ~ "truncation",
      alteration == "PHOSPHORYLATION"                                                                                ~ "phosphorylation",
      alteration == "LOSS"                                                                                           ~ "loss",
      alteration == "METHYLATION"                                                                                    ~ "methylation",
      alteration == "NUCLEAR EXPRESSION"                                                                             ~ "nuclear-expresssion",
      alteration == "CYTOPLASMIC EXPRESSION"                                                                         ~ "cytoplasmic-expression",
      alteration == "FRAMESHIFT MUTATION"                                                                            ~ "frameshift_ns",
      alteration == "Splicing alteration"                                                                            ~ "splicing_ns",
      TRUE                                                                                                           ~ "other"
      )
    ) %>%                              
    mutate(alteration = gsub("^((\\w|-)+) (FUSION|Fusion|fusion)(S|s)?$", "\\1", alteration)) # remove tailing "fusions"


  # Split gene fusions into two gene columns: gene and fusion.gene
  # "*" in fusion.gene represents any partner gene 
  civic <- civic %>% 
    rowwise() %>%
    mutate(fusion.gene = ifelse(type == "fusion", gsub("(-$|^-)", "", gsub(entrez_name, "", alteration)), NA),
           fusion.gene = ifelse(type == "fusion_ns", "*", fusion.gene)) 

  # Save parsed data to disk
  civic.transformed.file <- "civic.transformed.tsv"
  fwrite(civic, civic.transformed.file, sep="\t", quote=F)
  flog.info("    Transformed variant file saved as %s", civic.transformed.file)

  # Subset to only keep those GDC will use for annotations
  flog.info("  - Subsetting CIVIC variant sets")    
  civic.subset <- civic %>% 
    filter(type %in% allowedTypes) %>% 
    select(id, name, entrez_name, gene_id, transcript, transcript2, alteration, type) %>%
    setDT()

  # Save subset data to disk
  civic.subset.file <- "civic.subset.tsv"
  fwrite(civic.subset, civic.subset.file, sep="\t", quote=F)
  flog.info("    %s variants selected", nrow(civic.subset))
  flog.info("    Selected variant file saved as %s", civic.subset.file)


  # Add URLs
  flog.info("  - Constructing CIVIC urls")  
  civic.url <- civic.subset %>% 
    mutate(isoform = gsub("^(\\w+).\\d+$", "\\1", transcript)) %>%
    left_join(gene.model[["transcript"]], by = c("isoform" = "gdc.isoform")) %>%
    rename(gdc.gene.name.from.isoform = gdc.gene.name, 
           gdc.gene.id.from.isoform = gdc.gene.id) %>%
    left_join(gene.model[["gene"]], by = c("entrez_name" = "gdc.gene.name")) %>%
    rename(gdc.gene.id.from.gene = gdc.gene.id) %>% 
    rowwise() %>% 
    mutate(gdc.gene.name = ifelse(is.na(gdc.gene.id.from.gene), gdc.gene.name.from.isoform, entrez_name), 
           gdc.gene.id = ifelse(is.na(gdc.gene.id.from.gene),  gdc.gene.id.from.isoform, gdc.gene.id.from.gene),
           hgvsp = paste0("p.", alteration),
           civic.url = URLencode(paste0("https://civicdb.org/events/genes/", gene_id, "/summary/variants/", id, "/summary"))) %>% 
    setDT()

  # Check genes that do not have a match 
  left <- civic.url[is.na(gdc.gene.id)]
  if (nrow(left) > 0) {
    flog.warn("Can not find gene name match for %s", left$entrez_name)
  }

  # select columns
  civic.url <- civic.url %>%    
    select(gdc.gene.name, gdc.gene.id, gdc.transcript.id, hgvsp, type, civic.url, name) %>% 
    setDT()
 
  # Save final data to disk
  civic.final.file <- "civic.final.tsv"
  fwrite(civic.url, civic.final.file, sep="\t", quote=F)
  flog.info("    Final variant file saved as %s", civic.final.file)

  # return
  return(civic.url)
}


# CheckURL is a wrapper of url.exists 
# Check URL checks if an array of URLs exist, and return an array of booleans
CheckURL <- function(urls,  batch = 100) {
  n = length(urls)
  out = array(NA, n)
  flog.info("    Total %s URLs to check", n)

  # Split by batches to avoid "Error: Unable to establish connection with R session"
  for(i in 1: ceiling(n / batch)) {
    start = (i - 1) * batch + 1
    end = min(i * batch, n) 
    flog.info("    Checking URLs from %s to %s", start, end)
    out[start:end] = sapply(urls[start:end], url.exists)
  }

  # return results
  return(out)
}


suppressMessages(library(futile.logger))
suppressMessages(library(dplyr))
suppressMessages(library(RCurl))
suppressMessages(library(data.table))
suppressMessages(library(jsonlite))
suppressMessages(library(XML))
suppressMessages(library(tools))
suppressMessages(library(optparse))

# library(devtools)
# library(hgvsParseR)

# if gtf file is gz compressed, need zcat install in Linux or gunzip in OSX

# Init
options(scipen=999)
allowedTypes <- c("point_mutation", "duplication", "small_deletion", "insertion", "delins", "frameshift", "splicing")

# Read Options
opt = GetOptions()

# Create and set output dir
CreateWkDir(opt$wkdir)
flog.appender(appender.file(opt$log))

# Set log file
log.file <- opt$log
flog.appender(appender.file(log.file)) 
flog.info("Writing logs in %s in the working directory", log.file)

# extract gene model
flog.info("Preparing Gene Model ...")
gene.model <- GetGeneModel(opt$gtf)

# Process OncoKB
flog.info("Processing OncoKB variants ...")
oncokb <- ProcessOncokb(allowedTypes, gene.model)

# Process CIVIC
flog.info("Processing CIVIC variants ...")
civic <- ProcessCivic(allowedTypes, gene.model)

# Process MyCancerGenome
flog.info("Processing MyCancerGenome variants ...")
mcg <- ProcessMcg(allowedTypes, gene.model)

# Merge
# Please Note gene name are not unique, so ensembl id should be used
# For example, CRLF2 matches both ENSG00000205755.9 and ENSGR0000205755.9 on chrX and chrY
flog.info("Merging databases ...")
variants <- rbind(
  oncokb %>% 
    select(gdc.gene.name, gdc.gene.id, hgvsp, oncokb.url) %>%
    rename(url = oncokb.url) %>% 
    mutate(name = paste(gdc.gene.name, gsub("^p.", "", hgvsp)), 
           source = "OncoKB"), 
  civic %>% 
    select(gdc.gene.name, gdc.gene.id, hgvsp, civic.url, name) %>%
    rename(url = civic.url) %>% 
    mutate(source = "CIVIC"), 
  mcg %>% 
    select(gdc.gene.name, gdc.gene.id, hgvsp, mcg.url, name) %>%
    rename(url = mcg.url) %>% 
    mutate(source = "MyCancerGenome")) %>% 
  rename(gene.name = gdc.gene.name, 
         gene.id = gdc.gene.id) %>%
  setDT()

# Validate URLs
# This is extremely slow for CIVIC and MyCancerGenome (~ hour)
# Suggest to ignore during testing
if(opt$checkURL) {
  flog.info("Validating all URLs ...")
  validation <- CheckURL(variants$url)
  flog.info("    %s out of %s URLs are valid", sum(validation), nrow(variants))
  if(sum(!validation) > 0) {
    flog.warn("    The following URLs are invalid %s", variants[which(!validation)]$url)
  }
} else {
  flog.info("URL validation is skipped")
}

# write final data to disk
all.variants.file <- "all.variants.txt"
fwrite(variants, all.variants.file, sep="\t", quote=F)

# final summary output to the log file
flog.info("Summary")
flog.info("  - %s variants in total", nrow(variants))
flog.info("  - %s unique variants", length(unique(paste(variants$gene.id, variants$hgvsp))))
flog.info("  - %s variants from OncoKB", sum(variants$source == "OncoKB"))
flog.info("  - %s variants from CIVIC", sum(variants$source == "CIVIC"))
flog.info("  - %s variants from MyCancerGenome", sum(variants$source == "MyCancerGenome"))




