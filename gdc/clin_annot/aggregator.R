
# example command to run this script RScript aggregator.R --gtf gencode.v22.annotation.gtf.gz
# setwd("/Users/zhenyu/github/other/gdc/clin_annot")

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

  # normlize path
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
  flog.info("Working Directory: \t\t%s", opt$wkdir)
  flog.info("Log File: \t%s", opt$log)

  # return
  return(opt)
}

# CreateWkDir creates and sets the current working diretory to be ./runYYYYMMDD/
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


# OncoKB 
# URL http://oncokb.org/#/gene/{Gene}/alteration/{Alteration}  
# http://oncokb.org/api/v1/info
ProcessOncokb <- function(allowedTypes, gene.model) {
  oncokb.variant.url = "http://oncokb.org/api/v1/utils/allAnnotatedVariants.txt"

  # download oncokb raw variant txt file
  flog.info("  - Downloading OncoKB raw variant TXT file")  
  oncokb.raw.file <- "oncokb.raw.txt"
  download.file(oncokb.variant.url, oncokb.raw.file)
  flog.info("    Raw variant file saved as %s", oncokb.raw.file)

  # read into data.table, and save to disk
  flog.info("  - Reading OncoKB variant sets") 
  oncokb.complete <- fread(oncokb.raw.file)
  flog.info("    %s variants read", nrow(oncokb.complete))

  # make column name easy to use
  setnames(oncokb.complete, gsub(" ", "_", tolower(names(oncokb.complete))))
  # add index column for debug purpose
  oncokb.complete$index = 1:nrow(oncokb.complete)
  oncokb.complete.file <- "oncokb.complete.tsv"
  fwrite(oncokb.complete, oncokb.complete.file, sep="\t", quote=F)
  flog.info("    Complete variant file saved as %s", oncokb.complete.file)

  # Fix errors in OncoKB raw data 
  flog.info("  - Transforming OncoKB variant sets")    
  oncokb <- oncokb.complete %>% 
    mutate(protein_change = case_when(
      protein_change == "AGAP3?BRAF Fusion"      ~ "AGAP3-BRAF Fusion",
      protein_change == "CIITA?BX648577 Fusion"  ~ "CIITA-BX648577 Fusion",
      protein_change == "FGFR2?PPHLN1 Fusion"    ~ "FGFR2-PPHLN1 Fusion",
      protein_change == "PAX8-PPAR? Fusion"      ~ "PAX8-PPARG Fusion",
      protein_change == "ERLIN2?FGFR1 Fusion"    ~ "ERLIN2-FGFR1 Fusion", 
      TRUE                                       ~ protein_change
      )
    ) 

  # Categorize variants by variant types
  oncokb <- oncokb %>% 
    mutate(type = case_when(
      grepl("^[ACDEFGHIKLMNPQRSTVWY\\*]\\d+[ACDEFGHIKLMNPQRSTVWY\\*]$", protein_change)                                  ~ "point_mutation",
      grepl("^(\\w|-)+ Fusion$", protein_change)                                                                         ~ "fusion",
      grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+dup$", protein_change)                             ~ "duplication",
      grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+del$", protein_change)                             ~ "small_deletion",
      grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+ins[ACDEFGHIKLMNPQRSTVWY*]+$", protein_change)     ~ "insertion",
      grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+delins[ACDEFGHIKLMNPQRSTVWY*]+$", protein_change)  ~ "delins",
      grepl("^[ACDEFGHIKLMNPQRSTVWY]\\d+([ACDEFGHIKLMNPQRSTVWY])?fs(\\*\\d+)?$", protein_change)                         ~ "frameshift",
      grepl("^[ACDEFGHIKLMNPQRSTVWYX]\\d+_splice$", protein_change)                                                      ~ "splicing",
      protein_change == "Amplification"                                                                                  ~ "amplification",
      protein_change == "Deletion"                                                                                       ~ "large_deletion",
      protein_change == "Truncating Mutations"                                                                           ~ "truncation",
      protein_change == "Fusions"                                                                                        ~ "fusion_ns",
      protein_change == "Overexpression"                                                                                 ~ "overexpression", 
      TRUE                                                                                                               ~ "other"
      )
    )

  # Split gene fusions into two gene columns: gene and fusion.gene
  # "*" in fusion.gene represents any partner gene 
  oncokb <- oncokb %>% 
      rowwise() %>%
      mutate(fusion.gene = ifelse(type == "fusion", gsub("(-$|^-)", "", gsub(gene, "", gsub(" Fusion", "", protein_change))), NA), 
             fusion.gene = ifelse(type == "fusion_ns", "*", fusion.gene)) %>% 
      setDT()

  # Save parsed data to disk
  oncokb.transformed.file <- "oncokb.transformed.tsv"
  fwrite(oncokb, oncokb.transformed.file, sep="\t", quote=F)
  flog.info("    Transformed variant file saved as %s", oncokb.transformed.file)

  # Subset to only keep those GDC will use for annotations
  flog.info("  - Subsetting OncoKB variant sets")    
  oncokb.subset <- oncokb %>% 
    filter(type %in% allowedTypes) %>% 
    select(isoform, gene, alteration, protein_change, type) %>%
    setDT()

  # Save subset data to disk
  oncokb.subset.file <- "oncokb.subset.tsv"
  fwrite(oncokb.subset, oncokb.subset.file, sep="\t", quote=F)
  flog.info("    %s variants selected", nrow(oncokb.subset))
  flog.info("    Selected variant file saved as %s", oncokb.subset.file)


  # http://oncokb.org/#/gene/PAX8/alteration/PAX8-PPAR%CE%B3%20Fusion
  flog.info("  - Constructing OncoKB urls")    
  oncokb.url <- oncokb.subset %>% 
    left_join(gene.model[["transcript"]], by = c("isoform" = "gdc.isoform")) %>% 
    rename(gdc.gene.name.from.isoform = gdc.gene.name, 
           gdc.gene.id.from.isoform = gdc.gene.id) %>%
    left_join(gene.model[["gene"]], by = c("gene" = "gdc.gene.name")) %>%
    rename(gdc.gene.id.from.gene = gdc.gene.id) %>%
    rowwise() %>%
    mutate(gdc.gene.name = ifelse(is.na(gdc.gene.id.from.gene), gdc.gene.name.from.isoform, gene), 
           gdc.gene.id = ifelse(is.na(gdc.gene.id.from.gene),  gdc.gene.id.from.isoform, gdc.gene.id.from.gene),
           hgvsp = paste0("p.", protein_change),
           oncokb.url = URLencode(paste0("http://oncokb.org/#/gene/", gene, "/alteration/", alteration))) %>% 
    setDT()
 
   # Check genes that do not have a match 
  left <- oncokb.url[is.na(gdc.gene.id)]
  if (nrow(left) > 0) {
    flog.warn("Can not find gene name match for %s", left$gene)
  }

  # select columns
  oncokb.url <- oncokb.url %>%    
    select(gdc.gene.name, gdc.gene.id, gdc.transcript.id, hgvsp, type, oncokb.url) %>% 
    setDT()

  # Save final data to disk
  oncokb.final.file <- "oncokb.final.tsv"
  fwrite(oncokb.url, oncokb.final.file, sep="\t", quote=F)
  flog.info("    Final variant file saved as %s", oncokb.final.file)

  # return
  return(oncokb.url)
}



# CIVIC
# URL https://civicdb.org/events/genes/{gene_id}/summary/variants/{id}/summary
ProcessCivic <- function(allowedTypes, gene.model) {  

  civicMax <- 99999
  civic.variant.url <- paste0("https://civicdb.org/api/variants?count=", civicMax)
  flog.info("  - Downloading CIVIC raw variant JSON file")    
  civic.raw.file <- "civic.raw.json"
  download.file(civic.variant.url, civic.raw.file)
  flog.info("    Raw variant file saved as %s", civic.raw.file)  

  # extract variant record and check if numbers match
  flog.info("  - Reading CIVIC variant sets") 
  civic <- jsonlite::fromJSON(civic.raw.file)
  if (civic$`_meta`$total_count > civicMax ) {
    flog.warn("Number of variants is more than %s", civicMax)
  }
  if (civic$`_meta`$total_count != nrow(civic$records)) {
    flog.warn("Number of variants does not match total count")
  }
  civic <- civic$records

  # Flatten variant_types and extract "name"
  if(max(sapply(civic$variant_types, function(x) nrow(x))) > 2) {
    flog.error("Detected at least one variant with more than 2 variant types")
  }
  variant.type.names <- sapply(civic$variant_types, function(x) x$name)
  civic$variant.type.name <- unlist(lapply(variant.type.names, function(x) ifelse(is.null(x), NA, x[1])))
  civic$variant.type.name2 <- unlist(lapply(variant.type.names, function(x) ifelse(is.null(x), NA, x[2])))

  # Flatten coordinates and extract "representative_transcript"
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


# MyCancerGenome
# Note gene CRLF2 matches twice b/c existance of ENSG00000205755.9 and ENSGR0000205755.9 on chrX and chrY
ProcessMcg <- function(allowedTypes, gene.model) {  

  mcg.variant.url <- "https://www.mycancergenome.org/sitemap/"
  mcg.raw.file <- "mcg.raw.html"
  flog.info("  - Downloading MyCancerGenome raw variant HTML sitemap")    
  download.file(mcg.variant.url, mcg.raw.file)
  flog.info("    Raw variant file saved as %s", mcg.raw.file)

  # read html
  flog.info("  - Reading MyCancerGenome complete variant sets") 
  doc <- htmlTreeParse(mcg.raw.file, useInternal = T)

  # extract variant information
  variants <- getNodeSet(doc, "//a/h3[text() = 'Variants']/../../ul/li/a")
  mcg.complete <- data.table(name = unlist(lapply(variants, xmlValue)), 
                             url = unlist(lapply(variants, xmlAttrs)))
  # add index for easy debug and save 
  mcg.complete$index <- 1:nrow(mcg.complete)
  mcg.complete.file <- "mcg.complete.tsv"
  fwrite(mcg.complete, mcg.complete.file, sep="\t", quote=F)
  flog.info("    %s variants read", nrow(mcg.complete))
  flog.info("    Complete variant file saved as %s", mcg.complete.file)

  # extract disease and therapy
  flog.info("  - Transforming CIVIC variant sets") 
  mcg <- mcg.complete %>%
    mutate(is_association = grepl("Associated with", name),
           alteration = gsub("^(.+) (in|Associated with) .+$", "\\1", name),
           disease = ifelse(is_association, NA, gsub("^.+ in (.+)$", "\\1", name)),
           therapy = ifelse(is_association, gsub("^.+ Associated with (.+Therapy)$", "\\1", name), NA))

  # categorize, fix alias, and extract gene name
  mcg <- mcg %>% 
    mutate(type = case_when(
      grepl("Amplification$", alteration)   ~ "amplification", 
      grepl("Deletion$", alteration)        ~ "large_deletion", 
      grepl("Duplication$", alteration)     ~ "large_duplication", 
      grepl("Expression$", alteration)      ~ "expression", 
      grepl("Fusion(s)?$", alteration)      ~ "fusion", 
      grepl("Insertion$", alteration)       ~ "large_insertion", 
      grepl("Inversion$", alteration)       ~ "inversion", 
      grepl("[Mm]utation(s)?$", alteration) ~ "mutation_other", 
      grepl("Overexpression$", alteration)  ~ "overexpression", 
      grepl("Translocation$", alteration)   ~ "translocation",
      TRUE                                  ~ "other"
      )
    ) %>% 
    mutate(alteration = gsub(" (amplification|deletion|duplication|expression|fusion|fusions|insertion|inversion|mutation|mutations|overexpression|translocation)", "", alteration, ignore.case=T)) %>%
    mutate(alteration = case_when(
      alteration == "NTRK1 (TRKA)"                              ~ "NTRK1",
      alteration == "HER2 (ERBB2)"                              ~ "ERBB2",
      alteration == "ER (ESR1)"                                 ~ "ESR1",
      alteration == "PR (PGR)"                                  ~ "PGR",
      alteration == "PD-L1 (CD274)"                             ~ "CD274",
      alteration == "AR Splice Variant 7 (AR-V7)"               ~ "AR-V7",
      alteration == "HER2 (ERBB2) Exon 20"                      ~ "ERBB2 Exon 20",
      alteration == "RPN1-EVI1 (RPN1-MECOM) inv(3)(q21q26.2)"   ~ "RPN1-MECOM inv(3)(q21q26.2)",
      alteration == "RPN1-EVI1 (RPN1-MECOM) t(3;3)(q21;q26.2)"  ~ "RPN1-MECOM t(3;3)(q21;q26.2)",
      TRUE                                                      ~ alteration
      )) %>% 
    rowwise() %>%
    mutate(gene = ifelse(type=="other", NA, strsplit(alteration, " ")[[1]][1]),
           alteration = ifelse(type=="other", alteration, paste(strsplit(alteration, " ")[[1]][-1], collapse=" ")),
           alteration = gsub("^.+\\((.+)\\)$", "\\1", alteration)) %>% 
    mutate(gene = case_when(
      gene == "HER2"                  ~ "ERBB2", 
      gene == "MEK1"                  ~ "MAP2K1",
      gene == "CTNNB1(beta-catenin)"  ~ "CTNNB1", 
      TRUE                            ~ gene
      ))

  # further categorize mutation and split gene names if fusion/inversion/translocation
  mcg <- mcg %>%
    mutate(type = case_when(
      type == "mutation_other" & alteration == ""                                                                                               ~ "mutation_ns", 
      type == "mutation_other" & grepl("^[ACDEFGHIKLMNPQRSTVWY\\*]\\d+[ACDEFGHIKLMNPQRSTVWY\\*]$", alteration)                                  ~ "point_mutation", 
      type == "mutation_other" & grepl("^[ACDEFGHIKLMNPQRSTVWY]\\d+([ACDEFGHIKLMNPQRSTVWY])?fs(\\*\\d+)?$", alteration)                         ~ "frameshift",
      type == "mutation_other" & grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+dup$", alteration)                             ~ "duplication",
      type == "mutation_other" & grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+del$", alteration)                             ~ "small_deletion",
      type == "mutation_other" & grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+ins[ACDEFGHIKLMNPQRSTVWY*]+$", alteration)     ~ "insertion",
      type == "mutation_other" & grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+delins[ACDEFGHIKLMNPQRSTVWY*]+$", alteration)  ~ "delins",
      TRUE                                                                                                                                      ~ type                                                                                                                  
      )) %>% 
    mutate(type = ifelse(gene == "BCR-ABL1" & type == "point_mutation", "BCR-ABL1_point_mutation", type)) %>% 
    rowwise() %>%
    mutate(fusion.gene = ifelse(type %in% c("fusion", "inversion", "translocation", "BCR-ABL1_point_mutation") & grepl("-", gene), strsplit(gene, "-")[[1]][2], NA),
           fusion.gene = ifelse(type %in% c("fusion", "inversion", "translocation", "BCR-ABL1_point_mutation") & is.na(fusion.gene), "*", fusion.gene),
           gene = ifelse(!is.na(fusion.gene) & fusion.gene != "*", strsplit(gene, "-")[[1]][1], gene))

  # Save parsed data to disk
  mcg.transformed.file <- "mcg.transformed.tsv"
  fwrite(mcg, mcg.transformed.file, sep="\t", quote=F)
  flog.info("    Transformed variant file saved as %s", mcg.transformed.file)


  # Subset to only keep those GDC will use for annotations
  flog.info("  - Subsetting MyCancerGenome variant sets")  
  mcg.subset = mcg %>% 
    filter(type %in% allowedTypes) %>% 
    select(-is_association)

  # Save subset data to disk
  mcg.subset.file <- "mcg.subset.tsv"
  fwrite(mcg.subset, mcg.subset.file, sep="\t", quote=F)
  flog.info("    %s variants selected", nrow(mcg.subset))
  flog.info("    Selected variant file saved as %s", mcg.subset.file)

  # Add URLs
  mcg.url <- mcg.subset %>% 
    left_join(gene.model[["gene"]], by = c("gene" = "gdc.gene.name")) %>%
    mutate(gdc.gene.name = ifelse(is.na(gdc.gene.id), NA, gene), 
           hgvsp = paste0("p.", alteration),
           mcg.url = paste0("https://www.mycancergenome.org", url)) %>%
    setDT()

  # Check genes that do not have a match 
  left <- mcg.url[is.na(gdc.gene.id)]
  if (nrow(left) > 0) {
    flog.warn("Can not find gene name match for %s", left$gene)
  }

  # select columns
  mcg.url <- mcg.url %>%    
    select(gdc.gene.name, gdc.gene.id, hgvsp, type, mcg.url, name) %>% 
    setDT()
 
  # Save final data to disk
  mcg.final.file <- "mcg.final.tsv"
  fwrite(mcg.url, mcg.final.file, sep="\t", quote=F)
  flog.info("    Final variant file saved as %s", mcg.final.file)

  # return
  return(mcg.url)
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

# Get Options
opt = GetOptions()

# create and set output dir
CreateWkDir(opt$wkdir)
flog.appender(appender.file(opt$log))


# set log file
log.file <- opt$log
flog.appender(appender.file(log.file)) 
flog.info("Writing logs in %s in the working directory", log.file)

# extract gene model
flog.info("Preparing Gene Model ...")
gene.model <- GetGeneModel(opt$gtf)

flog.info("Processing OncoKB variants ...")
oncokb <- ProcessOncokb(allowedTypes, gene.model)

flog.info("Processing CIVIC variants ...")
civic <- ProcessCivic(allowedTypes, gene.model)

flog.info("Processing MyCancerGenome variants ...")
mcg <- ProcessMcg(allowedTypes, gene.model)

# Merge
flog.info("Merging databases ...")
variants <- rbind(
  oncokb %>% 
    select(gdc.gene.name, gdc.gene.id, hgvsp, oncokb.url) %>%
    rename(gene.name = oncokb.url, 
           gene.id = gdc.gene.id, 
           url = oncokb.url) %>% 
    mutate(name = paste(gdc.gene.name, gsub("^p.", "", hgvsp)), 
           source = "OncoKB"), 
  civic %>% 
    select(gdc.gene.name, gdc.gene.id, hgvsp, civic.url, name) %>%
    rename(gene.name = civic.url, 
           gene.id = gdc.gene.id, 
           url = civic.url) %>% 
    mutate(source = "CIVIC"), 
  mcg %>% 
    select(gdc.gene.name, gdc.gene.id, hgvsp, mcg.url, name) %>%
    rename(gene.name = mcg.url, 
           gene.id = gdc.gene.id, 
           url = mcg.url) %>% 
    mutate(source = "MyCancerGenome")) %>% 
  setDT()

# validate URLs
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

flog.info("Summary")
flog.info("  - %s variants in total", nrow(variants))
flog.info("  - %s unique variants", length(unique(paste(variants$gene.id, variants$hgvsp))))
flog.info("  - %s variants from OncoKB", sum(variants$source == "OncoKB"))
flog.info("  - %s variants from CIVIC", sum(variants$source == "CIVIC"))
flog.info("  - %s variants from MyCancerGenome", sum(variants$source == "MyCancerGenome"))




