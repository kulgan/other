library(data.table)
library(dplyr)
library(tidyr)


setwd("~/github/other/gdc/clin_annot/test/")
# input
# https://portal.gdc.cancer.gov/files/03652df4-6090-4f5a-a2ff-ee28a37f9301
maf.file <- "gunzip -c TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.gz"
dna.file <- "civic_gdcmaf_mapping_dna.tsv"
prot.file <- "civic_gdcmaf_mapping_prot.tsv"

# read (1389 variants in maf)
maf <- fread(maf.file)
dna <- fread(dna.file)
prot <- fread(prot.file)

# function to accept MAF and DNA mapping file, and generate an index mapping 
AnnotateDNA <- function(maf, dna.map) {
  # add index number for easier QC
  maf$maf_index <- 1:nrow(maf)
  dna$dna_index <- 1:nrow(dna)

  # Mapping of DNA Changes
  # 4 keys between bio DNA mapping file and GDC MAFs
  # - chromosome → Chromosome
  # - start_position → Start_Position
  # - reference_allele → Reference_Allele
  # - alternative_allele → Allele
  # 2 keys between bio DNA mapping file and CIVIC database
  # - civic_gene_id
  # - civic_var_id

  # rename columns
  dna <- dna %>%
    rename( Chromosome = chromosome, 
            Start_Position = start_position, 
            Reference_Allele = reference_allele, 
            Allele = alternative_allele) %>%
    select(-source)

  # join to get mapping
  map <- maf %>% 
    inner_join(dna, by = c("Chromosome", "Start_Position", "Reference_Allele", "Allele")) %>% 
    select(maf_index, dna_index, civic_var_id, civic_gene_id) %>% 
    distinct() %>%
    setDT()

  # return   
  return(map)
}


# function to accept MAF and Protein mapping file, and generate an index mapping 
AnnotateProtein <- function(maf, prot.map) {
  # add index number for easier QC
  maf$maf_index <- 1:nrow(maf)
  prot.map$prot_index <- 1:nrow(prot.map)

  # Mapping of Protein Changes
  # 2 keys between bio Protein mapping file and all_effects in GDC MAFs. 
  # - hugo_symbol → Hugo_Symbol
  # - hgvsp → hgvs.p_short
  # - The above two keys should match together. For example, (BRAF, p.V600E) should NOT match to BRAF,missense_variant,p.V601D,ENST00000546424,,c.*965T>C,,,,,-1;TP53,missense_variant,p.V600E,ENST00000546425,,c.*966T>C,,,,,-1
  # 2 keys between bio Protein mapping file and CIVIC database
  # - civic_gene_id
  # - civic_var_id

  # parsing hugo_symbol and hgvsp out of all_effects
  maf <- maf %>% 
    select(all_effects, maf_index) %>%
    mutate(all_effects = strsplit(all_effects, ";")) %>%
    unnest(all_effects) %>%
    setDT()

  maf$hugo_symbol <- sapply(maf$all_effects, function(x) strsplit(x, ",")[[1]][1])
  maf$hgvsp <- sapply(maf$all_effects, function(x) strsplit(x, ",")[[1]][3])

  prot.map <- prot.map %>% 
    rename(hgvsp = hgvs.p)

  # join to get mapping
  map <- maf %>% 
    inner_join(prot.map, by = c("hugo_symbol", "hgvsp")) %>% 
    select(maf_index, prot_index, civic_var_id, civic_gene_id) %>%
    distinct() %>%
    setDT()

  # return   
  return(map)
}

# get mappings
dna.mapping <- AnnotateDNA(maf, dna)
prot.mapping <- AnnotateProtein(maf, prot)

# check for duplicates
if (sum(duplicated(dna.mapping$maf_index)) > 0 | sum(duplicated(prot.mapping$maf_index)) > 0) {
  print("the same variants could have been annotated twice")
}

# Merge dna.mapping and protein.mapping
# We will trust DNA mapping more 
mapping <- merge(dna.mapping, prot.mapping, by = "maf_index", all = T, suffixes = c("", "_prot"))

# QC checks
if(nrow(mapping[!is.na(dna_index) & !is.na(prot_index)]) > 0) {
  print("the same variants have been annotated by both DNA and Protein Mapping")
} 




