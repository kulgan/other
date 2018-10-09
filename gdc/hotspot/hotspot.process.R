# The script combines cancer hotspot mutations from the following 4 sources
# 1. Hotspots v1 from Chang et al. 2016 (dx.doi.org/10.1038/nbt.3391)
# 2. Hotspots v2 from Chang et al. 2018 (dx.doi.org/10.1158/2159-8290.CD-17-0321)
# 3. Somatic mutations associated with clonal hematopoietic expansion from Xie et al. (dx.doi.org/10.1038/nm.3733)
# 4. Somatically mutable germline sites at MSH6:F1088, TP53:R290, TERT:E280, ASXL1:G645_G646 (listed in GENIE data manual)
# 1/3/4 have been merged in https://github.com/mskcc/vcf2maf/blob/v1.6.14/data/known_somatic_sites.bed 


library(data.table)
library(dplyr)
library(tidyr)


setwd("/Users/zhenyu/github/other/gdc/hotspot")

########################################
# Hotspots v2
########################################

# v2 hotspot exists from both paper supplementary, and http://www.cancerhotspots.org/
# unfortuantely, both sources contain unique information, as such comparison and merge might be necessary

hotspot.snp.v2.paper <- fread("hotspot.v2.from_website.csv")
hotspot.indel.v2.paper <- fread("hotspot.v2.recovered.from_website.csv")

# read and extract variants from v2 paper
hotspot.snp.v2.paper.reformat <- hotspot.snp.v2.paper %>% 
  select(Hugo_Symbol, Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid) %>% 
  rowwise() %>% 
  mutate(splice = grepl("splice", Amino_Acid_Position), 
         Reference_Amino_Acid = strsplit(Reference_Amino_Acid, ":")[[1]][1], 
         Variant_Amino_Acid = strsplit(Variant_Amino_Acid, ":")[[1]][1], 
         change = ifelse(splice, Amino_Acid_Position, paste0(Reference_Amino_Acid, Amino_Acid_Position, Variant_Amino_Acid)), 
         type = ifelse(splice, "splice", "snv")) %>% 
  select(Hugo_Symbol, change, type) 

hotspot.indel.v2.paper.reformat <- hotspot.indel.v2.paper %>%
  rowwise() %>%
  mutate(change = strsplit(Variant_Amino_Acid, ":")[[1]][1], 
         type = "indel") %>%
  select(Hugo_Symbol, change, type)

hotspot.v2 <- rbind(hotspot.snp.v2.paper.reformat, hotspot.indel.v2.paper.reformat) %>%
  setDT()




# we may want to compare variants in paper to that in the website
check.web = F

if(check.web) {
  # read and extract variants from v2 website
  hotspot.snp.v2.web <- fread("hotspot.snv.v2.from_paper.csv")
  hotspot.indel.v2.web <- fread("hotspot.indel.v2.from_paper.csv")

  hotspot.snp.v2.web.reformat <- hotspot.snp.v2.web %>% 
    rename(Hugo_Symbol = "Hugo Symbol", 
           Reference_Amino_Acid = "Reference Amino Acid", 
           Amino_Acid_Position = "Amino Acid Position", 
           Variant_Amino_Acid = "Variant Amino Acid") %>% 
    select(Hugo_Symbol, Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid) %>%
    separate_rows(., Variant_Amino_Acid, sep="\\|") %>% 
    rowwise() %>% 
    mutate(Variant_Amino_Acid = strsplit(Variant_Amino_Acid, ":")[[1]][1], 
           change = ifelse(Variant_Amino_Acid == "splice", Amino_Acid_Position, paste0(Reference_Amino_Acid, Amino_Acid_Position, Variant_Amino_Acid)), 
           type = ifelse(Variant_Amino_Acid == "splice", "splice", "snv")) %>% 
    select(Hugo_Symbol, change, type)

  hotspot.indel.v2.web.reformat <- hotspot.indel.v2.web %>% 
    rename(Hugo_Symbol = "Hugo Symbol", 
           Variant_Amino_Acid = "Variant Amino Acid") %>% 
    select(Hugo_Symbol, Variant_Amino_Acid) %>%
    separate_rows(., Variant_Amino_Acid, sep="\\|") %>% 
    rowwise() %>% 
    mutate(change = strsplit(Variant_Amino_Acid, ":")[[1]][1], 
           type = "indel") %>%
    select(Hugo_Symbol, change, type)

  # confirming website and paper have the same variants
  hotspot.v2.web <- rbind(hotspot.snp.v2.web.reformat, hotspot.indel.v2.web.reformat) %>% 
    setDT()

  # compare to confirm these two sets are the same
  anti_join(hotspot.v2.web, hotspot.v2)
  anti_join(hotspot.v2, hotspot.v2.web)
}


# we may want to extract coordinates, that requires merging snp variants from web and indel variants from paper
get.coord = F

if(get.coord & check.web) {
  hotspot.indel.v2.paper.reformat2 <- hotspot.indel.v2.paper %>%
    rowwise() %>%
    mutate(change = strsplit(Variant_Amino_Acid, ":")[[1]][1], 
           type = "indel", 
           chr = strsplit(Genomic_Position, ":")[[1]][1], 
           start = strsplit(Genomic_Position, ":")[[1]][2], 
           start = strsplit(start, "_")[[1]][1]) %>%
    select(Hugo_Symbol, chr, start)

  hotspot.snp.v2.web.reformat2 <- hotspot.snp.v2.web %>% 
    rename(Hugo_Symbol = "Hugo Symbol", 
           Reference_Amino_Acid = "Reference Amino Acid", 
           Amino_Acid_Position = "Amino Acid Position", 
           Variant_Amino_Acid = "Variant Amino Acid", 
           Genomic_Position = "Genomic Position") %>% 
    select(Hugo_Symbol, Amino_Acid_Position, Reference_Amino_Acid, Genomic_Position) %>%
    separate_rows(., Genomic_Position, sep="\\|") %>% 
    rowwise() %>% 
    mutate(chr = strsplit(Genomic_Position, ":")[[1]][1], 
           start = strsplit(Genomic_Position, ":")[[1]][2], 
           start = strsplit(start, "_")[[1]][1]) %>%
    select(Hugo_Symbol, chr, start)

  hotspot.v2.with.coord <- rbind(hotspot.indel.v2.paper.reformat2, hotspot.snp.v2.web.reformat2) %>% 
    mutate(start = as.integer(start),
           end = start) %>% 
    setDT() %>% 
    setkey(chr, start, end)
}



########################################
# Hotspots v1
########################################

# Build variant list from the start
# https://raw.githubusercontent.com/taylor-lab/hotspots/master/publication_hotspots.tsv
hotspot.v1 <- fread("publication_hotspots.tsv")

hotspot.v1.reformat <- hotspot.v1 %>% 
  rename(Hugo_Symbol = "Hugo Symbol", 
         Variant_Amino_Acid = "Variant Amino Acid") %>% 
  select(Hugo_Symbol,Codon, Variant_Amino_Acid) %>%
  separate_rows(., Variant_Amino_Acid, sep="\\|") %>% 
  rowwise() %>% 
  mutate(Variant_Amino_Acid = strsplit(Variant_Amino_Acid, ":")[[1]][1], 
         splice = grepl("splice", Codon), 
         Codon = ifelse(splice == T, paste0("X", substr(Codon, 2, nchar(Codon))), Codon),
         change = paste0(Codon, Variant_Amino_Acid),
         type = ifelse(splice == T, "splice", "snv")) %>% 
  select(Hugo_Symbol, change, type) %>% 
  setDT()


########################################
# hematopoietic expansion
########################################

hem <- fread("hematopoietic.nm.3733-S9.txt")

hem <- hem %>% 
  rename(Hugo_Symbol = "Gene") %>%
  mutate(first.char = substr(Annotation, 1, 1), 
         last2.char = substr(Annotation, nchar(Annotation)-1, nchar(Annotation)),
         type = ifelse(first.char == "e", "splice_region", "snv"), 
         type = ifelse(last2.char == "fs", "frameshift", type), 
         change = gsub("^p\\.", "", Annotation)) %>%
  select(Hugo_Symbol, change, type) %>%
  distinct() %>%
  setDT()

########################################
# Other Somatic Mutations indicated by GENIE
########################################
# the list is MSH6:F1088, TP53:R290, TERT:E280, ASXL1:G645_G646
# the "fs" status is obtained from MSK known_somatic_sites.bed 
# TERT:E280K is a freq mutation in TCGA
other <- data.table(Hugo_Symbol = c("MSH6", "TP53", "TERT", "ASXL1", "ASXL1"),
                    change = c("F1088fs", "R290fs", "E280K", "G645fs", "G646fs"), 
                    type = c("frameshift", "frameshift", "snv", "frameshift", "frameshift"))



########################################
# Merge
########################################
hotspot <- rbind(hotspot.v2, hotspot.v1.reformat, hem, other) %>% 
  distinct() %>%
  setDT()


########################################
# Gene Annotation
########################################
gencode <- fread("gencode.v22.genes.txt")

# find gene symbol not in data model
print(unique(hotspot[!Hugo_Symbol %in% gencode$gene_name]$Hugo_Symbol))
# [1] "WDR52"    "LST3"     "C16orf80"

# Swap gene names. Since the Hugo_Symbols are not unique, I am not going to get Ensemble IDs
mapping <- data.table(old = c("WDR52", "LST3", "C16orf80"), 
                     new = c("CFAP44", "SLCO1B3", "CFAP20"))
hotspot <- hotspot %>% 
  mutate(Hugo_Symbol = ifelse(Hugo_Symbol == "WDR52", "CFAP44", Hugo_Symbol), 
         Hugo_Symbol = ifelse(Hugo_Symbol == "LST3", "SLCO1B3", Hugo_Symbol), 
         Hugo_Symbol = ifelse(Hugo_Symbol == "C16orf80", "CFAP20", Hugo_Symbol)) %>%
  arrange(type, Hugo_Symbol) %>% 
  setDT()

write.table(hotspot, "GDC.mutation.hotspot.20181009.tsv", col.names=T, row.names=F, sep="\t", quote=F)





