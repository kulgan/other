library(data.table)
library(jsonlite)
# library(devtools)
# library(hgvsParseR)

current = format(Sys.time(), "%Y-%m-%d_%H-%M-%S")


# OncoKB
# URL http://oncokb.org/#/gene/{Gene}/alteration/{Alteration}  
# Note: need to replace space with %20
# http://oncokb.org/api/v1/info
oncokb.variant.dump = "http://oncokb.org/api/v1/utils/allAnnotatedVariants.txt"
oncokb = fread(oncokb.variant.dump)

# make column name easy to use
setnames(oncokb, gsub(" ", "_", tolower(names(oncokb)))) 

# Fix OncoKB raw data
# Manual corrections of some protein change values
oncokb = oncokb %>% 
         mutate(alteration = ifelse(alteration == "X129R", "*129R", alteration)) %>% # fix
         mutate(protein_change = ifelse(protein_change == "X129R", "*129R", protein_change)) %>% # fix
         mutate(alteration = ifelse(alteration == "AGAP3?BRAF Fusion", "AGAP3-BRAF Fusion", alteration)) %>% # fix
         mutate(protein_change = ifelse(protein_change == "AGAP3?BRAF Fusion", "AGAP3-BRAF Fusion", protein_change)) %>%  # fix
         mutate(alteration = ifelse(alteration == "CIITA?BX648577 Fusion", "CIITA-BX648577 Fusion", alteration)) %>% # fix
         mutate(protein_change = ifelse(protein_change == "CIITA?BX648577 Fusion", "CIITA-BX648577 Fusion", protein_change)) %>% # fix
         mutate(alteration = ifelse(alteration == "FGFR2?PPHLN1 Fusion", "FGFR2-PPHLN1 Fusion", alteration)) %>% # fix
         mutate(protein_change = ifelse(protein_change == "FGFR2?PPHLN1 Fusion", "FGFR2-PPHLN1 Fusion", protein_change)) %>% # fix
         mutate(alteration = ifelse(alteration == "PAX8-PPAR? Fusion", "PAX8-PPARG Fusion", alteration)) %>% # fix
         mutate(protein_change = ifelse(protein_change == "PAX8-PPAR? Fusion", "PAX8-PPARG Fusion", protein_change)) %>% # fix
         mutate(alteration = ifelse(alteration == "ERLIN2?FGFR1 Fusion", "ERLIN2-FGFR1 Fusion", alteration)) %>% # fix
         mutate(protein_change = ifelse(protein_change == "ERLIN2?FGFR1 Fusion", "ERLIN2-FGFR1 Fusion", protein_change))

oncokb = oncokb %>% 
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
oncokb = oncokb %>% 
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

# Split gene fusions into two gene columns: gene and fusion_gene
# "*" in fusion_gene represents any partner gene 
oncokb = oncokb %>% 
         rowwise() %>%
         mutate(fusion_gene = ifelse(type == "fusion", gsub("(-$|^-)", "", gsub(gene, "", gsub(" Fusion", "", protein_change))), NA)) %>%
         mutate(fusion_gene = ifelse(type == "fusion_ns", "*", fusion_gene)) %>% 
         setDT()

# Subset to only keep those GDC will use for annotations
type.allowed = c("point_mutation", "duplication", "small_deletion", "insertion", "delins", "frameshift", "splicing")
oncokb.subset = oncokb %>% 
                filter(type %in% type.allowed) %>% 
                select(isoform, gene, alteration, protein_change, type) %>%
                setDT()

# http://oncokb.org/#/gene/PAX8/alteration/PAX8-PPAR%CE%B3%20Fusion




# CIVIC
# URL https://civicdb.org/events/genes/{gene_id}/summary/variants/{id}/summary
CIVIC_MAX = 99999
civic.variant.api = paste0("https://civicdb.org/api/variants?count=", CIVIC_MAX)

civic = jsonlite::fromJSON(civic.variant.api)

if (civic$`_meta`$total_count > CIVIC_MAX | civic$`_meta`$total_count != nrow(civic$records)) {
  print("civic not completed")
}
civic = civic$records

# Flatten variant_types and extract "name"
if(max(sapply(civic$variant_types, function(x) nrow(x))) > 2) {
  print("error")
}
variant_type_names = sapply(civic$variant_types, function(x) x$name)
civic$variant_type_name = unlist(lapply(variant_type_names, function(x) ifelse(is.null(x), NA, x[1])))
civic$variant_type_name2 = unlist(lapply(variant_type_names, function(x) ifelse(is.null(x), NA, x[2])))

# Flatten coordinates and extract "representative_transcript"
civic$transcript = civic$coordinates$representative_transcript
civic$transcript2 = civic$coordinates$representative_transcript2

# Drop columns of variant_types, coordinates and description
civic.flat = subset(civic, select=-c(variant_types, coordinates, description)) %>%
             mutate(variant_type_name = ifelse(variant_type_name == "N/A", NA, variant_type_name)) 

# Manual corrections of some alteration values
civic.flat = civic.flat %>% 
             mutate(alteration = gsub(" *\\((\\w|\\.| |_|-|>|\\+|\\?)+\\)", "", name)) %>% # removing tailing bracket comments
             mutate(alteration = gsub("^p\\.", "", alteration)) %>% # remove p. from begining             
             mutate(alteration = gsub("^([ACDEFGHIKLMNPQRSTVWY]\\d+([ACDEFGHIKLMNPQRSTVWY])?fs)Ter(\\d+)$", "\\1*\\3", alteration)) %>% # change Ter into *
             mutate(alteration = gsub("^([ACDEFGHIKLMNPQRSTVWY]\\d+([ACDEFGHIKLMNPQRSTVWY])?)FS((\\*\\d+)?)$", "\\1fs\\3", alteration)) %>% # change FS into fs
             mutate(alteration = gsub("^(([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+)DUP$", "\\1dup", alteration)) %>% # change DUP into dup
             mutate(alteration = gsub("^(([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+)DEL$", "\\1del", alteration)) %>% # change DEL into del
             mutate(alteration = gsub("^(([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+)INS([ACDEFGHIKLMNPQRSTVWY*]+)$", "\\1ins\\3", alteration)) %>% # change INS into ins
             mutate(alteration = gsub("^(([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+)DELins([ACDEFGHIKLMNPQRSTVWY*]+)$", "\\1delins\\3", alteration)) %>% # change DELins into delins
             mutate(alteration = gsub("^([ACDEFGHIKLMNPQRSTVWY\\*]\\d+)X$", "\\1*", alteration)) %>% # change terminal X into *
             mutate(alteration = case_when(
                                           alteration == "Asn67fs"    ~ "N67fs", 
                                           alteration == "Glu34Lys"   ~ "E34K", 
                                           alteration == "EZH2 Y641F" ~ "Y641F", 
                                           alteration == "MLL-MLLT3"  ~ "KMT2A-MLLT3", 
                                           alteration == "BCR-ABL"    ~ "BCR-ABL1", 
                                           TRUE                       ~ alteration
                                 )
                    ) 

# Categorize variants by variant types
civic.flat = civic.flat %>%
             mutate(type = case_when(
                                     grepl("^[ACDEFGHIKLMNPQRSTVWY\\*]\\d+[ACDEFGHIKLMNPQRSTVWY\\*]$", alteration)                                  ~ "point_mutation", 
                                     grepl("^\\w+ (FUSIONS|Fusions|fusions)$", alteration)                                                          ~ "fusion_ns", 
                                     grepl("^(\\w|-)+ (FUSION|Fusion|fusion)$", alteration)                                                         ~ "fusion",
                                     grepl("^(\\w|-)+$", alteration) & ( grepl("fusion", variant_type_name) | grepl("fusion", variant_type_name2))  ~ "fusion", 
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
                                     alteration == "Splicing alteration"                                                                            ~ "splicing",
                                     TRUE                                                                                                           ~ "other"
                           )
                    ) %>%                              
              mutate(alteration = gsub("^((\\w|-)+) (FUSION|Fusion|fusion)(S|s)?$", "\\1", alteration)) # remove tailing "fusions"


# Split gene fusions into two gene columns: gene and fusion_gene
# "*" in fusion_gene represents any partner gene 
civic.flat = civic.flat %>% 
             rowwise() %>%
             mutate(fusion_gene = ifelse(type == "fusion", gsub("(-$|^-)", "", gsub(entrez_name, "", alteration)), NA)) %>%
             mutate(fusion_gene = ifelse(type == "fusion_ns", "*", fusion_gene)) %>% 
             setDT()

# Subset to only keep those GDC will use for annotations
# type.allowed = c("point_mutation", "small_deletion", "insertion", "delins", "frameshift", "splicing")
civic.subset = civic.flat %>% 
               filter(type %in% type.allowed) %>% 
               select(id, name, entrez_name, gene_id, transcript, transcript2, alteration, type) %>%
               setDT()




# MyCancerGenome

library(XML)
library(RCurl)

mcg.sitemap.url = "https://www.mycancergenome.org/sitemap/"

# read html
doc = htmlTreeParse(getURL(mcg.sitemap.url), useInternal = T)

# extract variant information
variants = getNodeSet(doc, "//a/h3[text() = 'Variants']/../../ul/li/a")
mcg = data.table(name = unlist(lapply(variants, xmlValue)), 
                 url = unlist(lapply(variants, xmlAttrs)))

# extract disease and therapy
mcg = mcg %>%
      mutate(is_association = grepl("Associated with", name)) %>%
      mutate(alteration = gsub("^(.+) (in|Associated with) .+$", "\\1", name)) %>%
      mutate(disease = ifelse(is_association, NA, gsub("^.+ in (.+)$", "\\1", name))) %>%
      mutate(therapy = ifelse(is_association, gsub("^.+ Associated with (.+Therapy)$", "\\1", name), NA)) %>% 
      setDT()

# categorize, fix alias, and extract gene name
mcg = mcg %>% 
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
                                    alteration=="NTRK1 (TRKA)"                              ~ "NTRK1",
                                    alteration=="HER2 (ERBB2)"                              ~ "ERBB2",
                                    alteration=="ER (ESR1)"                                 ~ "ESR1",
                                    alteration=="PR (PGR)"                                  ~ "PGR",
                                    alteration=="PD-L1 (CD274)"                             ~ "CD274",
                                    alteration=="AR Splice Variant 7 (AR-V7)"               ~ "AR-V7",
                                    alteration=="HER2 (ERBB2) Exon 20"                      ~ "ERBB2 Exon 20",
                                    alteration=="RPN1-EVI1 (RPN1-MECOM) inv(3)(q21q26.2)"   ~ "RPN1-MECOM inv(3)(q21q26.2)",
                                    alteration=="RPN1-EVI1 (RPN1-MECOM) t(3;3)(q21;q26.2)"  ~ "RPN1-MECOM t(3;3)(q21;q26.2)",
                                    TRUE ~ alteration
                                    )
            ) %>% 
      rowwise() %>%
      mutate(gene = ifelse(type=="other", NA, strsplit(alteration, " ")[[1]][1])) %>%
      mutate(alteration = ifelse(type=="other", alteration, paste(strsplit(alteration, " ")[[1]][-1], collapse=" "))) %>%
      mutate(alteration = gsub("^.+\\((.+)\\)$", "\\1", alteration)) 

# further categorize mutation and split gene names if fusion/inversion/translocation
mcg = mcg %>%
      mutate(type = case_when(
                              type=="mutation_other" & alteration=="" ~ "mutation_ns", 
                              type=="mutation_other" & grepl("^[ACDEFGHIKLMNPQRSTVWY\\*]\\d+[ACDEFGHIKLMNPQRSTVWY\\*]$", alteration)                                  ~ "point_mutation", 
                              type=="mutation_other" & grepl("^[ACDEFGHIKLMNPQRSTVWY]\\d+([ACDEFGHIKLMNPQRSTVWY])?fs(\\*\\d+)?$", alteration)                         ~ "frameshift",
                              type=="mutation_other" & grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+dup$", alteration)                             ~ "duplication",
                              type=="mutation_other" & grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+del$", alteration)                             ~ "small_deletion",
                              type=="mutation_other" & grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+ins[ACDEFGHIKLMNPQRSTVWY*]+$", alteration)     ~ "insertion",
                              type=="mutation_other" & grepl("^([ACDEFGHIKLMNPQRSTVWY]\\d+_)?[ACDEFGHIKLMNPQRSTVWY]\\d+delins[ACDEFGHIKLMNPQRSTVWY*]+$", alteration)  ~ "delins",
                              TRUE                                                                                                                                    ~ type                                                                                                                  
                              )
            ) %>% 
      rowwise() %>%
      mutate(gene2 = ifelse(type %in% c("fusion", "inversion", "translocation") & grepl("-", gene), strsplit(gene, "-")[[1]][2], NA)) %>%
      mutate(gene2 = ifelse(type %in% c("fusion", "inversion", "translocation") & is.na(gene2), "*", gene2)) %>%
      mutate(gene = ifelse(!is.na(gene2) & gene2 != "*", strsplit(gene, "-")[[1]][1], gene))

# Subset to only keep those GDC will use for annotations
mcg.subset = mcg %>% 
             filter(type %in% type.allowed) %>% 
             select(-is_association) %>%
             setDT()













