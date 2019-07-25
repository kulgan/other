

library(data.table)
library(dplyr)

g = fread("cmp_civic_transvar_gdcmaf_TCGA_gDNA_full_info_update_geneid_added.tsv")
c = fread("cmp_civic_transvar_gdcmaf_TCGA_cDNA_full_info_update_geneid_added.tsv")
p = fread("gdc_parsed_civic_variants_hgvs.p_aa1_format_gene_name_added.tsv")

################
# C2G
################
# keep for distinct records
c2 = c %>% 
  select(-c(aliases, db_sources)) %>% 
  distinct()

c2$nref = nchar(sapply(c2$parsed_var_hg38, function(x) strsplit(x, ":")[[1]][3])) 
c2 = c2 %>% 
  filter(nref <= 50) %>% 
  select(parsed_var_hg38, parsed_type, var_freq_in_gdc_maf, civic_var_id, gen_id, transvar_hgvs.g_hg19) %>%
  distinct() %>% 
  rename(hgvs.g = transvar_hgvs.g_hg19) %>%
  setDT()


# merge
dna = rbind(g[, -"entrez_id"], c2) %>% 
  rename(civic_gene_id = gen_id, 
         hgvs = hgvs.g) %>%  
  mutate(civic_url = paste0("https://civicdb.org/events/genes/", civic_gene_id, "/summary/variants/", civic_var_id, "/summary")) %>% 
  mutate(mapping = "DNA", 
         source = ifelse(parsed_type == "g", "gDNA", "cDNA"), 
         gdc_gene_name = NA, 
         gdc_gene_id = NA) %>% 
  select(civic_var_id, civic_gene_id, gdc_gene_name, gdc_gene_id, hgvs, parsed_var_hg38, source, mapping, civic_url) %>%
  setDT()

# recover failed hgvs.p 
protein = fread("/Users/zhenyu/github/other/gdc/clin_annot/run_20190129/civic.final.tsv")
protein$civic_var_id = sapply(protein$civic.url, function(x) strsplit(x, "\\/")[[1]][9])
protein$civic_gene_id = sapply(protein$civic.url, function(x) strsplit(x, "\\/")[[1]][6])
protein = protein[!civic_var_id %in% dna$civic_var_id] %>%
  rename(gdc_gene_name = gdc.gene.name, 
         gdc_gene_id = gdc.gene.id, 
         hgvs = hgvsp, 
         civic_url = civic.url) %>%
  mutate(mapping = "Protein", 
         source = "Protein", 
         parsed_var_hg38 = NA) %>%
  select(civic_var_id, civic_gene_id, gdc_gene_name, gdc_gene_id, hgvs, parsed_var_hg38, source, mapping, civic_url) %>%
  setDT()        

v = rbind(dna, protein)



protein = rbind(p.failed, old[, -c("gdc.transcript.id", "name")]) %>%
  mutate(mapping = "Protein") %>%
  setDT()
write.table(protein, "nam.parsed.civic.protein.tsv", col.names=T ,row.names=F, quote=F, sep="\t")



write.table(old[!civic_var_id %in% d$civic_var_id], "nam.missing.civic.tsv", col.names=T ,row.names=F, quote=F, sep="\t")


> names(g)
[1] "parsed_var_hg38"     "parsed_type"         "var_freq_in_gdc_maf"
[4] "civic_var_id"        "hgvs.g" 

> names(c2)
[1] "parsed_var_hg38"      "parsed_type"          "var_freq_in_gdc_maf" 
[4] "civic_var_id"         "transvar_hgvs.g_hg19" "civic_transcript"    
[7] "transvar_transcript"  "isoform_match"        "nref" 

