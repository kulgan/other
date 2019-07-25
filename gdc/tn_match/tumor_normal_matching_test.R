library(data.table)
library(dplyr)

setwd("~/Downloads")

both = fread("CPTAC-3-SARS_batches.tsv", na.strings = "None")
both = both %>% 
  mutate(used_variant_calling = ifelse(used_variant_calling == "Yes", T, F)) %>%
  rename(gdc.no_matched_normal = no_matched_normal_wxs, 
         gdc.tissue_type = tissue_type, 
         gdc.selected_normal = selected_normal_wxs) %>% 
  setDT()

both$bio.tissue_type = getSampleTissueType(tissue_type = both$gdc.tissue_type, sample_type = both$sample_type) 

both = both %>% 
  mutate(existing.bio.selected_normal = ifelse(sar_batch == 1 & bio.tissue_type == "Normal", T, NA)) %>%
  setDT()

both$case_has_blood_cancer = both$case_id %in% getCaseWithBloodCancer(case_id = both$case_id, sample_type = both$sample_type) 


meta = both

case.to_be_ignored = meta[existing.bio.selected_normal == T | gdc.no_matched_normal == T]$case_id


  meta = meta %>% 
    mutate(to_be_ignored = case_id %in% case.to_be_ignored) %>% 
    setDT()

  meta.to_be_ignored = meta[to_be_ignored == T]
  meta.to_be_ignored$bio.selected_normal = meta.to_be_ignored$existing.bio.selected_normal
  meta = meta[to_be_ignored == F]

  case.normal.summary = meta %>% 
    filter(bio.tissue_type == "Normal") %>% 
    group_by(case_id) %>% 
    summarise(normal_count = n()) %>% 
    setDT()


  case.unambiguous = unique(case.normal.summary[normal_count == 1]$case_id)

  meta.unambiguous = meta[case_id %in% case.unambiguous]
  meta.ambiguous = meta[!case_id %in% case.unambiguous]

  meta.unambiguous = meta.unambiguous %>% 
    mutate(bio.selected_normal = ifelse(bio.tissue_type == "Normal", T, NA))

solid.order = c("Blood Derived Normal", 
                "Bone Marrow Normal",
                "Mononuclear Cells from Bone Marrow Normal",
                "Fibroblasts from Bone Marrow Normal",
                "Lymphoid Normal",
                "Buccal Cell Normal",
                "Solid Tissue Normal",
                "EBV Immortalized Normal")

blood.order = c("Solid Tissue Normal",
                "Buccal Cell Normal",
                "Lymphoid Normal",
                "Fibroblasts from Bone Marrow Normal",
                "Mononuclear Cells from Bone Marrow Normal",
                "Bone Marrow Normal",
                "Blood Derived Normal",
                "EBV Immortalized Normal")

  meta.ambiguous.solid = meta.ambiguous[case_has_blood_cancer == F & bio.tissue_type == "Normal"] 
  meta.ambiguous.blood = meta.ambiguous[case_has_blood_cancer == T & bio.tissue_type == "Normal"] 

selected.normal = meta.ambiguous.solid[sample_type == "Blood Derived Normal"]$aliquot_id

  meta.ambiguous = meta.ambiguous %>% 
    mutate(bio.selected_normal = ifelse(aliquot_id %in% selected.normal, T, NA))

temp = rbind(meta.to_be_ignored, meta.ambiguous, meta.unambiguous)



temp = meta.ambiguous.solid %>% 
  mutate(sample_type = factor(sample_type, levels = solid.order)) %>% 
  arrange(case_id, sample_type) %>% 
  group_by(case_id) %>% 
  slice(1) 


    left_join(case.normal.summary, by = "case_id") 




getSelectedNormal = function() {

  # If selected_normals per strategy per case exist in bio graph, keep as is and ignore what's in the gdc graph (skip to the end)
  meta = meta %>% 
    mutate(bio.selected_normal = ifelse(existing.bio.selected_normal == T, T, NA)) %>% 
    setDT()

# If no_matched_normal per strategy per case exist in bio graph, keep as is and ignore what's in the gdc graph (skip to the end)
  case.to_be_ignored = meta[existing.bio.selected_normal == T | no_matched_normal == T]$case_id

  meta = meta %>% 
    mutate(to_be_ignored = case_id %in% case.to_be_ignored) %>% 
    setDT()

  meta.to_be_ignored = meta[to_be_ignored == T]
  meta = meta[to_be_ignored == F]

# For submitter populated selected_normals, remove any whose bio.aliquot.tissue_type is not "Normal" (see TT-656 about how to populate bio.aliquot.tissue_type)
  meta %>% 
    mutate(gdc.selected_normal = ifelse(bio.tissue_type == "Normal", gdc.selected_normal, F)) %>% 
    setDT()

# If there is only one normal aliquot per strategy per case in bio graph, it automatically becomes the selected_normal of that strategy of the case
# if only one selected_normal per strategy per case exists in gdc graph: if more than one, pick only one to populate bio graph (see how to set priority below)
# If no user defined selected normal detected from gdc graph, and no pre-exist selected_normal in bio graph, and there are more than one normal per strategy per case, pick only one to populate bio graph (see how to set priority below)

  case.normal.summary = meta %>% 
    filter(bio.tissue_type == "Normal") %>% 
    group_by(case_id) %>% 
    summarise(normal_count = n(), 
              selected_normal_count = sum(gdc.selected_normal, na.rm = T)) %>% 
    setDT()

  case.unambiguous = unique(case.normal.summary[normal_count == 1 | selected_normal_count == 1]$case_id)

  meta.unambiguous = meta[case_id %in% case.unambiguous]
  meta.unambiguous %>% 
    mutate(bio.selected_normal = ifelse(bio.tissue_type == "Normal", T, bio.selected_normal))

  meta.ambiguous = meta[!case_id %in% case.unambiguous]


    left_join(case.normal.summary, by = "case_id") 

     %>% 
    mutate(bio.selected_normal = ifelse((normal_count == 1 | selected_normal_count == 1) & bio.tissue_type == "Normal", T, bio.selected_normal)) %>% 


}

getCaseWithBloodCancer = function(case_id, sample_type) {
  meta = data.table(case_id, sample_type)
  # if aliquot is a blood cancer sample
  meta$aliquot_has_blood_cancer = meta$sample_type %in% c("Blood Derived Cancer - Bone Marrow, Post-treatment", 
                                                          "Blood Derived Cancer - Peripheral Blood, Post-treatment", 
                                                          "Primary Blood Derived Cancer - Peripheral Blood",
                                                          "Recurrent Blood Derived Cancer - Peripheral Blood",
                                                          "Primary Blood Derived Cancer - Bone Marrow", 
                                                          "Recurrent Blood Derived Cancer - Bone Marrow")
  # summarise to case
  case.summary = meta %>% 
    group_by(case_id) %>% 
    summarise(case_has_blood_cancer = sum(aliquot_has_blood_cancer) > 0)

  # return cases with blood cancer
  return(case.summary$case_id[which(case.summary$case_has_blood_cancer)])
}


getSampleTissueType = function(tissue_type, sample_type) {
  meta = data.table(tissue_type, sample_type)
  meta$bio.tissue_type = "Unknown"
 
  known_index = which(meta$tissue_type %in% c("Tumor", "Normal", "Abnormal", "Peritumoral"))
  meta$bio.tissue_type[known_index] = meta$tissue_type[known_index]

  index = which(meta$sample_type == "Tumor Adjacent Normal - Post Neo-adjuvant Therapy")
  index = setdiff(index, known_index)
  meta$bio.tissue_type[index] = "Peritumoral"


  index = which(meta$sample_type %in% c("Additional Metastatic",
                                        "Additional - New Primary",
                                        "Blood Derived Cancer - Bone Marrow, Post-treatment",
                                        "Blood Derived Cancer - Peripheral Blood, Post-treatment",
                                        "FFPE Recurrent",
                                        "Human Tumor Original Cells",
                                        "Metastatic",
                                        "Primary Blood Derived Cancer - Peripheral Blood",
                                        "Recurrent Blood Derived Cancer - Peripheral Blood",
                                        "Primary Blood Derived Cancer - Bone Marrow",
                                        "Primary Tumor",
                                        "Primary Xenograft Tissue",
                                        "Recurrent Blood Derived Cancer - Bone Marrow",
                                        "Recurrent Tumor",
                                        "Tumor"))
  index = setdiff(index, known_index)
  meta$bio.tissue_type[index] = "Tumor"

  index = which(meta$sample_type %in% c("Blood Derived Normal",
                                        "Bone Marrow Normal",
                                        "Buccal Cell Normal",
                                        "EBV Immortalized Normal",
                                        "Fibroblasts from Bone Marrow Normal",
                                        "Lymphoid Normal",
                                        "Mononuclear Cells from Bone Marrow Normal",
                                        "Solid Tissue Normal"))
  index = setdiff(index, known_index)
  meta$bio.tissue_type[index] = "Normal"

  return(meta$bio.tissue_type)
}

