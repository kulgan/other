library(data.table)
library(dplyr)

setwd("~/Downloads/target/")

# read AR and VCF uuids from API portal query
portal = fread("TARGET_WXS_UUIDS.txt", h=T)

# read sar
sar = fread("target.wxs.sar.tsv", h=T, na.strings="")
names(sar)[1] = "sar.idx" 
sar$sar.idx = sar$sar.idx + 1

# read ar
ar = fread("target.wxs.ar.tsv", h=T, na.strings="")
names(ar)[1] = "ar.idx" 
ar$ar.idx = ar$ar.idx + 1
ar$ar_in_portal = ar$ar %in% portal$id

# read vcf
vcf = fread("target.wxs.vcf.tsv", h=T, na.strings="")
names(vcf)[1] = "vcf.idx" 
vcf$vcf.idx = vcf$vcf.idx + 1
vcf$vcf_in_portal = vcf$vcf %in% portal$id


# merge ar and sar
bam = merge(ar, sar, all=T)
bam$tissue = "normal"
bam$tissue[which(bam$sample.sample_type_id %in% 1:9)] = "tumor"

# assume vcf.state == "submitted" old VarScan2 VCFs are removed
vcf = vcf[vcf.state == "released"]

# focus on existing ARs
bam2 = bam[!is.na(ar)]


vcf = vcf %>% 
  mutate(tumor_ar = ifelse(ar1 %in% bam2[tissue=="tumor"]$ar, ar1, ar2), 
         normal_ar = ifelse(ar1 %in% bam2[tissue=="normal"]$ar, ar1, ar2))%>% 
  setDT()

m = match(vcf$tumor_ar, bam$ar)
vcf$tumor_aliquot = bam$aliquot[m]
vcf$tumor_barcode = bam$aliquot.submitter_id
vcf$case = bam$case[m]


  
ar$ar.index = gsub(".+AlignedReadsIndexDerivedFromAlignedReads\\(\\(((\\w|-)+)\\).+", "\\1", ar$ar.index)
ar$ar.sar = gsub(".*\\(((\\w|-)+)\\)>]", "\\1", ar$ar.sar)
ar$ssm.wf = lapply(sapply(ar$ssm.wf, function(x) strsplit(x, ",")[[1]]), function(x) gsub(".+SomaticMutationCallingWorkflowPerformedOnAlignedReads\\(\\(((\\w|-)+)\\)-.+", "\\1", x))
ar$ssm.wf[which(ar$ssm.wf == "[]")] = NA
ar = ar[, -"V1"]
ar$ar_in_portal = ar$ar %in% portal$id

# read vcf
vcf = fread("tcga.raw_vcf.tsv", h=T)
vcf$ssm.wf = gsub(".+\\(((\\w|-)+)\\)>\\]", "\\1", vcf$ssm.wf)
vcf$vcf.index = gsub(".+SomaticMutationIndexDerivedFromSimpleSomaticMutation\\(\\(((\\w|-)+)\\).+", "\\1", vcf$vcf.index)
vcf$vcf.index[which(vcf$vcf.index == "[]")] = NA
vcf = vcf[, -"V1"]
vcf$vcf_in_portal = vcf$vcf %in% portal$id

# read sample
sample = fread("sample.tsv", h=T)
sample = sample[, -"V1"]
sample$case = gsub("\\)>\\]", "", gsub("\\[<\\Case\\(", "", sample$case))

# match AR and SAR, and found AR "58d4d9ed-2aae-45b9-b7df-7732987fdc25" does not link to any SAR
w = which(ar$ar.sar == "[]")
ar$comment = NA
if(ar$ar[w] == "58d4d9ed-2aae-45b9-b7df-7732987fdc25"){
  ar$ar.sar[w] = "a9f011aa-0513-42ce-90cb-e8ba48bb3698"
  ar$comment[w] = "missing_ar_to_sar_edge"
}

# remove the previous AR for now
ar = ar[-w]

# merge AR and SAR
d = merge(sar, ar, by.x="sar", by.y="ar.sar", all=T)
# identify new data (188 TGCT WXS)
d[is.na(ar)]$comment = "sar without ar"






wf fb39f86f-e96c-4f4e-bcf2-925ee945f434



vcf 205861a9-01c2-47a0-a3b3-52fb07c61bf3 does not have batch id


In [80]: ar = g.nodes(AlignedReads).get('58d4d9ed-2aae-45b9-b7df-7732987fdc25')
    ...:

In [81]: ar._AlignedReadsDataFromAlignmentCocleaningWorkflow_out
    ...:
Out[81]: [<AlignedReadsDataFromAlignmentCocleaningWorkflow((58d4d9ed-2aae-45b9-b7df-7732987fdc25)-[data_from]->(fb39f86f-e96c-4f4e-bcf2-925ee945f434)>]

In [82]: ar._AlignedReadsMatchedToSubmittedAlignedReads_out
    ...:
Out[82]: []

In [83]: wf = g.nodes(AlignmentCocleaningWorkflow).get('fb39f86f-e96c-4f4e-bcf2-925ee945f434')
    ...:

In [84]: wf._AlignmentCocleaningWorkflowPerformedOnSubmittedAlignedReads_out
Out[84]: [<AlignmentCocleaningWorkflowPerformedOnSubmittedAlignedReads((fb39f86f-e96c-4f4e-bcf2-925ee945f434)-[performed_on]->(a9f011aa-0513-42ce-90cb-e8ba48bb3698)>]

In [85]: sar = g.nodes(SubmittedAlignedReads).get('a9f011aa-0513-42ce-90cb-e8ba48bb3698')
    ...:

In [86]: sar.aligned_reads_files
Out[86]: [<AlignedReads(66832985-f999-4020-8f55-50853e3d86dd)>]

In [87]: sar._AlignedReadsMatchedToSubmittedAlignedReads_in
Out[87]: [<AlignedReadsMatchedToSubmittedAlignedReads((66832985-f999-4020-8f55-50853e3d86dd)-[matched_to]->(a9f011aa-0513-42ce-90cb-e8ba48bb3698)>]


['_AlignedReadsIndexDerivedFromSubmittedAlignedReads_in',
 '_AlignedReadsMatchedToSubmittedAlignedReads_in',
 '_AlignmentCocleaningWorkflowPerformedOnSubmittedAlignedReads_in',
 '_AlignmentWorkflowPerformedOnSubmittedAlignedReads_in',
 '_AnalysisMetadataDerivedFromSubmittedAlignedReads_in',
 '_AnnotationAnnotatesSubmittedAlignedReads_in',
 '_ReadGroupQcDataFromSubmittedAlignedReads_in',
 '_RnaExpressionWorkflowPerformedOnSubmittedAlignedReads_in',
 '_SubmittedAlignedReadsDataFromReadGroup_out',
 '_SubmittedAlignedReadsRelatesToCase_out',






