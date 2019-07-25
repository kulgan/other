from collections import Counter as table

projects = ['TCGA-ACC', 
'TCGA-BLCA','TCGA-BRCA','TCGA-CESC','TCGA-CHOL','TCGA-COAD','TCGA-ESCA','TCGA-GBM','TCGA-HNSC','TCGA-KIRC','TCGA-KIRP','TCGA-LGG','TCGA-LIHC','TCGA-LUAD','TCGA-LUSC','TCGA-MESO','TCGA-PCPG','TCGA-PRAD','TCGA-SKCM','TCGA-STAD','TCGA-TGCT','TCGA-THCA','TCGA-UCEC']

surs = g.nodes(SubmittedUnalignedReads).props(experimental_strategy='ATAC-Seq').all()

rgs = []
als = []
for sur in surs:
  rgs.append(sur._SubmittedUnalignedReadsDataFromReadGroup_out[0].dst.node_id)
  als.append(sur._SubmittedUnalignedReadsDataFromReadGroup_out[0].dst._ReadGroupDerivedFromAliquot_out[0].dst.node_id)

