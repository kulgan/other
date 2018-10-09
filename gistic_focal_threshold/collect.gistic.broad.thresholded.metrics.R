library(data.table)

setwd("~/SCRATCH/gistic/data/all_thresholded.by_genes")


####################################################################################################################
# convert focal_data_by_genes to scores, and only keep protein-coding genes in chr1-22/X
####################################################################################################################
diseases = fread("../d", h=F)$V1

temp = fread("ACC.all_thresholded.by_genes.txt")
ngene = nrow(temp) # 59852

gene.summary = data.table( gene = temp$"Gene Symbol", 
                           count = rep(0, ngene))

aliquot.summary = data.table( aliquot = character(0), 
                              count = numeric(0))

project.summary = data.table(project = diseases, del = 0, loss = 0, neutral = 0, gain = 0, amp = 0, aliquot = 0)


for (i in 1:length(diseases)) {
  disease = diseases[i]
  print(disease)
  
  # read file
  broad.file = paste0("./", disease, ".all_thresholded.by_genes.txt")
  broad = fread(broad.file)
  data = broad[, 4:ncol(broad)]

  # collect metrics
  gene.summary$count = gene.summary$count + apply(data, 1, function(x) sum(x != 0))

  my.aliquot.summary = apply(data, 2, function(x) sum(x != 0))
  aliquot.summary = rbind(aliquot.summary, data.table(aliquot = names(my.aliquot.summary), count = my.aliquot.summary))

  tbl = table(unlist(data))
  project.summary$del[i] = tbl[which(names(tbl) == "-2")] 
  project.summary$loss[i] = tbl[which(names(tbl) == "-1")]
  project.summary$neutral[i] = tbl[which(names(tbl) == "0")]
  project.summary$gain[i] = tbl[which(names(tbl) == "1")]
  project.summary$amp[i] = tbl[which(names(tbl) == "2")]
  project.summary$aliquot[i] = ncol(data) 

}

# write summary output
# 11,368 aliquots and 19,729 protein coding genes in chr1-22/X

write.table(gene.summary, "gistic.broad.gene.summary.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(aliquot.summary, "gistic.broad.aliquot.summary.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(project.summary, "gistic.broad.project.summary.txt", col.names=T, row.names=F, sep="\t", quote=F)

temp = apply(project.summary[, 2:6], 2, sum)

# 
gene.info.file = "~/SCRATCH/gencode/gencode.v22.genes.txt"
gene.info = fread(gene.info.file)
contigs = paste0("chr", c(1:22, "X"))
pc.genes = gene.info[seqname %in% contigs & gene_type == "protein_coding"]$gene_id









