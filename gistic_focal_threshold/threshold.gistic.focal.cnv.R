library(data.table)


setwd("~/SCRATCH/gistic/data/pc_thresholded_03.by_genes")
noise = 0.3

##########################################################
# generate protein-coding gene list
##########################################################

gene.info.file = "~/SCRATCH/gencode/gencode.v22.genes.txt"
gene.info = fread(gene.info.file)
contigs = paste0("chr", c(1:22, "X"))

pc.genes = gene.info[seqname %in% contigs & gene_type == "protein_coding"]$gene_id


####################################################################################################################
# convert focal_data_by_genes to scores, and only keep protein-coding genes in chr1-22/X
####################################################################################################################
diseases = fread("../d", h=F)$V1

gene.gain.summary = array(0, length(pc.genes))
names(gene.gain.summary) = pc.genes
gene.loss.summary = gene.gain.summary

aliquot.gain.summary = numeric(0)
aliquot.loss.summary = numeric(0)

project.summary = data.table(project = diseases, loss = 0, neutral = 0, gain = 0, aliquot = 0)


for (i in 1:length(diseases)) {
  disease = diseases[i]
  print(disease)
  
  # read file
  focal.file = paste0("../", disease, "/focal_data_by_genes.txt")
  focal = fread(focal.file)

  # filter for protein coding genes
  focal = focal[focal$"Gene Symbol" %in% pc.genes]

  # separate data into metadata and copy number values, and generate thresholed values
  meta = focal[, 1:3]
  data = focal[, 4:ncol(focal)]
  data = replace(data, data > noise, 1)
  data = replace(data, data < -noise, -1)

  # collect metrics
  gene.gain.summary = gene.gain.summary + apply(data, 1, function(x) sum(x == 1))
  gene.loss.summary = gene.loss.summary + apply(data, 1, function(x) sum(x == -1))
  aliquot.gain.summary = c(aliquot.gain.summary, apply(data, 2, function(x) sum(x == 1)))
  aliquot.loss.summary = c(aliquot.loss.summary, apply(data, 2, function(x) sum(x == -1)))
  tbl = table(unlist(data))
  project.summary$loss[i] = tbl[which(names(tbl) == "-1")]
  project.summary$neutral[i] = tbl[which(names(tbl) == "0")]
  project.summary$gain[i] = tbl[which(names(tbl) == "1")]
  project.summary$aliquot[i] = ncol(data) 

  # merge metadata and scores  
  data = cbind(meta, data)
  outfile = paste(disease, "focal_score_by_genes.txt", sep=".")
  write.table(data, outfile, col.names=T, row.names=F, sep="\t", quote=F)
}

# write summary output
# 11,368 aliquots and 19,729 protein coding genes in chr1-22/X

gene.gain.summary = data.table(gene = names(gene.gain.summary), count = gene.gain.summary)
write.table(gene.gain.summary, "gene.gain.summary.txt", col.names=T, row.names=F, sep="\t", quote=F)

gene.loss.summary = data.table(gene = names(gene.loss.summary), count = gene.loss.summary)
write.table(gene.loss.summary, "gene.loss.summary.txt", col.names=T, row.names=F, sep="\t", quote=F)

gene.event.summary = gene.gain.summary
gene.event.summary$count = gene.gain.summary$count + gene.loss.summary$count
write.table(gene.event.summary, "gene.event.summary.txt", col.names=T, row.names=F, sep="\t", quote=F)

aliquot.gain.summary = data.table(aliquot = names(aliquot.gain.summary), count = aliquot.gain.summary)
write.table(aliquot.gain.summary, "aliquot.gain.summary.txt", col.names=T, row.names=F, sep="\t", quote=F)

aliquot.loss.summary = data.table(aliquot = names(aliquot.loss.summary), count = aliquot.loss.summary)
write.table(aliquot.loss.summary, "aliquot.loss.summary.txt", col.names=T, row.names=F, sep="\t", quote=F)

aliquot.event.summary = aliquot.gain.summary
aliquot.event.summary$count = aliquot.gain.summary$count + aliquot.loss.summary$count
write.table(aliquot.event.summary, "aliquot.event.summary.txt", col.names=T, row.names=F, sep="\t", quote=F)

write.table(project.summary, "project.summary.txt", col.names=T, row.names=F, sep="\t", quote=F)












