library(data.table)
library(dplyr)


setwd("~/SCRATCH/gistic/data/")
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
diseases = fread("./d", h=F)$V1

gene.gain.summary = array(0, length(pc.genes))
names(gene.gain.summary) = pc.genes
gene.loss.summary = gene.gain.summary

aliquot.gain.summary = numeric(0)
aliquot.loss.summary = numeric(0)

project.summary = data.table(project = diseases, t1=0, t2=0, t3=0, total=0)


for (i in 1:length(diseases)) {
  disease = diseases[i]
  print(disease)

  # read file
  all.file = paste0("./", disease, "/all_data_by_genes.txt")
  all = fread(all.file)

  # filter for protein coding genes
  all = all[all$"Gene Symbol" %in% pc.genes]

  # separate data into metadata and copy number values, and generate thresholed values
  data = unlist(all[, 4:ncol(all)])
  project.summary$total[i] = length(data)
  project.summary$t1[i] = sum(data < 0.1 & data > -0.1)
  project.summary$t2[i] = sum(data < 0.2 & data > -0.2)
  project.summary$t3[i] = sum(data < 0.3 & data > -0.3)
 
}






threshold = function(data, t1, t2, t3, t4) {
  out = array(0, length(data))
  out[which(data < t2)] = -1
  out[which(data < t1)] = -2
  out[which(data > t3)] = 1
  out[which(data > t4)] = 2
  return(out)
}

for (i in 1:length(diseases)) {
  disease = diseases[i]
  print(disease)

  # sample cutoffs
  sample.cutoff.file = paste0("./", disease, "/sample_cutoffs.txt")
  sample.cutoff = fread(sample.cutoff.file)
  sample.cutoff = sample.cutoff %>% 
    mutate(t1 = Low - noise, 
           t2 = -noise,
           t3 = noise,
           t4 = High + noise) %>%
    setDT()

  temp = data$"154e9edd-b354-4e62-be0f-8b26ebddffbc"
  out = threshold(temp, sample.cutoff[Sample == "154e9edd-b354-4e62-be0f-8b26ebddffbc"]$t1, sample.cutoff[Sample == "154e9edd-b354-4e62-be0f-8b26ebddffbc"]$t2, sample.cutoff[Sample == "154e9edd-b354-4e62-be0f-8b26ebddffbc"]$t3, sample.cutoff[Sample == "154e9edd-b354-4e62-be0f-8b26ebddffbc"]$t4)
  temp1 = th$"154e9edd-b354-4e62-be0f-8b26ebddffbc"

  # read file
  broad.file = paste0("./", disease, "/all_data_by_genes.txt")
  broad = fread(broad.file)



  data = broad[, 4:ncol(broad)]

  th = fread("ACC/all_thresholded.by_genes.txt")
  th = th[, 4:ncol(th)]

  # filter for protein coding genes
  broad = broad[broad$"Gene Symbol" %in% pc.genes]

  # separate data into metadata and copy number values, and generate thresholed values
  meta = focal[, 1:3]
  data = focal[, 4:ncol(focal)]
  data = replace(data, data > noise, 1)
  data = replace(data, data < -noise, -1)
  data = replace(data, data >= -noise & data <= noise, 0)  

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
ncase = length(aliquot.gain.summary)
ngene = length(pc.genes)

gene.summary = data.table(gene = names(gene.gain.summary), gain = gene.gain.summary, loss = gene.loss.summary) %>%
  mutate(neutral = ncase - loss - gain)
write.table(gene.summary, "gene.event.summary.txt", col.names=T, row.names=F, sep="\t", quote=F)

aliquot.summary = data.table(aliquot = names(aliquot.gain.summary), gain = aliquot.gain.summary, loss = aliquot.loss.summary) %>%
  mutate(neutral = ngene - loss - gain)
write.table(aliquot.summary, "aliquot.event.summary.txt", col.names=T, row.names=F, sep="\t", quote=F)

write.table(project.summary, "project.summary.txt", col.names=T, row.names=F, sep="\t", quote=F)












