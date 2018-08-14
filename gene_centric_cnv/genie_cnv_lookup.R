library(data.table)
library(dplyr)
library(ggplot2)
library(refGenome)


gm = fread("Homo_sapiens.gene_info.txt")

# read and parse hg19 Gencode v85 gtf
ens = ensemblGenome()
read.gtf(ens, "Homo_sapiens.GRCh37.85.gtf")
gene = data.table(
		getGenePositions(ens) %>%
		select(seqid, start, end, gene_name, gene_biotype, gene_id)
		)
