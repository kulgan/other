library(data.table)
library(dplyr)

setwd("/Users/zhenyu/github/other/gdc/gene/")


gene = fread("gencode.v22.genes.txt")
refseqgene = fread("gene_RefSeqGene") # 6000+ RefSeqGenes
refseq = fread("mart_export.txt")
cosmic = fread("cancer_census_gene_20181017.txt")
bkgd = fread("DNMR-average.txt")

# clean ensembl <-> refseq protein_coding gene conversion table from biomart
names(refseq) = c("ensg", "enst", "refseq") 
refseq = refseq[refseq != ""]

# cleanup mutsig genes
# mutsig = fread("TCGA.MutSig.summary.txt") # 3516
# mutsig$gene[which(mutsig$gene == "C19orf10")] = "MYDGF"
# mutsig$gene[which(mutsig$gene == "WDR16")] = "CFAP52"
# mutsig$gene[which(mutsig$gene == "C4orf6")] = "LINC01587"
# mutsig$gene[which(mutsig$gene == "CBWD6")] = "CBWD7"
# mutsig$gene[which(mutsig$gene == "FAM123B")] = "AMER1"

# mutsig2 = mutsig[min.q < 0.01] # 1865
# mutsig3 = mutsig[min.q < 0.001] # 1141
# mutsig10 = mutsig[min.q < 1e-10] # 555

# more annotation on gene model
gene = gene %>% 
  mutate(ensg = gsub("^(\\w+)\\.\\d+$", "\\1", gene_id), 
         cosmic = ensg %in% cosmic$ensg, # excluding IG and TR genes, 
         pc_gene = gene_type == "protein_coding", 
         pc_curated = gene_type == "protein_coding" & gene_status == "KNOWN" & source == "HAVANA", 
         refseqgene = gene_name %in% refseqgene$Symbol, # only 6719 out of 6817 has name matches
         refseq_pc = ensg %in% refseq$ensg, # only 19023 out of 20930 ensg name matches
         refseq_curated_protein_coding = refseq_pc == T & pc_gene == T,  # 19003 out of 19814 protein_coding genes 
         mutsig = gene_name %in% mutsig[["gene"]], 
         mutsig2 = gene_name %in% mutsig2[["gene"]], 
         mutsig3 = gene_name %in% mutsig3[["gene"]], 
         mutsig10 = gene_name %in% mutsig10[["gene"]]) %>% 
  setDT()
m = match(gene$gene_name, bkgd$Gene)
gene$background_rate = bkgd$Averge[m]

write.table(gene, "gencode.v22.gene.with_some_partitial_refseq_cosmic_annotations.txt", col.names=T, row.names=F, sep="\t", quote=F)



mutsig = fread("mutsig.top20.per.project.txt", h=F)$V1 # 402
mutsig[which(mutsig == "C19orf10")] = "MYDGF"
mutsig[which(mutsig == "WDR16")] = "CFAP52"
mutsig[which(mutsig == "C4orf6")] = "LINC01587"
mutsig[which(mutsig == "CBWD6")] = "CBWD7"
mutsig[which(mutsig == "FAM123B")] = "AMER1"
mutsig.genes = mutsig

# more annotation on gene model
gene = gene %>% 
  mutate(ensg = gsub("^(\\w+)\\.\\d+$", "\\1", gene_id), 
         cosmic = ensg %in% cosmic$ensg, # excluding IG and TR genes, 
         pc_gene = gene_type == "protein_coding", 
         pc_curated = gene_type == "protein_coding" & gene_status == "KNOWN" & source == "HAVANA", 
         refseqgene = gene_name %in% refseqgene$Symbol, # only 6719 out of 6817 has name matches
         refseq_pc = ensg %in% refseq$ensg, # only 19023 out of 20930 ensg name matches
         refseq_curated_protein_coding = refseq_pc == T & pc_gene == T,  # 19003 out of 19814 protein_coding genes 
         mutsig = gene_name %in% mutsig.genes) %>% 
  setDT()
m = match(gene$gene_name, bkgd$Gene)
gene$background_rate = bkgd$Averge[m]

tcga = fread("gunzip -c tcga.gene_mutation_count.txt.gz")
m = match(tcga$gene_id, gene$gene_id)
tcga$cosmic = gene$cosmic[m]
tcga$pc_curated = gene$pc_curated[m]
tcga$refseqgene = gene$refseqgene[m]
tcga$background_rate = gene$background_rate[m]
tcga$refseq = gene$refseq_curated_protein_coding[m]
tcga$mutsig = gene$mutsig[m]




write.table(gene, "gencode.v22.gene.with_some_partitial_refseq_cosmic_annotations.txt", col.names=T, row.names=F, sep="\t", quote=F)



# annotate count
tcga = fread("gunzip -c tcga.gene_mutation_count.txt.gz")
m = match(tcga$gene_id, gene$gene_id)
tcga$cosmic = gene$cosmic[m]
tcga$pc_curated = gene$pc_curated[m]
tcga$refseqgene = gene$refseqgene[m]
tcga$background_rate = gene$background_rate[m]
tcga$refseq = gene$refseq_curated_protein_coding[m]
tcga$mutsig = gene$mutsig[m]
tcga$mutsig2 = gene$mutsig2[m]
tcga$mutsig3 = gene$mutsig3[m]
tcga$mutsig10 = gene$mutsig10[m]


min_length = seq(0, 50000, 100)
n = length(min_length)
out = data.table(min_length = min_length, 
                 cosmic20 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                            sum(tcga[order][refseq == T][1:20]$cosmic)}), 
                 cosmic100 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                             sum(tcga[order][refseq == T][1:100]$cosmic)}), 
                 cosmic200 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                             sum(tcga[order][refseq == T][1:200]$cosmic)}), 
                 cosmic1000 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                              sum(tcga[order][refseq == T][1:1000]$cosmic)}), 
                 mutsig_20 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                            sum(tcga[order][refseq == T][1:20]$mutsig)}), 
                 mutsig_100 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                             sum(tcga[order][refseq == T][1:100]$mutsig)}), 
                 mutsig_200 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                             sum(tcga[order][refseq == T][1:200]$mutsig)}), 
                 mutsig_1000 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                              sum(tcga[order][refseq == T][1:1000]$mutsig)}), 
                 mutsig2_20 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                            sum(tcga[order][refseq == T][1:20]$mutsig2)}), 
                 mutsig2_100 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                             sum(tcga[order][refseq == T][1:100]$mutsig2)}), 
                 mutsig2_200 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                             sum(tcga[order][refseq == T][1:200]$mutsig2)}), 
                 mutsig2_1000 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                              sum(tcga[order][refseq == T][1:1000]$mutsig2)}), 
                 mutsig3_20 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                            sum(tcga[order][refseq == T][1:20]$mutsig3)}), 
                 mutsig3_100 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                             sum(tcga[order][refseq == T][1:100]$mutsig3)}), 
                 mutsig3_200 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                             sum(tcga[order][refseq == T][1:200]$mutsig3)}), 
                 mutsig3_1000 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                              sum(tcga[order][refseq == T][1:1000]$mutsig3)}), 
                 mutsig10_20 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                            sum(tcga[order][refseq == T][1:20]$mutsig10)}), 
                 mutsig10_100 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                             sum(tcga[order][refseq == T][1:100]$mutsig10)}), 
                 mutsig10_200 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                             sum(tcga[order][refseq == T][1:200]$mutsig10)}), 
                 mutsig10_1000 = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                              sum(tcga[order][refseq == T][1:1000]$mutsig10)}))

d = data.table(minimum_length = rep(min_length, 4), 
               num_cosmic = c(out$cosmic20, out$cosmic100, out$cosmic200, out$cosmic1000), 
               num_mutsig = c(out$mutsig_20, out$mutsig_100, out$mutsig_200, out$mutsig_1000),
               num_mutsig2 = c(out$mutsig2_20, out$mutsig2_100, out$mutsig2_200, out$mutsig2_1000),
               num_mutsig3 = c(out$mutsig3_20, out$mutsig3_100, out$mutsig3_200, out$mutsig3_1000),
               num_mutsig10 = c(out$mutsig10_20, out$mutsig10_100, out$mutsig10_200, out$mutsig10_1000),
               gdc_top_genes = factor(c(rep("top20", n), rep("top100", n), rep("top200", n), rep("top1000", n)), levels = c("top1000", "top200", "top100", "top20")))
write.table(d, "data.to.plot.txt", col.names=T, row.names=F, sep="\t", quote=F)

library(ggplot2)

ggplot(d, aes(minimum_length, num_cosmic, col = gdc_top_genes)) + 
  geom_line() + 
  ggtitle("Number of COSMIC Genes in GDC Top Genes") + 
  xlab("Min Length (bps)") + 
  ylab("Number of Cancer Genes (COSMIC)") +
  theme(legend.position=c(0.8, 0.8), text = element_text(size=20)) 

library(gridExtra)

p1 = ggplot(d, aes(minimum_length, num_mutsig, col = gdc_top_genes)) + 
  geom_line() + 
  ggtitle("Number of MutSig Genes (FDR 0.1)") + 
  xlab("Min Length (bps)") + 
  ylab("") +
  theme(legend.position="none", text = element_text(size=20)) 

p2 = ggplot(d, aes(minimum_length, num_mutsig2, col = gdc_top_genes)) + 
  geom_line() + 
  ggtitle("Number of MutSig Genes (FDR 0.01)") + 
  xlab("Min Length (bps)") + 
  ylab("") +
  theme(legend.position="none", text = element_text(size=20)) 

p3 = ggplot(d, aes(minimum_length, num_mutsig3, col = gdc_top_genes)) + 
  geom_line() + 
  ggtitle("Number of MutSig Genes (FDR 0.001)") + 
  xlab("Min Length (bps)") + 
  ylab("") +
  theme(legend.position="none", text = element_text(size=20)) 

p10 = ggplot(d, aes(minimum_length, num_mutsig10, col = gdc_top_genes)) + 
  geom_line() + 
  ggtitle("Number of MutSig Genes (FDR 1e-10)") + 
  xlab("Min Length (bps)") + 
  ylab("") +
  theme(legend.position="none", text = element_text(size=20))     

grid.arrange(p1, p2, p3, p10, nrow=2)


len = c(0, 500, 1000, 5000, 10000, 15000, 20000, 30000)
top20 = sapply(len, function(x) {order = with(tcga, order(-count/pmax(exon_length, x))); temp = tcga[order][refseq == T][1:20]; c(temp$gene_name, sum(temp$cosmic))})
rownames(top20) = c(1:20, "cosmic_score")
colnames(top20) = len

# test ranking of cosmic genes in GDC 
min_length = 10000
ranking = tcga[order(-count/pmax(exon_length, min_length))][refseq == T]
ranking$ensg = sapply(ranking$gene_id, function(x) strsplit(x, "\\.")[[1]][1])


cosmic.all = fread("gunzip -c CosmicMutantExport.tsv.gz")
setnames(cosmic.all, gsub("-", "_", gsub(" ", "_", tolower(names(cosmic.all)))))
cosmic.gw = cosmic.all %>% 
  filter(genome_wide_screen == "y", 
         mutation_somatic_status %in% c("Confirmed somatic variant", "Reported in another cancer sample as somatic"), 
         ! mutation_description %in% c("Substitution - coding silent", "Unknown")) %>%
  select(gene_name, accession_number) %>%
  mutate(gene_name = gsub("^(.+)_(ENST|NM_).+$", "\\1", gene_name)) %>% 
  group_by(gene_name) %>% 
  summarise(count = n()) %>%
  arrange(-count) %>%
  setDT()





names(cosmic) = gsub(" ", "_", tolower(names(cosmic)))
cosmic2 = cosmic %>% 
  filter(tier == "1", 
         grepl("oncogene", role_in_cancer), 
         gene_type == "protein_coding", 
         somatic == "yes", 
         grepl("(F|Mis|N|S)", mutation_types)) %>%
  setDT()


ranking$cos = ranking$ensg %in% cosmic2$ensg
ranking$hallmark = ranking$ensg %in% cosmic2[hallmark=="Yes"]$ensg
ranking$index = 1:nrow(ranking)
p1 = ggplot(ranking[cos==T], aes(index)) + 
  geom_histogram(binwidth = 200) + 
  xlab("Gene Ranking") + 
  ggtitle("Histogram of 125 COSMIC gene rankings") +
  theme(text = element_text(size=20)) 

p2 = ggplot(ranking[mutsig==T], aes(index)) + 
  geom_histogram(binwidth = 200) + 
  xlab("Gene Ranking") + 
  ggtitle("Histogram of 402 MutSig gene rankings") +
  theme(text = element_text(size=20)) 

ggplot(ranking[mutsig==T], aes(index)) + 
  geom_histogram(binwidth = 200) + 
  xlab("Gene Ranking") + 
  ggtitle("Histogram of 125 COSMIC gene rankings") +
  theme(text = element_text(size=20)) 


ggplot(ranking[cos==T], aes(index, fill="blue")) + 
  geom_histogram(binwidth = 200) + 
  xlab("Gene Ranking") + 
  ggtitle("Histogram of 125 COSMIC gene rankings") +
  theme(text = element_text(size=20)) 

ggplot(ranking[cosmic==T], aes(index, fill="#33AADE")) + 
  geom_histogram(binwidth = 200) + 
  xlab("Gene Ranking") + 
  ggtitle("Histogram of 125 COSMIC gene rankings") +
  theme(text = element_text(size=20)) 


ranking2 = ranking[order(-count)]
ranking2$index = 1:nrow(ranking2)
ggplot(ranking2[cos==T], aes(index)) + 
  geom_histogram(binwidth = 500) + 
  xlab("Gene Ranking") + 
  ggtitle("Histogram of 125 COSMIC gene rankings") +
  theme(text = element_text(size=20)) 



# count only
out1 = data.table(order = 1:20)
out1$count = tcga[order(-count/background_rate)]$gene_name[1:20]
out1$cosmic = tcga[order(-count/background_rate)][cosmic == T]$gene_name[1:20]
out1$curated = tcga[order(-count/background_rate)][pc_curated == T]$gene_name[1:20]
out1$refseqgene = tcga[order(-count/background_rate)][refseqgene == T]$gene_name[1:20]
out1$refseqgene = tcga[order(-count/background_rate)][refseq == T]$gene_name[1:20]

# normalization by length
out2 = data.table(order = 1:20)
out2$count = tcga[order(-count/exon_length)]$gene_name[1:20]
out2$cosmic = tcga[order(-count/exon_length)][cosmic == T]$gene_name[1:20]
out2$curated = tcga[order(-count/exon_length)][pc_curated == T]$gene_name[1:20]
out2$refseqgene = tcga[order(-count/exon_length)][refseqgene == T]$gene_name[1:20]
out2$refseq = tcga[order(-count/exon_length)][refseq == T]$gene_name[1:20]

# normalization by length (with min_length)
out3 = data.table(order = 1:20)
out3$count = tcga[order(-count/pmax(exon_length, 2000))]$gene_name[1:20]
out3$cosmic = tcga[order(-count/pmax(exon_length, 2000))][cosmic == T]$gene_name[1:20]
out3$curated = tcga[order(-count/pmax(exon_length, 2000))][pc_curated == T]$gene_name[1:20]
out3$refseqgene = tcga[order(-count/pmax(exon_length, 2000))][refseqgene == T]$gene_name[1:20]

out4 = data.frame(len = rep(seq(1000, 40000, 1000), 4), cosmic = NA, method = c(rep("min_curated", 40), rep("min_refseqgene", 40), rep("padding_curated", 40), rep("padding_refseqgene", 40)))
for(i in 1:40) {
  min_length = 1000 * i
  out4 = out4 %>% 
    mutate(cosmic = ifelse(len == min_length & method == "min_curated", sum(tcga[order(-count/pmax(exon_length, min_length))][pc_curated == T][1:20]$cosmic), cosmic), 
           cosmic = ifelse(len == min_length & method == "min_refseqgene", sum(tcga[order(-count/pmax(exon_length, min_length))][refseqgene == T][1:20]$cosmic), cosmic),
           cosmic = ifelse(len == min_length & method == "padding_curated", sum(tcga[order(-count/(exon_length + min_length))][pc_curated == T][1:20]$cosmic), cosmic),
           cosmic = ifelse(len == min_length & method == "padding_refseqgene", sum(tcga[order(-count/(exon_length + min_length))][refseqgene == T][1:20]$cosmic), cosmic)) %>%
    setDT()
}


out4 = data.frame(len = rep(seq(1000, 40000, 1000), 4), cosmic = NA, method = c(rep("min_curated", 40), rep("min_refseqgene", 40), rep("padding_curated", 40), rep("padding_refseqgene", 40)))
for(i in 1:40) {
  min_length = 1000 * i
  out4 = out4 %>% 
    mutate(cosmic = ifelse(len == min_length & method == "min_curated", sum(tcga[order(-count/pmax(exon_length, min_length))][pc_curated == T][1:20]$cosmic), cosmic), 
           cosmic = ifelse(len == min_length & method == "min_refseqgene", sum(tcga[order(-count/pmax(exon_length, min_length))][refseqgene == T][1:20]$cosmic), cosmic),
           cosmic = ifelse(len == min_length & method == "padding_curated", sum(tcga[order(-count/(exon_length + min_length))][pc_curated == T][1:20]$cosmic), cosmic),
           cosmic = ifelse(len == min_length & method == "padding_refseqgene", sum(tcga[order(-count/(exon_length + min_length))][refseqgene == T][1:20]$cosmic), cosmic))
}



min_length = seq(100, 50000, 100)
out6 = data.table(min_length = min_length, 
                  cosmic = sapply(min_length, function(x) {order = order(with(tcga, -count/pmax(exon_length, x))); 
                                                                 sum(tcga[order][refseq == T][1:20]$cosmic)}))

len = c(0, 500, 1000, 10000, 12000, 13000, 20000, 30000)
out7 = sapply(len, function(x) {order = with(tcga, order(-count/pmax(exon_length, x))); temp = tcga[order][refseq == T][1:20]; c(temp$gene_name, sum(temp$cosmic))})
rownames(out7) = c(1:20, "cosmic_score")
colnames(out7) = len


min_length = 10000
out5 = data.frame(min_curated = tcga[order(-count/pmax(exon_length, min_length))][pc_curated == T]$gene_name[1:20],
                  min_refseqgene = tcga[order(-count/pmax(exon_length, min_length))][refseqgene == T]$gene_name[1:20],
                  padding_curated = tcga[order(-count/(exon_length + min_length))][pc_curated == T]$gene_name[1:20],
                  padding_refseqgene = tcga[order(-count/(exon_length + min_length))][refseqgene == T]$gene_name[1:20])

min_length = 11000
out6 = data.frame(min_curated = tcga[order(-count/pmax(exon_length, min_length))][pc_curated == T]$gene_name[1:20],
                  min_refseqgene = tcga[order(-count/pmax(exon_length, min_length))][refseqgene == T]$gene_name[1:20],
                  padding_curated = tcga[order(-count/(exon_length + min_length))][pc_curated == T]$gene_name[1:20],
                  padding_refseqgene = tcga[order(-count/(exon_length + min_length))][refseqgene == T]$gene_name[1:20])


min_length = 2000
out7 = data.frame(min_curated = tcga[order(-count/pmax(exon_length, min_length))][pc_curated == T]$gene_name[1:20],
                  min_refseqgene = tcga[order(-count/pmax(exon_length, min_length))][refseqgene == T]$gene_name[1:20],
                  padding_curated = tcga[order(-count/(exon_length + min_length))][pc_curated == T]$gene_name[1:20],
                  padding_refseqgene = tcga[order(-count/(exon_length + min_length))][refseqgene == T]$gene_name[1:20])


# input attribute field from GTF, and an attribute feature name, return attribute values
extractAttr = function(gtf_attr, feature) {
  pattern = paste0(".*", feature, " \"((\\w|\\.)+)\".*")
  value = gsub(pattern, "\\1", gtf_attr)

  # when no pattern exist, set value to NA
  value[which(value == gtf_attr)] = NA
  return(value)
}

# read gtf
gtf = fread("gunzip -c gencode.v22.annotation.gtf.gz")

# label gtf and extract properties
# set ensg as ensembl gene_id without version number
# set pc_curated = T for KNOWN HAVANA protein_coding genes
colnames(gtf) = c("seqname", "source", "feature", "start", "end",  "score", "strand", "frame", "attribute")
gene = gtf[which(gtf$feature=="gene"), ] %>% 
       mutate(gene_id = extractAttr(attribute, "gene_id"), 
              gene_type = extractAttr(attribute, "gene_type"), 
              gene_status = extractAttr(attribute, "gene_status"), 
              gene_name = extractAttr(attribute, "gene_name"), 
              havana_gene = extractAttr(attribute, "havana_gene"), 
              level = gsub(".+level (\\d).+", "\\1", attribute), 
              ensg = gsub("^(\\w+)\\.\\d+$", "\\1", gene_id), 
              pc_curated = gene_type == "protein_coding" & gene_status == "KNOWN" & source == "HAVANA") %>% 
       setDT()       
}