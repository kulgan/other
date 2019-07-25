setwd("/Users/zhenyu/Downloads/working/Steph/ov3")
# setwd("~/SCRATCH/ov3")
library(data.table)
library(survival)
library(ggplot2)
library(GGally)
library(RANN)

surv = fread("surv.txt")
drug = fread("drug.use.table")
# tsne = fread("combined.tsne3.PC200.dim3.3_3_1_1.perm1000.txt")
load("TCGA.rna.filter1.allgene.pca.rda")
# load("~/SCRATCH/tsne/rna/TCGA.rna.filter1.allgene.pca.rda")
expr = data.frame(pca$x[, 1:3])
expr$patient = substr(rownames(expr), 1, 12)

d = merge(merge(surv, drug, by="patient"), expr, by="patient")
# remove NA
d = d[which(!is.na(apply(d[, 7:10, with=F], 1, sum))), ]

d$tax_and_carbo = (d$Taxol | d$Taxotere) & d$Carboplatin

# Test on Cisplatin on Tax/Carboplatin users
d2 = d[(d$Taxol | d$Taxotere) & d$Carboplatin==TRUE]

# k nearest neighbour
k = 10
nn2 = nn2(d2[, c("PC1", "PC2", "PC3"), with=F], k=k)

d2$k.cis.month = 0
d2$k.nocis.month = 0
for(i in 1:nrow(d2)) {
	idx = nn2$nn.idx[i, 2:k]
	cis = d2$Cisplatin[idx]
	d2$k.cis.month[i] = mean(d2$months[idx[cis]])
	d2$k.nocis.month[i] = mean(d2$months[idx[!cis]])
}

d2$cis = "ambiguous"
d2$cis[d2$k.cis.month - d2$k.nocis.month >= 1] = "better"
d2$cis[d2$k.cis.month - d2$k.nocis.month <= -1] = "worse"

judge = d2$cis
match = with(d2, (cis=="better" & Cisplatin == T) | (cis=="worse" & Cisplatin == F))
amb.index = which(judge=="ambiguous")
match[amb.index ] = NA

diff = survdiff(Surv(d2$months[-amb.index], d2$death[-amb.index]) ~ match[-amb.index])
p0 = pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)

nperm = 1000000
perm = array(0, nperm)
for (i in 1:nperm){
	print(i)
	set.seed(i)
	new.match = sample(match, length(match))
	amb.index = which(is.na(new.match))
	diff = survdiff(Surv(d2$months[-amb.index], d2$death[-amb.index]) ~ match[-amb.index])
	perm[i] = pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
}

sum(p0>=perm)/nperm
[1] 0.274547
not significant

# Test on Carboplatin on Tax/Cisplatin users
d3 = d[(d$Taxol | d$Taxotere) & d$Cisplatin==TRUE]

# k nearest neighbour
k = 10
nn3 = nn2(d3[, c("PC1", "PC2", "PC3"), with=F], k=k)

d3$k.carb.month = 0
d3$k.nocarb.month = 0
for(i in 1:nrow(d3)) {
	idx = nn3$nn.idx[i, 2:k]
	carb = d3$Carboplatin[idx]
	d3$k.carb.month[i] = mean(d3$months[idx[carb]])
	d3$k.nocarb.month[i] = mean(d3$months[idx[!carb]])
}

d3$carb = "ambiguous"
d3$carb[d3$k.carb.month - d3$k.nocarb.month >= 1] = "better"
d3$carb[d3$k.carb.month - d3$k.nocarb.month <= -1] = "worse"

judge = d3$carb
match = with(d3, (carb=="better" & Carboplatin == T) | (carb=="worse" & Carboplatin == F))
amb.index = which(judge=="ambiguous")
match[amb.index ] = NA

diff = survdiff(Surv(d3$months[-amb.index], d3$death[-amb.index]) ~ match[-amb.index])
p0 = pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)

nperm = 10000
perm = array(0, nperm)
for (i in 1:nperm){
	print(i)
	set.seed(i)
	new.match = sample(match, length(match))
	amb.index = which(is.na(new.match))
	diff = survdiff(Surv(d3$months[-amb.index], d3$death[-amb.index]) ~ match[-amb.index])
	perm[i] = pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
}


sum(p0>=perm)/nperm
[1] 0.274547
not significant







ggsurv(survfit(with(data, Surv(months, death) ~ match)))

apply(data[, c("Taxol", "Taxotere", "Carboplatin", "Cisplatin"), with=F], 1, function(x) paste(x, collapse=" "))


final boolean germlineAtRisk = !dbsnpVC.isEmpty() && cosmicVC.isEmpty();
normalLodFilterThreshold = germlineAtRisk ? MTAC.NORMAL_DBSNP_LOD_THRESHOLD : MTAC.NORMAL_LOD_THRESHOLD;
 

for(i in 1:33){
file = f[i]
clin = fread(file, h=F)
w = which(grepl("drug_name", clin$V1))
temp = table(unlist(clin[w, 2:ncol(clin), with=F]))
cat(file, sum(temp), "\n")
print(tail(sort(temp)))
}

for(j in 1:33){
file = f[j]
clin = fread(file, h=F)
w = which(grepl("drug_name|patient.bcr_patient_barcode", clin$V1))
temp = clin[w, ]
out = NULL
for(i in 2:ncol(temp)) {
	out = c(out, paste(sort(unlist(temp[2:nrow(temp), i, with=F])), collapse=" "))
}
cat(file, ncol(temp), "\n")
print(tail(sort(table(out))))
}



file="BRCA.clin.merged.txt"
clin = fread(file, h=F)
w = which(grepl("drug_name|patient.bcr_patient_barcode", clin$V1))
temp = clin[w, ]
