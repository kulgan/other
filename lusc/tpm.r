library(data.table)

a = fread("luad.lusc.csv")

# extract genes
f = a$file[1]
gene = fread(paste0("zcat ./data/", f), h=F)$V1

# extract count
out = matrix(0, length(gene), nrow(a))
for(i in 1:nrow(a)){
	f = a$file[i]
	out[, i] = fread(paste0("zcat ./data/", f), h=F)$V2
}
colnames(out) = a$barcode
rownames(out) = gene
count = out

# make tpm
total = apply(count, 2, sum)
tpm = t(t(count)/ total * 1e6)

# low expression filter
w1 = which(apply(tpm, 1, function(x) sum(x) < 1 ) > ncol(tpm)*0.8)
w2 = which(apply(tpm, 1, function(x) mean(x) < 1))
tpm2 = tpm[-c(w1, w2), ]

# mad filter
mad = apply(tpm2, 1, function(x) mad(x, constant=1))
w3 = which(mad > mad[order(-mad)][1501])
tpm3 = tpm2[w3, ]

# log transformation and remove dup
tpm4 = log2(tpm3 + 1)
w4 = which(duplicated(substr(colnames(tpm4), 1, 12)))
tpm5 = tpm4[, -w4]

# euro
# pc2 > 0.05 -15*pc1
euro = fread("euro.pop",h=F)$V1
w5 = which(substr(colnames(tpm5), 1, 12) %in% euro)
tpm6 = tpm5[, w5]

# PCA
pca = prcomp(t(tpm6), center = TRUE, scale. = TRUE) 
x = data.frame(pca$x)
x$aliquot = rownames(x)
x$patient = substr(x$aliquot, 1, 12)
x$disease = "LUSC"
x$disease[which(x$aliquot %in% a$barcode[which(a$project =="LUAD")])] = "LUAD"



library("Rtsne")
disease = a$project[match(colnames(tpm2), a$barcode)]
pca2 = prcomp(t(tpm2), center = TRUE, scale. = TRUE) 
x2 = pca2$x[, 1:200]
set.seed(1)
tsne <- Rtsne(x2, dims = 3, pca = F)
pdata = data.frame(tsne$Y)
colnames(pdata) = c("feature_1", "feature_2", "feature_3")
pdata$disease = disease
ggplot(pdata, aes(x=feature_1, y=feature_2, col=disease)) + geom_point()


library("shiney")
sdata = pdata
sdata$disease = paste0(sdata$disease, ":sample", 1:nrow(sdata))
sdata = sdata[, c("disease", "feature_1", "feature_2", "feature_3")]
colnames(sdata) = c("", "PC1", "PC2", "PC3")
sdata$PC1 = sdata$PC1 / 1e4
sdata$PC2 = sdata$PC2 / 1e4
sdata$PC3 = sdata$PC3 / 1e4
write.table(sdata, "/Users/zhenyu/github/Shiny_fun/LUAD_LUSC.6-15-16/LUSC_LUAD.top3pc.csv", col.names=T, row.names=F, quote=F, sep=",")
