library(data.table)

expr <- fread("zcat CPTAC-3.rnaseq.expression.htseq_counts.tsv.gz")

# TPM normalization
genes <- expr$ensembl_gene_id
expr <- as.matrix(expr[, -1])
colsum <- apply(expr, 2, sum) # need to transform to numeric first if count too large
expr <- t(t(expr)/ colsum) * 1e6

# filter for protein coding
gencode <- fread("../segs/gencode.v22.genes.txt")
expr <- expr[which(genes %in% gencode[gene_type=="protein_coding"]$gene_id), ]
expr <- log2(expr + 1)


bm.se <- function(expr, nsample, top) {
  # sampling to archieve the right size for benchmark
  sample.index <- sample(ncol(expr), nsample, replace = T)
  expr <- expr[, sample.index]

  # start benchmark
  start.time <- proc.time()

  # collect average expression and filter out low expression genes
  genestats <- data.table(mean = rowMeans(expr))
  high.expr.index <- which(rowMeans(expr) >= 1) 
  expr2 <- expr[high.expr.index, ]
  genestats <- genestats[high.expr.index]

  # collect sd expression and filter out and only keep genes with high expression variations
  genestats$sd <- apply(expr2, 1, sd)
  genestats$se <- genestats$sd / genestats$mean
  high.var.index <- order(-genestats$se)
  genestats <- genestats[high.var.index][1:top]
  expr3 <- expr2[high.var.index, ][1:top, ]

  # z-transformation
  expr4 <- (expr3 - genestats$mean) / genestats$sd 
  mid.time <- proc.time()

  # headmap
  png("temp.png")
  heatmap(expr4, scale="none")
  dev.off()

  # benchmark output
  end.time <- proc.time()
  ellapse.time <- c((mid.time- start.time)[3], (end.time - mid.time)[3])
  return(ellapse.time)
}

d <- expand.grid(ngene = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 15000), 
                 nsample = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000))

d$se.data.elapsed <- NA
d$se.plot.elapsed <- NA

for(i in 1:nrow(d)) {
  print(i)
  if(d$nsample[i] * d$ngene[i] <= 1e7) {
    d[i, 3:4] <- bm.se(expr, d$nsample[i], d$ngene[i])
  }
}


bm.mad <- function(expr, nsample, top) {
  # sampling to archieve the right size for benchmark
  sample.index <- sample(ncol(expr), nsample, replace = T)
  expr <- expr[, sample.index]

  # start benchmark
  start.time <- proc.time()

  # collect average expression and filter out low expression genes
  expr2 <- expr[rowMeans(expr) >= 1, ]

  # collect mad and select top genes
  m <- apply(expr2, 1, mad)
  expr3 <- expr2[order(-m)[1:top], ]


  # z-transformation
  m <- rowMeans(expr3)
  s <- apply(expr3, 1, sd)
  expr4 <- (expr3 - m) / s
  mid.time <- proc.time()

  # headmap
  png("temp.png")
  heatmap(expr4, scale="none")
  dev.off()

  # benchmark output
  end.time <- proc.time()
  ellapse.time <- c((mid.time- start.time)[3], (end.time - mid.time)[3])
  return(ellapse.time)
}


d$mad.data.elapsed <- NA
d$mad.plot.elapsed <- NA

for(i in c(1:39, 41:48)) {
  print(i)
  d[i, 5:6] <- bm.mad(expr, d$nsample[i], d$ngene[i])
}


d <- expand.grid(ngene = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 15000), 
                 nsample = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000))

expr2 <- expr[, sample(ncol(expr), 10000, replace = T)]

bm.heatmap <- function(expr, nsample, ngene) {
  if(nsample * ngene > 1e7) {
    return(NA)
  } else {
    expr <- expr[1:ngene, 1:nsample]
    mid.time <- proc.time()
    # headmap
    pdf(file=NULL)
    heatmap(expr, scale="none")
    dev.off()

    # benchmark output
    end.time <- proc.time()
    ellapse.time <- (end.time - mid.time)[3]
    print(ellapse.time)
    return(ellapse.time)
  }
}

for(i in 1:nrow(d)) {
  print(i)
  d$time[i] <- bm.heatmap(expr2, d$nsample[i], d$ngene[i])
  print(d$time[i])
  print("")
}



