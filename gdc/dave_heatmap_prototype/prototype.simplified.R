library(data.table)


# define the sample size and number of selected genes we want to test
nsample <- 1000
ngene <- 50


##########################################################
# Step 1. Load data from TSV file
##########################################################
expr <- fread("zcat cptac3.expr.tsv.gz")
# we cam also load data with the following command without data.table packages, but it will be very slow
# expr <- read.table(gzfile("cptac3.expr.tsv.gz"), stringsAsFactors=F)

# let's ignore gene names for now
expr <- as.matrix(expr[, 2:ncol(expr), with=F])

# Random sampling to archieve the right sample size for benchmark
expr <- expr[, sample(ncol(expr), nsample, replace = T)]


##########################################################
# Alternative of Step 1. Load data from RDS file
##########################################################
expr <- readRDS("cptac3.expr.rds")



##########################################################
# Step 2. Data Processing
##########################################################

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
genestats <- genestats[high.var.index][1:ngene]
expr3 <- expr2[high.var.index, ][1:ngene, ]

# z-transformation
expr4 <- (expr3 - genestats$mean) / genestats$sd 
mid.time <- proc.time()


##########################################################
# Step 3. Clustering
##########################################################
pdf(file=NULL)
heatmap(expr4, scale="none")
dev.off()

# benchmark output
end.time <- proc.time()


##########################################################
# Step 4. Benchmark Results
##########################################################

processing.time <- mid.time- start.time
clustering.time <- end.time - mid.time
