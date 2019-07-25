library(data.table)

# set confidence threshold
confidence.threshold = 0.1
cohorts = fread("cohort", h=F)$V1

# Set filenames
method1 = "birdseed"
call1.file = paste0("grep -v ^# ", method1, ".calls.txt")
conf1.file = paste0("grep -v ^# ", method1, ".confidences.txt")

method2 = "birdseed-v2"
call2.file = paste0("grep -v ^# ", method2, ".calls.txt")
conf2.file = paste0("grep -v ^# ", method2, ".confidences.txt")

# method3 = "brlmm-p"
# call3.file = paste0("grep -v ^# ", method3, ".calls.txt")
# conf3.file = paste0("grep -v ^# ", method3, ".confidences.txt")

for(cohort in cohorts) {
  print(cohort)
  setwd(cohort)

  # read BirdSeed v1 and v2
  call1 = fread(cmd = call1.file)
  conf1 = fread(cmd = conf1.file)
  probeid = call1$probeset_id
  call1 = as.matrix(call1[, 2:ncol(call1)])
  conf1 = as.matrix(conf1[, 2:ncol(conf1)])

  call2 = fread(cmd = call2.file)
  conf2 = fread(cmd = conf2.file)
  call2 = as.matrix(call2[, 2:ncol(call2)])
  conf2 = as.matrix(conf2[, 2:ncol(conf2)])

  # get consensus call 
  w = unique(c(which(call1 != call2), which(call1 == -1), which(conf1 > confidence.threshold), which(conf2 > confidence.threshold)))
  call1[w] = NA

  # collect metrix
  geno = apply(call1, 1, function(x) sum(x, na.rm=T))
  missing = apply(call1, 1, function(x) sum(is.na(x)))
  qc = data.table(probeid, geno, missing)

  call1 = cbind(probeid, call1)

  fwrite(call1, "consensus.call.tsv", col.names=T, row.names=F, sep="\t", quote=F)
  fwrite(qc, "consensus.call.metric.tsv", col.names=T, row.names=F, sep="\t", quote=F)

  setwd("../")
}


