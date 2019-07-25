library(data.table)
library(CNTools)

# read data
seg <- fread("gunzip -c account_axiom.segs.gz")



names(seg) <- c("ID", "chrom", "arm", "loc.start", "loc.end", "num.mark", "seg.mean")
cnseg <- CNSeg(seg)

cnseg <- CNSeg(seg[, -"arm", with=F])


rdseg <- getRS(cnseg, by = "region", imput = FALSE, XY = FALSE, what = "mean")

