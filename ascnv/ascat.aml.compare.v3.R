# looking for genotype == AB and LOH (nb = 0 & na > 0)
library(data.table)
library(dplyr)

# define chrs
chrs <- paste0("chr", 1:22)

# read metadata
meta <- fread("~/SCRATCH/GDC.SNP6.pairs.10878.02192024.tsv")
cels <- fread("~/SCRATCH/GDC.CEL.file.tsv")
my.meta <- meta %>% 
  filter(project %in% c("TCGA-LAML", "TARGET-AML")) %>% 
  rename(type = project, aliquot = tumor_aliquot) %>% 
  setDT()
my.meta$normal_cel_filename <- cels$filename[match(my.meta$normal_cel, cels$cel)]
my.meta$tumor_cel_filename <- cels$filename[match(my.meta$tumor_cel, cels$cel)]


# read probe location information
probe <- fread("zcat ~/SCRATCH/snp6.na35.remap.hg38.txt.gz")
probe <- probe %>% 
  filter(chr %in% chrs & grepl("^SNP_A-", probeid)) %>%
  mutate(start = pos, 
         end = pos) %>% 
  select(probeid, chr, start, end, pos) %>% 
  setDT() %>%
  setkey(chr, start, end)

# read segment
seg <- fread("zcat ~/SCRATCH/TCGA.TARGET.SNP6.ASCAT.11582.1_off_modified.segment.txt.gz")
names(seg) <- c("aliquot", "chr", "start", "end", "na", "nb", "cn")
my.seg <- seg %>%
  filter(aliquot %in% my.meta$aliquot & chr %in% chrs) %>% 
  merge(my.meta[, c("aliquot", "type"), with=F], by = "aliquot") %>% 
  mutate(chr = factor(chr, levels=chrs)) %>%
  setDT() %>%
  setkey(chr, start, end)


# find overlapped probes with segments
ov <- foverlaps(my.seg, probe, mult="all", type="any", nomatch=0L) 
names(ov)[6:7] <- c("seg_start", "seg_end")


# find boundaries
bound <- melt(my.seg[, c("chr", "start", "end")], id.var="chr", variable.name="type", value.name="position") %>%
  distinct(chr, type, position)
bound2 <- bound %>% 
  mutate(position = ifelse(type == "start", position - 1, position + 1), 
         type = ifelse(type == "start", "end", "start"))  
bound3 <- data.frame(chr = chrs, type = "start", position = 1)
bound4 <- data.frame(chr = chrs, type = "end", position = 1e9)
temp <- rbind(bound, bound2, bound3, bound4) %>% 
  arrange(chr, position) %>% 
  distinct(chr, type, position) 
w <- seq(1, nrow(temp), 2)
bound <- data.table(chr = factor(temp$chr[w], levels = chrs), 
                   start = temp$position[w], 
                   end = temp$position[-w], 
                   index = 1:length(w)) 
setkey(bound, chr, start, end)

ov <- foverlaps(probe, bound, mult="all", nomatch=0L) 
bound <- ov %>% 
  group_by(index) %>% 
  summarise(chr = chr[1], 
          start = min(pos), 
          end = max(pos)) %>% 
  mutate(chr = factor(chr, levels = chrs)) %>%
  arrange(chr, start, end) %>%
  mutate(index = 1:n()) %>%
  setDT() %>%
  setkey(chr, start, end)
fwrite(bound, "aml.boundaries.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# overlap with segment
my.seg$loh <- my.seg$nb == 0 & my.seg$na > 0
ov <- foverlaps(my.seg, bound, mult="all", type="any", nomatch=0L)

block.summary <- ov%>%
  group_by(index) %>%
  summarise(adult.loh = sum(type == "TCGA-LAML" & loh==T), 
            adult.other = sum(type == "TCGA-LAML" & loh==F), 
            kid.loh = sum(type == "TARGET-AML" & loh==T), 
            kid.other = sum(type == "TARGET-AML" & loh==F)) %>%
  setDT()


block.summary$pval <- apply(block.summary[, c("adult.loh", "adult.other", "kid.loh", "kid.other"), with=F], 1, function(x) {
                  dim(x) = c(2, 2); 
                  fisher.test(x, conf.int = F, alternative = "two.sided")$p.value})

block.summary$pval[which(block.summary$pval > 1)] = 1
qobj = qvalue(block.summary$pval)
block.summary$qval = qobj$qvalues
block.summary$lfdr = qobj$lfdr
fwrite(block.summary, "block.summary.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# potential block with FDR <= 0.1
block <- block.summary %>%
  filter(qval <= 0.1) %>%
  left_join(bound, by="index") %>% 
  rename(block.index = index) %>%
  setDT() %>%
  setkey(chr, start, end)
fwrite(block, "block.fdr01.tsv", col.names=T, row.names=F, sep="\t", quote=F)

ov <- foverlaps(probe, block, mult="all", type="any", nomatch=0)
my.probe <- ov %>%
  rename(seg.start = start,
         seg.end = end) %>%
  select(-c("i.start", "i.end")) %>%
  setDT()
fwrite(my.probe, "snp.fdr01.tsv", col.names=T, row.names=F, sep="\t", quote=F)


# read and subset genotypes
geno <- fread("zcat ~/SCRATCH/genotype/AML2/consensus.call.tsv.gz")

common.probes = intersect(my.probe$probeid, geno$probeid)
my.probe = my.probe[match(common.probes, probeid)]
my.geno <- geno[match(common.probes, probeid)]

m <- match(my.meta$normal_cel_filename, names(geno))
my.geno <- my.geno[, m, with=F]

w <- which(apply(my.geno, 1, function(x) sum(x==1)) < 5)
my.geno <- my.geno[-w]
my.probe <- my.probe[-w]
common.probes <- common.probes[-w]

# read and subset baf
baf.probes <- fread("~/SCRATCH/baf/probes")
names(baf.probes) <- "probeid"
baf.files <- fread("~/SCRATCH/baf/baf.files", h=F)$V1
w <- which(baf.probes$probeid %in% common.probes)
baf.file.line <- data.table(baf.file = baf.files[floor(w/ 5000) + 1], 
                           line = w %% 5000 + 1)
fwrite(baf.file.line, "~/SCRATCH/baf/baf.file.line", col.names=F, row.names=F, sep="\t", quote=F)
# while read -r a b; do sed -ne ''"$b"'p' $a >> selected.baf; done< "baf.file.line"
baf <- fread("~/SCRATCH/baf/selected.baf", h=F)
baf.header <- fread("~/SCRATCH/baf/tcga_snp6_baf_header")
names(baf) <- names(baf.header)
baf <- baf[match(common.probes, baf$probeid), ]
normal.baf <- baf[, match(my.meta$normal_cel, names(baf)), with=F]
tumor.baf <- baf[, match(my.meta$tumor_cel, names(baf)), with=F]

# data aggregation
names(my.geno) <- my.meta$patient
my.geno$probeid <- common.probes
geno.long <- melt(my.geno, id.var="probeid", variable.name="patient", value.name="normal.genotype")
names(normal.baf) <- my.meta$patient
normal.baf$probeid <- common.probes
normal.baf.long <- melt(normal.baf, id.var="probeid", variable.name="patient", value.name="normal.baf")
names(tumor.baf) <- my.meta$patient
tumor.baf$probeid <- common.probes
tumor.baf.long <- melt(tumor.baf, id.var="probeid", variable.name="patient", value.name="tumor.baf")

d <- geno.long
d$normal.baf <- normal.baf.long$normal.baf
d$tumor.baf <- tumor.baf.long$tumor.baf 
d$baf.diff <- d$tumor.baf - d$normal.baf
m = match(d$patient, my.meta$patient)
d$type <- my.meta$type[m]
d$gender <- my.meta$gender[m]
d$purity <- my.meta$ascat_purity[m]
d$ploidy <- my.meta$ascat_ploidy[m]

d <- d %>% 
  inner_join(my.probe, by="probeid") %>%
  select(-c("block.index", "adult.loh", "adult.other", "kid.loh", "kid.other", "pval", "qval", "lfdr", "seg.start", "seg.end")) %>%
  mutate(start = pos) %>%
  rename(end = pos) %>%
  setDT()

seg.info = my.seg %>% 
  inner_join(my.meta, by="aliquot") %>% 
  select("patient", "chr", "start", "end", "loh", "na", "nb") %>%
  setDT() %>%
  setkey(patient, chr, start, end)

ov <- foverlaps(d, seg.info, by.x=c("patient", "chr", "start", "end"), by.y=c("patient", "chr", "start", "end"), mult="all", type="any", nomatch=0L)
d <- ov %>%
  rename(position = i.start) %>%
  select(-c("start", "end", i.end)) %>%
  mutate(chr = factor(chr, levels = chrs)) %>%
  arrange(chr, position) %>%
  setDT() 

fwrite(d, "aml.loh.258case.260probe.tsv", col.names=T, row.names=F, sep="\t", quote=F)

d2 <- rbind(d %>% 
            select(-c("baf.diff", "normal.baf")) %>%
            rename(baf = tumor.baf) %>%
            mutate(tissue = "tumor"), 
            d %>% 
            select(-c("baf.diff", "tumor.baf")) %>%
            rename(baf = normal.baf) %>%
            mutate(tissue = "normal")) %>%
      setDT()


ggplot(d2[type=="TARGET-AML" & normal.genotype==1 & loh==T], aes(baf, col=tissue)) + 
  geom_density() + 
  facet_wrap( ~ probeid, nrow=2)

d3 <- d2[type=="TARGET-AML" & normal.genotype==1 & loh==T]

temp <- d3 %>% 
  group_by(probeid) %>% 
  summarise(target_count = n())

d3 <- d3 %>%
  left_join(temp, by="probeid") %>%
  mutate(key = paste0(probeid, "\n", target_count, " patients")) %>% 
  setDT() 

d4 = d3[probeid %in% head(ps)]

ggplot(d3[type=="TARGET-AML" & normal.genotype==1], aes(baf, col=tissue)) + 
  geom_density() + 
  facet_wrap( ~ key, nrow=10)



