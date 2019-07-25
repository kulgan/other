library(IRanges)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)


projects = fread("project", h=F)$V1

p = projects[1]

################################################ 
# Comparing broad_significance_results.txt
################################################
bsr = NULL

for(p in projects) {
  old = fread(paste0("./v1/", p, "/broad_significance_results.txt"))
  names(old) = tolower(paste0("old_", gsub("-", "_", gsub(" ", "_", names(old)))))

  new = fread(paste0("./v2/", p, "/broad_significance_results.txt"))
  names(new) = tolower(paste0("new_", gsub("-", "_", gsub(" ", "_", names(new)))))

  both = cbind(old, new[, 3:ncol(new)]) %>% 
    rename(arm = old_arm, 
           ngene = "old_#_genes") %>%
    mutate(amp_z_ratio = new_amp_z_score / old_amp_z_score, 
           del_z_ratio = new_del_z_score / old_del_z_score, 
           amp_q_significant = new_amp_q_value < 0.1 | old_amp_q_value < 0.1, 
           del_q_significant = new_del_q_value < 0.1 | old_del_q_value < 0.1, 
           project = p)

  bsr = rbind(bsr, both)
}

# summarize 
bsr %>% 
  filter(((amp_z_ratio > 2 | amp_z_ratio < 0.5) & amp_q_significant == T) | 
         ((del_z_ratio > 2 | del_z_ratio < 0.5) & del_q_significant == T)) %>%
  select(project, arm, old_del_q_value, new_del_q_value, old_amp_q_value, new_amp_q_value)

#   project arm old_del_q_value new_del_q_value old_amp_q_value new_amp_q_value
# 1    BLCA  2p        1.000000          1.0000           0.465         0.04920
# 2    BRCA  1p        0.016300          0.0651           0.589         0.01990
# 3    LIHC  1p        0.000712          0.0073           0.224         0.00389
# 4    READ  9p        1.000000          1.0000           0.623         0.00920
# 5    READ 20p        1.000000          0.0090           0.000         0.00000
# 6     UVM  8p        0.422000          0.0169           0.000         0.00000


################################################ 
# Comparing genes.conf_99.txt
################################################

gc = NULL

for(p in projects) {
  for(v in c("v1", "v2")) {
    for(t in c("del", "amp")) {
      file = paste0("./", v, "/", p, "/", t, "_genes.conf_99.txt")
      data = data.table(t(fread(file, sep="\t", h=F))) %>% 
        mutate(project = p, 
               type = t, 
               version = v)
      gc = rbind(gc, data)  
    }
  }
}  

gc = data.table(gc)
names(gc)[1:5] = c("cytoband", "q_value", "residual_q_value", "wide_peak_boundaries", "gene")
gc = gc[!(is.na(q_value) | q_value == "q value")]


diff = gc %>% 
  group_by(project, type, gene) %>% 
  summarize(count = n()) %>% filter(count == 1) %>% 
  mutate(list = paste(project, type, gene))

diff.gene = gc %>% 
  filter(paste(project, type, gene) %in% diff$list)

temp = diff.gene %>% select(project, version, type, q_value, gene, cytoband, wide_peak_boundaries) %>% arrange(project, cytoband)





################################################ 
# Comparing all_lesions.conf_99.txt
################################################

lesion = NULL

for(p in projects) {
  for(v in c("v1", "v2")) {
     file = paste0("./", v, "/", p, "/all_lesions.conf_99.txt")
     data = fread(file)
     data = fread(file)[, 1:9] %>% 
       mutate(project = p, 
             version = v)
    lesion = rbind(lesion, data)  
  }
}  
lesion = data.table(lesion)
names(lesion) = c("name", "descriptor", "wide_peak_limits", "peak_limits", "region_limits", 
                  "q_value", "residual_q_value", "broad_or_focal", "amplitude_threshold", 
                  "project", "version")
lesion = lesion[!grepl("CN values", name), -c("amplitude_threshold", "broad_or_focal")] 


# extracting limits
lesion.parsed = lesion %>%
  rowwise() %>%
  mutate(type = strsplit(name, " ")[[1]][1], 
         type = factor(type, levels = c("Deletion", "Amplification")), 
         wide_peak_limits = gsub("^(chr\\w+):(\\d+)-(\\d+)\\(.*", "\\1:\\2:\\3", wide_peak_limits), 
         wpl_chr = strsplit(wide_peak_limits, ":")[[1]][1], 
         wpl_start = as.integer(strsplit(wide_peak_limits, ":")[[1]][2]), 
         wpl_end = as.integer(strsplit(wide_peak_limits, ":")[[1]][3]), 
         peak_limits = gsub("^(chr\\w+):(\\d+)-(\\d+)\\(.*", "\\1:\\2:\\3", peak_limits),  
         pl_chr = strsplit(peak_limits, ":")[[1]][1], 
         pl_start = as.integer(strsplit(peak_limits, ":")[[1]][2]), 
         pl_end = as.integer(strsplit(peak_limits, ":")[[1]][3]), 
         region_limits = gsub("^(chr\\w+):(\\d+)-(\\d+)\\(.*", "\\1:\\2:\\3", region_limits), 
         rl_chr = strsplit(region_limits, ":")[[1]][1], 
         rl_start = as.integer(strsplit(region_limits, ":")[[1]][2]), 
         rl_end = as.integer(strsplit(region_limits, ":")[[1]][3])) %>% 
  select(-c(wide_peak_limits, peak_limits, region_limits)) %>% 
  setDT()


# get chromosome offsets
cytoband = fread("gunzip -c cytoBand.txt.gz")
contigs = paste0("chr", c(1:22, "X"))
chr = cytoband %>% 
  mutate(chr = V1) %>% 
  filter(chr %in% contigs) %>% 
  mutate(chr = factor(chr, levels = contigs)) %>% 
  group_by(chr) %>%
  summarise(len = max(V3)) %>% 
  mutate(offset = 0)
for(i in 2: nrow(chr)) {
  chr$offset[i] = chr$offset[i-1] + chr$len[i-1]
}
chr = chr %>%
  mutate(chr = as.character(chr)) %>% 
  select(chr, offset)

# prepare data input for ggplot, using region_limit ranges
levels = as.vector(t(outer(sort(unique(lesion$project)), c("v1", "v2"), paste, sep=".")))
lesion.grid = lesion.parsed %>% 
  mutate(chr = rl_chr) %>% 
  left_join(chr, by="chr") %>%
  rowwise() %>%
  mutate(start = offset + rl_start,
         end = offset + rl_end, 
         range = paste(seq(floor(start/1e6), ceiling(end/1e6)), collapse=",")) %>% 
  separate_rows(range, sep=",") %>%
  mutate(pos = as.numeric(range), 
         project.version = factor(paste(project, version, sep="."), levels = levels)) %>%
  select(project.version, project, version, type, pos, chr, start, end) 

ggplot(lesion.grid, aes(x = pos, y = reorder(project.version, desc(project.version)))) + 
  geom_raster(aes(fill=type)) + 
  scale_fill_manual(values = c("blue", "red")) +
  labs(x="Position (MB)", y="") +
  theme_bw() + 
  theme(legend.position="top", text = element_text(size=20))


# identify differences between v1 and v2 in grid plot
grid.diff = lesion.grid %>% 
  group_by(project, pos, type, start, end) %>%
  summarise(count = n(), 
            version = paste(version, collapse=" ")) %>% 
  filter(count == 1) %>% 
  select(-count) %>%
  arrange(project, pos, type)

#############################################
# fine limit overlaps between v1 and v2
#############################################
# region limit
lesion.rl = lesion.parsed %>%
  mutate(chr = rl_chr, 
         start = rl_start, 
         end = rl_end) %>% 
  select(project, version, type, chr, start, end, q_value, residual_q_value) %>% 
  group_by(project, version, type, chr, start, end) %>% 
  summarise(q_value = min(q_value), 
            residual_q_value = min(residual_q_value)) %>% 
  mutate(width = end - start + 1, 
         unique = width) %>%
  setDT()

lesion.v1 = lesion.rl[version == "v1"] 
setkey(lesion.v1, project, type, chr, start, end)

lesion.v2 = lesion.rl[version == "v2"] 
setkey(lesion.v2, project, type, chr, start, end)

ov1 = foverlaps(lesion.v1, lesion.v2, type="any", mult="all", which=T, nomatch=0) 
for(i in unique(ov1$xid)) {
  overlaps = lesion.v2[ov1[xid == i]$yid]
  lesion.v1$unique[i] = sum(width(IRanges::setdiff(IRanges(start = lesion.v1$start[i], end = lesion.v1$end[i]), IRanges(start = overlaps$start, end = overlaps$end))))
}
ov2 = foverlaps(lesion.v2, lesion.v1, type="any", mult="all", which=T, nomatch=0) 
for(i in unique(ov2$xid)) {
  overlaps = lesion.v1[ov2[xid == i]$yid]
  lesion.v2$unique[i] = sum(width(IRanges::setdiff(IRanges(start = lesion.v2$start[i], end = lesion.v2$end[i]), IRanges(start = overlaps$start, end = overlaps$end))))
}
lesion.overlap = rbind(lesion.v1, lesion.v2)
lesion.overlap$limit = "region"

# peak limit
lesion.pl = lesion.parsed %>%
  mutate(chr = pl_chr, 
         start = pl_start, 
         end = pl_end) %>% 
  select(project, version, type, chr, start, end, q_value, residual_q_value) %>% 
  group_by(project, version, type, chr, start, end) %>% 
  summarise(q_value = min(q_value), 
            residual_q_value = min(residual_q_value)) %>% 
  mutate(width = end - start + 1, 
         unique = width) %>%
  setDT()

lesion.v1 = lesion.pl[version == "v1"] 
setkey(lesion.v1, project, type, chr, start, end)

lesion.v2 = lesion.pl[version == "v2"] 
setkey(lesion.v2, project, type, chr, start, end)

ov1 = foverlaps(lesion.v1, lesion.v2, type="any", mult="all", which=T, nomatch=0) 
for(i in unique(ov1$xid)) {
  overlaps = lesion.v2[ov1[xid == i]$yid]
  lesion.v1$unique[i] = sum(width(IRanges::setdiff(IRanges(start = lesion.v1$start[i], end = lesion.v1$end[i]), IRanges(start = overlaps$start, end = overlaps$end))))
}
ov2 = foverlaps(lesion.v2, lesion.v1, type="any", mult="all", which=T, nomatch=0) 
for(i in unique(ov2$xid)) {
  overlaps = lesion.v1[ov2[xid == i]$yid]
  lesion.v2$unique[i] = sum(width(setdiff(IRanges::IRanges(start = lesion.v2$start[i], end = lesion.v2$end[i]), IRanges(start = overlaps$start, end = overlaps$end))))
}
temp = rbind(lesion.v1, lesion.v2)
temp$limit = "peak"
lesion.overlap = rbind(lesion.overlap, temp)


# wide peak limit
lesion.wpl = lesion.parsed %>%
  mutate(chr = wpl_chr, 
         start = wpl_start, 
         end = wpl_end) %>% 
  select(project, version, type, chr, start, end, q_value, residual_q_value) %>% 
  group_by(project, version, type, chr, start, end) %>% 
  summarise(q_value = min(q_value), 
            residual_q_value = min(residual_q_value)) %>% 
  mutate(width = end - start + 1, 
         unique = width) %>%
  setDT()

lesion.v1 = lesion.wpl[version == "v1"] 
setkey(lesion.v1, project, type, chr, start, end)

lesion.v2 = lesion.wpl[version == "v2"] 
setkey(lesion.v2, project, type, chr, start, end)

ov1 = foverlaps(lesion.v1, lesion.v2, type="any", mult="all", which=T, nomatch=0) 
for(i in unique(ov1$xid)) {
  overlaps = lesion.v2[ov1[xid == i]$yid]
  lesion.v1$unique[i] = sum(width(IRanges::setdiff(IRanges(start = lesion.v1$start[i], end = lesion.v1$end[i]), IRanges(start = overlaps$start, end = overlaps$end))))
}
ov2 = foverlaps(lesion.v2, lesion.v1, type="any", mult="all", which=T, nomatch=0) 
for(i in unique(ov2$xid)) {
  overlaps = lesion.v1[ov2[xid == i]$yid]
  lesion.v2$unique[i] = sum(width(IRanges::setdiff(IRanges(start = lesion.v2$start[i], end = lesion.v2$end[i]), IRanges(start = overlaps$start, end = overlaps$end))))
}
temp = rbind(lesion.v1, lesion.v2)
temp$limit = "wide_peak"
lesion.overlap = rbind(lesion.overlap, temp)

lesion.overlap = lesion.overlap %>%
  rename(cna_type = type) %>%
  mutate(type = NA, 
         portion = unique/ width,
         type = ifelse(unique == 0, "same", type), 
         type = ifelse(portion == 1, "new", type), 
         type = ifelse(portion < 0.5 & portion > 0 & unique < 1e6, "similar", type), 
         type = ifelse((portion >= 0.5 & portion < 1) | unique >= 1e6, "extended", type)) %>%
  select(-portion) %>%
  setDT()

save(lesion, lesion.parsed, lesion.grid, grid.diff, lesion.overlap, file="lesion.rda")

lesion.overlap.summary = lesion.overlap %>% 
  group_by(limit, version, type) %>% 
  summarise(count = n()) %>%
  mutate(type = factor(type, levels = rev(c("same", "similar", "extended", "new"))))

ggplot(lesion.overlap.summary, aes(x = version, y = count, fill = type)) +
  geom_bar(position = "stack", stat="identity") + 
  facet_wrap( ~ limit, nrow=1)


########################################
# Comparing to Cancer Census Gene
########################################
census = fread("cancer_census_gene_20181017.txt")
names(census) = gsub(" ", "_", names(census))
# Previously confirmd no genes with both amp and del
census = census %>% 
  mutate(cna_type = ifelse(grepl("D", Mutation_Types), "Deletion", NA), 
         cna_type = ifelse(grepl("A", Mutation_Types), "Amplification", cna_type)) %>%
  filter(!is.na(cna_type)) %>% 
  select(ensg, cna_type) %>%
  setDT()

gencode = fread("gencode.v22.genes.txt")
gencode = gencode %>% 
  rowwise() %>% 
  mutate(chr = seqname, 
         ensg = strsplit(gene_id, "\\.")[[1]][1]) %>%
  select(ensg, chr, start, end) %>% 
  left_join(census, by="ensg") %>% 
  setDT() %>%
  setkey(chr, start, end)


segs = lesion.overlap 
segs$index = 1:nrow(segs)

ov = foverlaps(segs, gencode, type="any", mult="all", nomatch=NA) 
ov.summary = ov %>% 
  group_by(index) %>%
  summarise(project = project[1],
            version = version[1],
            type = type[1],
            chr = chr[1],
            start = i.start[1], 
            end = i.end[1],
            ngene = sum(!is.na(ensg)),
            ncensus = sum(cna_type == i.cna_type, na.rm=T), 
            limit = limit[1], 
            q_value = q_value[1],
            residual_q_value = residual_q_value[1]) %>% 
  setDT()



            

ov.summary %>% 
  group_by(limit, version, type) %>% 
  summarise(all = sum(ngene), census = sum(ncensus))



temp = ov.summary %>% 
  group_by(version, limit, type) %>% 
  summarise(nseg = n(), 
            nseg_census = sum(ncensus > 0), 
            average_seg_length = mean(end - start + 1), 
            ngene = sum(ngene), 
            ncensus = sum(ncensus)) %>% 
  setDT()


new.v2 = ov.summary[ncensus > 0 & version == "v2" & type == "new"] %>% setkey(project, chr, start, end)
ov2 = foverlaps(segs, new.v2)


r = sum(temp$all) / sum(temp$census)
temp$ratio = temp$census / temp$all * r
temp = data.table(temp)


diff = gc %>% 
  group_by(project, type, gene) %>% 
  summarize(count = n()) %>% filter(count == 1) %>% 
  mutate(list = paste(project, type, gene))

diff.gene = gc %>% 
  filter(paste(project, type, gene) %in% diff$list)

temp = diff.gene %>% select(project, version, type, q_value, gene, cytoband, wide_peak_boundaries) %>% arrange(project, cytoband)
































































