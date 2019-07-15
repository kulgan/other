
ReadAndFilterSegmentationFile = function(file, contigs, header = T) {
  # permissible values
  chr.contigs = paste0("chr", contigs)

  # read file
  tsv = fread(file, h = header)

  setnames(tsv, c("sample", "chr", "start", "end", "nprobes", "segmean"))

  # reformat Chromosome names and sort
  tsv = tsv %>%  
    filter(chr %in% c(contigs, chr.contigs)) %>% 
    mutate( chr = gsub("chr", "", as.character(chr)), 
        chr = replace(chr, chr == "23", "X"), 
        chr = replace(chr, chr == "24", "Y"), 
        chr = factor(chr, levels = contigs)) %>%
    arrange(chr, start, end) %>% 
    setDT()

  return(tsv)
}


CentromereToArm <- function(chr.info) {

  contigs <- levels(chr.info$chr)
  all.arms <- c(t(outer(contigs, c("p", "q"), FUN = paste0)))

  p.arms <- chr.info %>% 
    mutate(start = 1, 
      end = centromere, 
      arm = paste0(chr, "p")) %>% 
    select(chr, start, end, arm)

  q.arms <- chr.info %>% 
    mutate(start = centromere + 1, 
      end = chr_len, 
      arm = paste0(chr, "q")) %>% 
    select(chr, start, end, arm)

  arms <- rbind(p.arms, q.arms) %>%
    mutate(arm = factor(arm, levels = all.arms)) %>%
    setDT() %>% 
    setkey(., chr, start, end)

  return(arms)
}

CleanUpSegmentations = function(segs, chr.arms) {
  # find sample median of autosomal segmeans  
  sample.segmean.median = segs %>% 
    filter(chr %in% c(1:22)) %>%
    group_by(sample) %>% 
    summarise(sample_segmean_median = weightedMedian(segmean, end - start + 1, interpolate = T))

  # normalize segmeans by sample median of autosomal segmeans
  normlized.segs = segs %>% 
    left_join(., sample.segmean.median, by = "sample") %>% 
    mutate(segmean = segmean - sample_segmean_median) %>% 
    select(-sample_segmean_median) %>% 
    setDT()

  # split segments by centromere and sort
  splitted.segs = SplitSegByCentromere(normlized.segs, chr.arms) %>%
    arrange(sample, chr, start, end) %>% 
    setDT()

  return(splitted.segs)
}

SegmeanToCN = function(segmean, cap = 1.5) {
  # limit max and min segmean by cap, can calculate CN from Segmean
  segmean = replace(segmean, segmean > cap, cap)
  segmean = replace(segmean, segmean < -cap, -cap)
  cn = 2^(segmean + 1)
  return(cn)
}


ReadGeneFile = function(file, header = T) {
  # read file
  tsv = fread(file, h = header)
  setnames(tsv, c("chr", "start", "end", "strand", "gene"))
  return(tsv)
}

AnnotateGenes = function(genes, contigs, chr.arms) {

  ov = foverlaps(chr.arms, genes, mult="all", type="any", nomatch=0)

  annotated.genes <- ov %>% 
    rename(arm_start = i.start, 
           arm_end = i.end) %>% 
    arrange(chr, start, end) %>% 
    setDT() %>% 
    setkey(chr, start, end) 


}



ReadAndFilterGeneFile = function(file, contigs, header = T) {
  # permissible values
  chr.contigs = paste0("chr", contigs)

  # read file
  tsv = fread(file, h = header)

  setnames(tsv, c("chr", "start", "end", "strand", "gene"))

  # reformat Chromosome names and sort
  tsv = tsv %>%  
    filter(chr %in% c(contigs, chr.contigs)) %>% 
    mutate( chr = gsub("chr", "", as.character(chr)), 
        chr = replace(chr, chr == "23", "X"), 
        chr = replace(chr, chr == "24", "Y"), 
        chr = factor(chr, levels = contigs)) %>%
    arrange(chr, start, end) %>% 
    setDT() %>%
    setkey(chr, start, end)

  return(tsv)
}



# ReadGenomicSegment read a TSV file, and return a data.table of file content
# if header = T, the following 3 column names must exist: "Chromosome", "Start" and "End"; or the first 3 columns will be named so
# the function will also remove "chr" string from the contig names, and change 23 to X, 24 to Y
ReadGenomicSegment = function(file, header = T, keep.chry = F) {
  # permissible values
  contigs = c(1:24, "X", "Y")
  chr.contigs = paste0("chr", contigs)

  # read file
  tsv = fread(file, h = header)

  if (!header) {
    colnames(tsv)[1:3] = c("Chromosome", "Start", "End")
  }

  tsv = tsv %>%  
    filter(Chromosome %in% c(contigs, chr.contigs)) %>% 
    mutate( Chromosome = gsub("chr", "", Chromosome), 
        Chromosome = replace(Chromosome, Chromosome == "23", "X"), 
        Chromosome = replace(Chromosome, Chromosome == "24", "Y"), 
        Chromosome.Factor = factor(Chromosome, levels = contigs)) %>%
    arrange(Chromosome.Factor, Start, End) %>% 
    select(-Chromosome.Factor) %>% 
    setDT()

  # remove Y chromosome is needed
  if(keep.chry) {
    return(tsv) 
  } else {
    return(tsv[Chromosome != "Y"])
  }
}




ReadAndFilterChrInfo <- function(file, contigs, header = T) {
  # permissible values
  chr.contigs <- paste0("chr", contigs)

  # read file
  tsv <- fread(file, h = header)

  setnames(tsv, c("chr", "centromere", "chr_len"))

  # reformat Chromosome names and sort
  tsv <- tsv %>%  
    filter(chr %in% c(contigs, chr.contigs)) %>% 
    mutate( chr = gsub("chr", "", as.character(chr)), 
        chr = replace(chr, chr == "23", "X"), 
        chr = replace(chr, chr == "24", "Y"), 
        chr = factor(chr, levels = contigs)) %>%
    arrange(chr) %>% 
    setDT()

  return(tsv)
}



# ReadGenomicLocation read a TSV file, and return a data.table of file content
# The first 2 columns must be "Chromosome" and "Position"
# the function will also remove "chr" string from the contig names, and change 23 to X, 24 to Y
ReadAndFilterCentromereFile = function(file, contigs, header = T) {
  # permissible values
  chr.contigs = paste0("chr", contigs)

  # read file
  tsv = fread(file, h = header)

  names(tsv)[1:2] = c("chr", "centromere")

  # reformat Chromosome names and sort
  tsv = tsv %>%  
    filter(chr %in% c(contigs, chr.contigs)) %>% 
    mutate( chr = gsub("chr", "", as.character(chr)), 
        chr = replace(chr, chr == "23", "X"), 
        chr = replace(chr, chr == "24", "Y"), 
        chr_factor = factor(chr, levels = contigs)) %>%
    arrange(chr_factor) %>% 
    select(-chr_factor) %>% 
    setDT()

  return(tsv)
}



NormalizeSegmean = function(segs) {
  # only based on autosomes
  contigs = c(1:22) 

  # calculate weighted median 
  weighted.median = with(segs[Chromosome %in% contigs], weightedMedian(Segment_Mean, Num_Probes))
  segs$Segment_Mean = segs$Segment_Mean - weighted.median
  return(segs)
}


# GetArm will return arm values of segment in GenomicTSV data.table based on the centromere position provided
GetArm = function(genomictsv, centromere) {
  # get arm information of each segment and return
  arm.info = genomictsv %>% 
    left_join(., centromere, by = "chr") %>% 
    mutate(arm = ifelse(start > centromere, paste0(chr,"q"), paste0(chr, "p"))) 
  return(arm.info$arm)  
}


# SplitSegByCentromere will split segment by centromere, and adjust Num_Probes and Star or End values
SplitSegByCentromere <- function(segs, chr.arms) {

  ov <- foverlaps(segs, chr.arms, type = "any", mult = "all", nomatch = 0)
  split.segs <- ov %>% 
    mutate( start = pmax(start, i.start), 
            end = pmin(end, i.end), 
            nprobes = ceiling((end - start) / (i.end - i.start) * nprobes)) %>% 
    select(sample, chr, start, end, nprobes, segmean, arm) %>% 
    setDT()

  return(split.segs)  
} 

# GetArmSegmean returns a data.table of chromosome arms and their weighted-average segmeans
GetArmSegmean = function(segs, contigs, min.arm.evidence = 100, centromere, bed = NULL, min.bed.segment.len = 50) {
  # get all possible arms
  arms = c(t(outer(contigs, c("p", "q"), FUN = paste0)))

  sample.arm.combinations = expand.grid(sample = unique(segs$sample), arm = arms, stringsAsFactors = F)

  # calcualte arm-level segment means
  # if number of probes on an arm is less than min.arm.evidence, set segmean = 0
  arm.summary = segs %>% 
    mutate(arm = GetArm(., centromere)) %>% 
    group_by(sample, arm) %>%
    summarise(  arm_nprobes = sum(nprobes), 
          median_arm_segmean = weightedMedian(segmean, nprobes, interpolate = T)) %>% 
  # mutate(median_arm_segmean = ifelse(arm_nprobes < min.arm.evidence, 0, median_arm_segmean)) %>% 
    right_join(., sample.arm.combinations, by = c("sample", "arm")) %>% 
    mutate( median_arm_segmean = replace(median_arm_segmean, is.na(median_arm_segmean), 0) ,
        enough_arm_support = ! (is.na(arm_nprobes) | arm_nprobes < min.arm.evidence)) %>% 
    mutate(arm_factor = factor(arm, levels = arms)) %>% 
    arrange(sample, arm_factor) %>% 
    select(-c(arm_factor))

  return(arm.summary)
}

# GetArmSegmean returns a data.table of chromosome arms and their weighted-average segmeans
GetArmCopyNumber = function(segs, contigs, min.arm.evidence = 100, bed = NULL, min.bed.segment.len = 50) {
  # get all possible arms
  # arms = c(t(outer(contigs, c("p", "q"), FUN = paste0)))
  # sample.arm.combinations = expand.grid(sample = unique(segs$sample), arm = arms, stringsAsFactors = F)

  # calcualte arm-level segment means
  # if number of probes on an arm is less than min.arm.evidence, set segmean = 0
  arm.summary = segs %>% 
    group_by(sample, arm) %>%
    summarise(arm_nprobes = sum(nprobes), 
              median_arm_copy_number = weightedMedian(copy_number, nprobes, interpolate = T)) %>% 
    complete(arm, fill = list(arm_nprobes = 0, median_arm_copy_number = 2)) %>% 
    mutate(enough_arm_support = arm_nprobes >= min.arm.evidence) %>% 
    arrange(sample, arm) 


  # if number of bed segment on an arm is less than min.arm.evidence, set segmean = 0
  # here, we assume bed file is in good shape (not cover centromeres)
  if(!is.null(bed)) {
    # annotate bed by arm
    bed.arm.summary = bed %>% 
      mutate(arm = GetArm(., centromere)) %>% 
      filter(end - start >= min.bed.segment.len) %>%
      group_by(arm) %>%
      summarise(nsegment = length(arm))

    # set nsegment < min.arm.evidence arms to have segmean = 0  
    arm.summary[arm %in% bed.arm.summary$arm[which(bed.arm.summary$nsegment < min.arm.evidence)]]$segmean = 0
  }

  return(arm.summary)
}


# GetCutoffs return a new cutoffs array that has modified values
# the 1st element of cutoffs: low_cutoff = min(segmean_of_each_arm, default_low_cutoff)
# the 4th element of cutoffs: high_cutoff = max(segmean_of_each_arm, default_high_cutoff) 
GetSampleCopyNumber = function(arm.cn) {
  sample.summary = arm.cn %>% 
    filter(enough_arm_support) %>% 
    group_by(sample) %>%
    summarise(low = min(median_arm_copy_number),
        high = max(median_arm_copy_number))
  return(sample.summary)
}


GetSampleSegmean = function(arm.segmean) {
  sample.summary = arm.segmean %>% 
    filter(enough_arm_support) %>% 
    group_by(sample) %>%
    summarise(low = min(median_arm_segmean),
        high = max(median_arm_segmean))
  return(sample.summary)
}



GetSampleThreshold = function(sample.cn, noise) {
  # assign -2/-1/0/1/2 to each gene
  cutoffs = sample.cn %>%
    mutate(t1 = low - noise, 
           t2 = 2 - noise, 
           t3 = 2 + noise,
           t4 = high + noise) %>%
    setDT()
  return(cutoffs)
}


GetThresholedValue = function(gene.cn, sample.th) {
  # assign -2/-1/0/1/2 to each segmean
  gene.score = gene.cn %>% 
    left_join(sample.th, by="sample") %>% 
    mutate(score = 0, 
           score = replace(score, copy_number > t3, 1),
           score = replace(score, copy_number > t4, 2), 
           score = replace(score, copy_number < t2, -1), 
           score = replace(score, copy_number < t1, -2))
  return(gene.score$score)
}




# AddFilter add a new filter to existing filter string
# if existing filter is "PASS", replaced with new filter
# else, append new filter to the end, and separate with semicolon
AddFilter = function(old, new) {
  return(ifelse(old == "PASS", new, paste(old, new, sep=";")))
}

# GetGemeSegmean find all potential segmeans for each gene:
# 1. overlapped segment segmeans
# 2. arm-level segmean if no overlapped segments or the only overlapped segment is mall (less than min.seg.evidence probe support)
# 3. segmean of left and right neighbor segments (within neighbor.distance) if no overlapped segments
# It returns a data.frame of 3 columns Gene, potential_segmean and source.  
# "source" indicates the source of the segmean: overlap, arm, left_neighbor, right_neighbor. 
# GetGeneSegmean = function(segs, arm.segmean, gene, neighbor.distance = 1e6, min.seg.evidence = 4) {
  # sort segs, define small_segment, and setkey for later foverlaps 
  segs = segs %>% 
    mutate(is_small_segment = Num_Probes < min.seg.evidence) %>% 
    rename(seg_segmean = Segment_Mean) %>%
    select(Chromosome, Start, End, seg_segmean, is_small_segment) %>%
    arrange(Chromosome, Start) %>% 
    setDT() %>% 
    setkey(., Chromosome, Start, End)


  # get arm segmean of each gene  
  gene.info = gene %>% 
    mutate(arm = GetArm(., centromere)) %>% 
    inner_join(., arm.segmean, by="arm") %>%
    rename(arm_segmean = segmean) %>% 
    select(Chromosome, Start, End, Gene, arm_segmean) %>% 
    setDT()

  # overlap between gene and segment
  ov = foverlaps(gene.info, segs, type = "any", mult = "all", nomatch = 0) 

  # get seg_segmean as potential_segmean for genes overlapped with segments
  gene.segmean = ov %>% 
    mutate( potential_segmean = seg_segmean, 
        source = "overlap") %>%
    select(Chromosome, i.Start, i.End, Gene, potential_segmean, source) %>%
    rename(End = i.End, Start = i.Start)

  # get arm_segmean as potential_segmean for genes not overlaps with any segments
  gene.info.in_gap = gene.info[!Gene %in% ov$Gene] 
  gene.segmean0 = gene.info.in_gap %>%
    mutate( potential_segmean = arm_segmean, 
        source = "arm") %>% 
    select(Chromosome, Start, End, Gene, potential_segmean, source)

  # get left neighbor segmean as potential_segmean for genes not overlaps with any segments
  gene.segmean0.left = gene.info.in_gap %>% 
    mutate(Start = Start - neighbor.distance) %>% 
    setDT() %>%
    foverlaps(., segs, type = "any", mult = "last", nomatch = 0) %>% 
    mutate( potential_segmean = seg_segmean, 
        source = "left_neighbor") %>% 
    select(Chromosome, Start, End, Gene, potential_segmean, source)

  # get right neighbor segmean as potential_segmean for genes not overlaps with any segments    
  gene.segmean0.right = gene.info.in_gap %>% 
    mutate(End = End + neighbor.distance) %>% 
    setDT() %>%
    foverlaps(., segs, type = "any", mult = "first", nomatch = 0) %>% 
    mutate( potential_segmean = seg_segmean, 
        source = "right_neighbor") %>% 
    select(Chromosome, Start, End, Gene, potential_segmean, source) 

  # get arm_segmean as potential_segmean for genes that only overlaps with a small segment
  count = table(ov$Gene)
  gene.segmean1 = ov %>% 
    filter(Gene %in% names(count)[count == 1] & is_small_segment) %>% 
    mutate( potential_segmean = arm_segmean, 
        source = "arm") %>% 
    select(Chromosome, Start, End, Gene, potential_segmean, source)

  # merge all potential segmeans. (this is the most important intermediate dataset) 
  gene.segmean = rbind(gene.segmean, gene.segmean0, gene.segmean0.left, gene.segmean0.right, gene.segmean1) 

  return(gene.segmean)
  
}





GetGeneSegmentOverlaps = function(segs, arm.cn, gene, min.seg.evidence = 4) {

  # sort segs, define small_segment, and setkey for later foverlaps 
  segs = segs %>% 
    mutate(is_small_segment = nprobes < min.seg.evidence) %>% 
    select(sample, chr, start, end, copy_number, is_small_segment) %>%
    setDT() %>% 
    setkey(., sample, chr, start, end)


  # get arm segmean of each gene  
  gene.arm.info = genes %>% 
    inner_join(., arm.cn, by="arm") %>%
    # rename(arm_segmean = segmean) %>% 
    # select(sample, chr, start, end, gene, median_arm_copy_number, enough_arm_support, centromere) %>% 
    setDT()

  # overlap between gene and segment
  ov = foverlaps(gene.arm.info, segs, type = "any", mult = "all", nomatch = NA) %>%
    rename( gene_start = i.start, 
        gene_end = i.end) %>% 
    setDT()



  # get segment copy_number as copy_number for genes overlapped with segments
  gene.overlap = ov[!is.na(start)] %>% 
    mutate( overlap = pmin(end, gene_end) - pmax(start, gene_start) + 1,
            type = "overlap") %>%
    select(sample, gene, copy_number, overlap, is_small_segment, enough_arm_support, median_arm_copy_number, type)



  # add median_arm_copy_number as copy_number for genes not overlaps with any segments
  gene.in_gap = ov[is.na(start)]

  gene.in_gap.left = gene.in_gap %>% 
    mutate(start = arm_start, 
           end = gene_start) %>% 
    select(sample, chr, start, end, gene, gene_start, gene_end, median_arm_copy_number, enough_arm_support) %>% 
    setDT() %>%
    foverlaps(., segs, type = "any", mult = "last", nomatch = NA) %>% 
    mutate(overlap = gene_start - end, 
           type = ifelse(is.na(overlap), "arm", "left_neighbor"), 
           copy_number = ifelse(is.na(overlap), median_arm_copy_number, copy_number)) %>%  
    select(sample, gene, copy_number, overlap, is_small_segment, enough_arm_support, median_arm_copy_number, type) %>%
    setDT()

  gene.in_gap.right = gene.in_gap %>% 
    mutate(start = gene_end, 
           end = arm_end) %>% 
    select(sample, chr, start, end, gene, gene_start, gene_end, median_arm_copy_number, enough_arm_support) %>% 
    setDT() %>%
    foverlaps(., segs, type = "any", mult = "last", nomatch = NA) %>% 
    mutate(overlap = start - gene_end, 
          type = ifelse(is.na(overlap), "arm", "right_neighbor"), 
          copy_number = ifelse(is.na(overlap), median_arm_copy_number, copy_number)) %>%  
    select(sample, gene, copy_number, overlap, is_small_segment, enough_arm_support, median_arm_copy_number, type) %>%
    setDT()


  gene.overlap = data.table(rbind(gene.overlap, gene.in_gap.right, gene.in_gap.left))

  return(gene.overlap)
}








AggregateCopyNumber = function(gene.cn, sample.cn, noise) {
  gene.cn$score = GetThresholedValue(gene.overlap, sample.th)

  gene.overlap = gene.cn[type=="overlap"]
  gene.overlap %>% 
    group_by(sample, gene) %>%
    filter(n() == 1)
}



# SelectOneValue accept an array of numeric values, and return only one of them based on method
# In the default settings, it returns the item with the largest absolute value, and favor negative number when tie
# The function is not used for performance reason 
SelectOneValue = function(arr, method = "extreme", favor.positive = F) {
  # sort arr first to favor selection of positive or negative values with the same absolute value
  arr = sort(arr, decreasing = favor.positive)
  return(arr[which.max(abs(arr))])
}
 
# AggregateSegmean select only one segmean if multiple exist for each gene
# The default method is "extreme", and favor negative values
# filter values are also added to indicate the reliability of the segmeans
AggregateSegmean = function(gene.segmean = gene.segmean, cutoffs = cutoffs) {

  setDT(gene.segmean)
  imputed.genes = setdiff(gene.segmean[source == "arm"]$Gene, gene.segmean[source == "overlap"]$Gene)

  # aggregate, get min and max, find extreme value, and thresholded to integer scores by cutoffs
  score = gene.segmean %>%
    group_by(Gene) %>% 
    summarise(  Chromosome = Chromosome[1], 
          Start = Start[1], 
          End = End[1],
          min_segmean = min(potential_segmean), 
          max_segmean = max(potential_segmean)) %>% 
    mutate( min_score = GetThresholedValue(min_segmean, cutoffs), 
        max_score = GetThresholedValue(max_segmean, cutoffs), 
        similar = abs(min_score - max_score) <= 1, 
        one_segmean = ifelse(max_segmean > 0 & abs(max_segmean) > abs(min_segmean), max_segmean, min_segmean), 
        one_score = GetThresholedValue(one_segmean, cutoffs), 
        filter = "PASS", 
        filter = ifelse(Gene %in% imputed.genes, AddFilter(filter, "impute"), filter), 
        filter = ifelse(similar, filter, AddFilter(filter, "non_det"))) %>%
    select(Chromosome, Start, End, Gene, one_segmean, one_score, filter)

  return(data.table(score))
}

# add off-target filter to gene.score
FilterByBed = function(gene.score, bed) {
  # normalize bed coordinates with 1-indexed interval coordinates
  bed$Start = bed$Start + 1
  # overlap between gene and bed
  setkey(bed, Chromosome, Start, End)
  ov = foverlaps(gene.score, bed, type = "any", mult = "first", nomatch = NA, which = T) 
  
  # get genes without overlaps, and add filter "off_target"
  off.target.index = which(is.na(ov))
  gene.score$filter[off.target.index] = AddFilter(gene.score$filter[off.target.index], "off_target")

  return(gene.score)
}

# reformat and write  
WriteScores = function(gene.score, gene, filters, out.file) {
  # re-oder genes and change colnames
  gene.score = gene.score[match(gene$Gene, gene.score$Gene)] %>%
    mutate(Chromosome = paste0("chr", Chromosome)) %>%
    rename( Segment_Mean = one_segmean, 
        Copy_Number_Score = one_score, 
        Filter = filter)

  # write date
  date = format(Sys.Date(), "%Y%m%d")
  cat("##fileDate=", date, "\n", sep ="", file = out.file)

  # write filters 
  for(i in 1:nrow(filters)) {
    cat("##Filter=<ID=", filters$ID[i], ",Description=\"", filters$Description[i], "\">", "\n", sep = "", file = out.file, append = T)
  }

  # write content
  suppressWarnings(write.table(gene.score, file = out.file, append = T, col.names = T, row.names = F, sep="\t", quote=F))
}


