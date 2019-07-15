
ReadAndFilterSegmentation <- function(file, contigs, header = T) {
  # permissible values
  chr.contigs <- paste0("chr", contigs)

  # read file and set names
  tsv <- fread(file, h = header)
  setnames(tsv, c("sample", "chr", "start", "end", "nprobes", "segmean"))

  # reformat chromosome names, and set keys
  tsv <- tsv %>%  
    filter(chr %in% c(contigs, chr.contigs)) %>% 
    mutate( chr = gsub("chr", "", as.character(chr)), 
        chr = replace(chr, chr == "23", "X"), 
        chr = replace(chr, chr == "24", "Y"), 
        chr = factor(chr, levels = contigs)) %>%
    setDT() %>%
    setkey(chr, start, end)

  # return segmentations as data.table
  return(tsv)
}


ReadAndFilterArmInfo <- function(file, contigs, header = T) {
  # permissible values
  chr.contigs <- paste0("chr", contigs)
  arms = c(t(outer(contigs, c("p", "q"), FUN = paste0)))


  # read file and setnames
  tsv <- fread(file, h = header)
  setnames(tsv, c("chr", "arm", "start", "end"))

  # reformat chromosome names and setkeys
  tsv <- tsv %>%  
    filter(chr %in% c(contigs, chr.contigs)) %>% 
    mutate( chr = gsub("chr", "", as.character(chr)), 
        chr = replace(chr, chr == "23", "X"), 
        chr = replace(chr, chr == "24", "Y"), 
        chr = factor(chr, levels = contigs),
        arm = factor(arm, levels = arms)) %>%
    setDT() %>% 
    setkey(chr, start, end)

  # return chromosome arm information as data.table 
  return(tsv)
}




CleanUpSegmentations <- function(segs, arm.info, min.seg.evidence) {
  # find sample median of autosomal segmeans  
  sample.segmean.median <- segs %>% 
    filter(chr %in% c(1:22)) %>%
    group_by(sample) %>% 
    summarise(median = weightedMedian(segmean, end - start + 1, interpolate = T))

  # normalize segmeans by sample median of autosomal segmeans
  normalized.segs <- segs %>% 
    left_join(., sample.segmean.median, by = "sample") %>% 
    mutate(segmean = segmean - median) %>% 
    setDT()

  # split segments by arms, and define dubious segment as nprobes < min.seg.evidence
  splitted.segs <- foverlaps(normalized.segs, chr.arms, type = "any", mult = "all", nomatch = 0) %>% 
    mutate(start = pmax(start, i.start), 
           end = pmin(end, i.end), 
           nprobes = ceiling((end - start) / (i.end - i.start) * nprobes)) %>%
    select(sample, chr, start, end, nprobes, segmean, arm) %>% 
    arrange(sample, chr, start, end) %>%
    setDT() %>%
    setkey(chr, start, end)

  # return splitted segs as data.table
  return(splitted.segs)
}




SegmeanToCopyNumber <- function(segmean, cap = 1.5) {
  # limit max and min segmean by cap, can calculate CN from Segmean
  segmean <- replace(segmean, segmean > cap, cap)
  segmean <- replace(segmean, segmean < -cap, -cap)

  # return copy number calculated from segment mean
  return(2^(segmean + 1))
}


GetArmCopyNumber <- function(segs, min.arm.evidence = 100) {
  # get all possible arms
  # arms = c(t(outer(contigs, c("p", "q"), FUN = paste0)))
  # sample.arm.combinations = expand.grid(sample = unique(segs$sample), arm = arms, stringsAsFactors = F)

  # calcualte arm-level segment means
  # if number of probes on an arm is less than min.arm.evidence, set segmean = 0
  arm.copy_number <- segs %>% 
    group_by(sample, arm) %>%
    summarise(arm.nprobes = sum(nprobes), 
              arm.copy_number = weightedMedian(copy_number, nprobes, interpolate = T)) %>% 
    complete(arm, fill = list(arm.nprobes = 0, arm.copy_number = 2)) %>% 
    mutate(dubious_arm = arm.nprobes < min.arm.evidence) %>% 
    arrange(sample, arm) 

  # return 
  return(arm.copy_number)
}


# GetCutoffs return a new cutoffs array that has modified values
# the 1st element of cutoffs: low_cutoff = min(segmean_of_each_arm, default_low_cutoff)
# the 4th element of cutoffs: high_cutoff = max(segmean_of_each_arm, default_high_cutoff) 
GetSampleCopyNumber <- function(arm.copy_number) {
  # set sample copy number low or high as min or max of arm level copy_number
  sample.copy_number <- arm.copy_number %>% 
    filter(!dubious_arm) %>% 
    group_by(sample) %>%
    summarise(low = min(arm.copy_number),
              high = max(arm.copy_number))

  # return sample.copy_number as data.frame  
  return(sample.copy_number)
}




GetSampleThreshold <- function(sample.copy_number, noise) {
  # assign -2/-1/0/1/2 to each gene
  cutoffs <- sample.copy_number %>%
    mutate(t1 = low - noise, 
           t2 = 2 - noise, 
           t3 = 2 + noise,
           t4 = high + noise) %>%
    setDT()
  return(cutoffs)
}




ReadAndFilterGene <- function(file, contigs, header = T) {
  # permissible values
  chr.contigs <- paste0("chr", contigs)

  # read file and set names
  tsv <- fread(file, h = header)
  setnames(tsv, c("chr", "start", "end", "strand", "gene"))

  # reformat Chromosome names, sort and setkeys
  tsv <- tsv %>%  
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



AnnotateGeneByArm <- function(gene.info, arm.info) {

  # annotate gene information by chromosome arm, arm_start and arm_end
  annotated.gene <- foverlaps(arm.info, gene.info, mult="all", type="any", nomatch=0) %>% 
    rename(arm_start = i.start, 
           arm_end = i.end) %>% 
    arrange(chr, start, end) %>% 
    setDT() %>% 
    setkey(chr, start, end) 

  return(annotated.gene)
}





GetGeneSegmentOverlaps <- function(segs, gene.info, arm.copy_number, arm.info, min.seg.evidence = 4, noise = 0.2) {

  # annotated gene.info by arm information
  annotated.gene <- AnnotateGeneByArm(gene.info = gene.info, arm.info = arm.info) 


  # sort segs, define dubious segment, and setkey for later foverlaps 
  segs = segs %>% 
    inner_join(arm.copy_number, ., by = c("sample", "arm"))

    mutate(dubious_seg = nprobes < min.seg.evidence) %>% 
    select(sample, chr, start, end, copy_number, dubious_seg) %>%


    setDT() %>% 
    setkey(chr, start, end)




  # get arm segmean of each gene  
  gene.arm.info = annotated.gene %>% 
    inner_join(., arm.copy_number, by="arm") %>%
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
    mutate(overlap = round((gene_start + gene_end)/2) - end, 
           type = ifelse(is.na(overlap), "arm_to_the_left", "left_neighbor"), 
           copy_number = ifelse(is.na(overlap), median_arm_copy_number, copy_number)) %>%  
    select(sample, gene, copy_number, overlap, is_small_segment, enough_arm_support, median_arm_copy_number, type) %>%
    setDT()

  gene.in_gap.right = gene.in_gap %>% 
    mutate(start = gene_end, 
           end = arm_end) %>% 
    select(sample, chr, start, end, gene, gene_start, gene_end, median_arm_copy_number, enough_arm_support) %>% 
    setDT() %>%
    foverlaps(., segs, type = "any", mult = "first", nomatch = NA) %>% 
    mutate(overlap = start - round((gene_start + gene_end)/2), 
          type = ifelse(is.na(overlap), "arm_to_the_right", "right_neighbor"), 
          copy_number = ifelse(is.na(overlap), median_arm_copy_number, copy_number)) %>%  
    select(sample, gene, copy_number, overlap, is_small_segment, enough_arm_support, median_arm_copy_number, type) %>%
    setDT()


  gene.overlap = data.table(rbind(gene.overlap, gene.in_gap.right, gene.in_gap.left))

  return(gene.overlap)
}





















ReadGeneFile = function(file, header = T) {
  # read file
  tsv = fread(file, h = header)
  setnames(tsv, c("chr", "start", "end", "strand", "gene"))
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









GetSampleSegmean = function(arm.segmean) {
  sample.summary = arm.segmean %>% 
    filter(enough_arm_support) %>% 
    group_by(sample) %>%
    summarise(low = min(median_arm_segmean),
        high = max(median_arm_segmean))
  return(sample.summary)
}






GetThresholedValue = function(sample, copy.number, sample.copy_number.thresholds) {
  # assign -2/-1/0/1/2 to each segmean
  data = data.frame(sample = sample, copy_number = copy.number)
  gene.score = data %>% 
    left_join(sample.copy_number.thresholds, by="sample") %>% 
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










CalcualteCopyNumberEventModel = function(segs) {
  segs = segs[is_small_segment == F] 
  segs$left_segment_start = c(NA, segs$start[-nrow(segs)])
  segs$left_segment_end = c(NA, segs$end[-nrow(segs)])
  dups = duplicated(segs, by=c("sample", "arm"))
  segs = segs[dups==T] %>%
    mutate(distance_start = start - left_segment_start, 
           distance_end = end - left_segment_end)
  samples = unique(segs$sample)
  models = with(segs, tapply(c(rep(1, length(samples)), distance_start, distance_end), 
                             c(samples, rep(sample, 2)), ecdf))
  return(models)
}


# as.integer(gsub("[^0-9]", "", attributes(summary(temp))$header))


AggregateCopyNumber <- function(gene.seg.overlaps, sample.thresholds, noise, event.model) {
  # add copy number scores
  gene.seg.overlaps$score = GetThresholedValue(gene.seg.overlaps$sample, gene.seg.overlaps$copy_number, sample.thresholds)

  # A. Subset overlapped segments
  overlap = gene.seg.overlaps[type=="overlap"]

  # find those that overlapped with more than 1 segment
  dups = duplicated(overlap, by=c("sample", "gene")) | 
         duplicated(overlap, by=c("sample", "gene"), fromLast=T)

  # A.1 Subset genes that overlap with only one segment
  one.overlap = overlap[!dups] %>%
    mutate(type = "single", 
           has_small_segment = is_small_segment) %>%
    select(sample, gene, copy_number, score, type, has_small_segment, enough_arm_support) %>%
    setDT() 

  # A.2 Subset genes that overlap with multiple segments
  mult.overlap = overlap[dups]
  gain = mult.overlap[copy_number >= 2] %>%  # process gain segments
    group_by(sample, gene) %>%
    summarise(aggregate_copy_number = weighted.mean(copy_number, overlap), 
              enough_arm_support = enough_arm_support[1], 
              has_small_segment = sum(is_small_segment) > 0)

  loss = mult.overlap[copy_number <= 2] %>%  # process loss segments
    group_by(sample, gene) %>%
    summarise(aggregate_copy_number = weighted.mean(copy_number, overlap), 
              enough_arm_support = enough_arm_support[1], 
              has_small_segment = sum(is_small_segment) > 0)

  both = rbind(gain, loss) %>%  # merge both segments
    mutate(copy_number = aggregate_copy_number, 
           score = GetThresholedValue(sample, aggregate_copy_number, sample.thresholds), 
           score_factor = factor(score, levels = c(-2, 2, -1, 1, 0)), 
           type = "multiple") %>% 
    arrange(score_factor) %>% 
    select(sample, gene, copy_number, score, type, has_small_segment, enough_arm_support) %>%
    setDT() 

  # detect and remove duplicates (of lower rankings)
  not.preferred = duplicated(both, by=c("sample", "gene"))
  mult.overlap = both[!not.preferred]

  # B. Subset genes that do not overlap with any segments
  no.overlap = gene.cn[type != "overlap"]

  # calculate probability of copy number events based on distance between gene and 
  # segment using sample-specific event models
  no.overlap = no.overlap %>% 
    rowwise() %>% 
    mutate(p = event.model[[sample]](overlap)) %>% 
    setDT()

  # B.1 Subset those genes that have segments in both neibhborhood
  both.side = no.overlap[type == "right_neighbor"] %>%
    inner_join(., no.overlap[type == "left_neighbor"], by=c("sample", "gene")) %>%
    mutate(copy_number = (p.x * (1 - p.y) * copy_number.y + 
                         p.y * (1 - p.x) * copy_number.x + 
                         p.x * p.y * median_arm_copy_number.x) /
                         (p.x + p.y - p.x * p.y), 
           score = GetThresholedValue(sample, copy_number, sample.thresholds), 
           type = "impute", 
           has_small_segment = is_small_segment.x | is_small_segment.y, 
           enough_arm_support = enough_arm_support.x) %>%  
    select(sample, gene, copy_number, score, type, has_small_segment, enough_arm_support) %>% 
    setDT()

  # B.2 Subset those genes that have only one segment in any one neibhborhood
  one.side = no.overlap[type %in% c("right_neighbor", "left_neighbor")] %>%
    inner_join(., no.overlap[type %in% c("arm_to_the_right", "arm_to_the_left")], by=c("sample", "gene")) %>%
    mutate(copy_number = ((1 - p.x) * copy_number.x + p.x * copy_number.y), 
           score = GetThresholedValue(sample, copy_number, sample.thresholds), 
           type = "impute", 
           has_small_segment = is_small_segment.x, 
           enough_arm_support = enough_arm_support.x) %>%  
    select(sample, gene, copy_number, score, type, has_small_segment, enough_arm_support) %>% 
    setDT()

  # B.3 Subset those genes that have no segment in any neibhborhood
  no.side = no.overlap[type == "arm_to_the_right"] %>%
    inner_join(., no.overlap[type == "arm_to_the_left"], by=c("sample", "gene")) %>%
    mutate(copy_number = median_arm_copy_number.x, 
           score = score.x, 
           type = "impute", 
           has_small_segment = F, 
           enough_arm_support = enough_arm_support.x) %>%  
    select(sample, gene, copy_number, score, type, has_small_segment, enough_arm_support) %>% 
    setDT()

  # merge all A and B together
  all = rbind(one.overlap, mult.overlap, both.side, one.side, no.side)

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


