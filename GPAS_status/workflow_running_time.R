library(data.table)

setwd("~/Downloads/")

data = fread("GDC-Automation-Metrics-20181026.tsv")

d = data %>% 
  mutate(start = as.integer(as.POSIXct(strptime(start_time, "%a, %d %b %Y %H:%M:%S"))), 
         end = as.integer(as.POSIXct(strptime(end_time, "%a, %d %b %Y %H:%M:%S"))), 
         experimental_strategy = case_when(experimental_strategy == "Targeted Sequencing" ~ "TS", 
                                           experimental_strategy == "RNA-Seq" ~ "RNA", 
                                           TRUE  ~ experimental_strategy), 
         workflow_name = case_when(grepl("align", workflow_name) ~ "Align", 
                              grepl("purecn", workflow_name) ~ "PureCN", 
                              grepl("filtration", workflow_name) ~ "Filter", 
                              workflow_name == "rnaseq-htseq-count" ~ "HTSeq",
                              workflow_name == "somatic-mutation-annotation" ~ "Annotation", 
                              workflow_name == "gatk4-mutect2-tumor-only" ~ "GATK4 Calling", 
                              workflow_name == "somatic-variant-calling.internal_chunk" ~ "Somatic Calling"), 
         workflow = paste(experimental_strategy, workflow_name), 
         workflow = case_when(grepl("Filter", workflow) ~ "Variant Filtering", 
                              TRUE  ~ workflow)) %>% 
  setDT()

# last file vs q90 ratio 
ratio = d %>% 
  group_by(workflow) %>% 
  summarise(median = median(elapsed_seconds),
            q90 = quantile(elapsed_seconds, 0.90), 
            last = max(elapsed_seconds), 
            ratio = last/q90) %>%
  setDT()

# boxplot with changed definition (q10 and q90)
f = function(x) {
  r = quantile(x, probs = c(0.0, 0.10, 0.5, 0.90, 1.00))
  names(r) = c("ymin", "lower", "middle", "upper", "ymax")
  r
}

d$workflow = factor(d$workflow, levels = ratio[order(median)]$workflow)

ggplot(d, aes(x = workflow, y = elapsed_seconds/3600, fill=workflow)) + 
  stat_summary(fun.data = f, geom="boxplot", position="dodge") +
  coord_flip() + 
  ylab("Hours") + 
#  expand_limits(x = 0, y = 0) +
#  scale_x_discrete(expand = c(0.05, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme_bw() + 
  theme(legend.position="none", 
        text = element_text(size=20), 
        axis.title.y = element_blank()) 


# density plot
ggplot(d, aes(x = workflow, y = elapsed_seconds/3600, fill=workflow)) + 
  geom_violin() + 
  coord_flip() + 
  ylab("Hours") + 
  theme_bw() + 
  theme(legend.position = "top", axis.title.y = element_blank()) + 
  scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500)) + 
  theme(legend.position="none") +


# reference time is 10/10/2018 00:00:00
reference = 1539129600
# filter for jobs run after reference time
d.subset = d %>% 
  mutate(start_hr = floor((start - reference) / 3600), 
         end_hr = floor((end - reference) / 3600)) %>% 
  filter(start_hr > 0) %>%
  mutate(range_hr = end_hr - start_hr) %>%
  setDT()
# convert running ranges to multiple individual record
d.expanded = d.subset[rep(1:nrow(d.subset), range_hr + 1)]
d.expanded$hr = unlist(apply(d.subset[, c("start_hr", "end_hr")], 1, function(x) seq(x[1], x[2])))
d.subset2 = d.expanded[project=="CPTAC-3" & !grepl("RNA", workflow)] %>%
  mutate(workflow = factor(workflow, levels = c("WXS Align", "WXS Somatic Calling", "Variant Filtering", "WXS Annotation")))

ggplot(d.subset2, aes(hr, fill = workflow)) + 
  geom_density(alpha = 0.7) + 
  scale_x_continuous(breaks=seq(0, 400, 24), labels = as.character(10:26)) +
  xlab("October Date") + 
  ylab("Density") + 
  theme_bw() + 
  theme(text = element_text(size=20), 
        axis.title.y = element_blank(),
        legend.position = c(0.7, 0.8), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  guides(fill = guide_legend(title = "CPTAC-3 WXS Workflows")) 











