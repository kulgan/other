
# filter report by epoch period
FilterByPeriod = function(data = report, start_epoch = 0, end_epoch = 1e99) {
  return(data[epoch >= start_epoch & epoch <= end_epoch])
}

# get start epoch of month in format YYYY-MM
GetMonthStart = function(month = month) {
  start_date = paste(month, "01", sep="-")
  return(as.integer(as.POSIXct(strptime(start_date, "%Y-%m-%d"))))
}

# get start epoch of the next month in format YYYY-MM
GetNextMonthStart = function(month = month) {
  start_date = paste(month, "01", sep="-")
  start_date = strptime(start_date, "%Y-%m-%d")
  start_date$mon = start_date$mon + 1
  return(as.integer(as.POSIXct(start_date)))
}


library(data.table)
library(dplyr)

report.file = "GPAS-file-report.20181109.tsv"
month = "2018-10"

report = fread(report.file)
names(report)[which(names(report) == "date_type")] = "data_type"
report$epoch = as.integer(as.POSIXct(strptime(report$created_datetime, "%Y-%m-%dT%H:%M:%S")))
report = FilterByPeriod(report, start_epoch = GetMonthStart(month), end_epoch = GetNextMonthStart(month) - 1)
 
# read rules
rule = fread("file.category.tsv")
rule$index = 1:nrow(rule)

# apply rule
report2  = report %>% 
  full_join(rule[data_type == "*" & strategy == "*", c("node_type", "index")], by = "node_type") %>% 
  full_join(rule[data_type != "*" & strategy == "*", c("node_type", "data_type", "index")], by = c("node_type", "data_type")) %>% 
  full_join(rule[data_type == "*" & strategy != "*", c("node_type", "strategy", "index")], by = c("node_type", "strategy")) %>% 
  full_join(rule[data_type != "*" & strategy != "*", c("node_type", "data_type", "strategy", "index")], by = c("node_type", "data_type", "strategy")) %>%
  filter(!is.na(node_id))

# QC checks
if(nrow(report) != nrow(report2)) {
  cat("Number doesn't match\n")}

index = report2[, which(grepl("^index", colnames(report2)))]
num_index = apply(index, 1, function(x) sum(!is.na(x)))
if(sum(num_index == 1) != nrow(report)) {
  cat("Duplicated rules\n")
}

# apply rule and make summary
rule.index = apply(index, 1, function(x) x[!is.na(x)])
report2$file_category = rule$file_category[rule.index]

summary = report2 %>% 
  group_by(project, strategy, file_category) %>%
  summarise(count = n())

# join export summary and import summary
summary.export = summary %>% 
  filter(file_category == "to_be_exported") %>%
  rename(num_file_to_export = count) %>% 
  select(-file_category) 

summary.io = summary %>% 
  filter(file_category == "imported") %>%
  rename(num_file_imported = count) %>% 
  select(-file_category) %>%
  full_join(summary.export, by = c("project", "strategy")) %>% 
  mutate(num_file_imported = ifelse(is.na(num_file_imported), 0, num_file_imported), 
         num_file_to_export = ifelse(is.na(num_file_to_export), 0, num_file_to_export)) %>% 
  arrange(project, strategy)


# output
out.file = paste0("GPAS_file_import_generation_report_", month, ".tsv")
write.table(summary.io, out.file, col.names=T, row.names=F, sep="\t, quote=F")



