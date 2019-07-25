# https://jira.opensciencedatacloud.org/rest/api/2/issue/SPT-516?expand=changelog
# save as SPT-516.json


library(rjson)
library(data.table)

date <- '2019-05-10'
date <- as.POSIXct(strptime(date, "%Y-%m-%d"))

# read json
json <- fromJSON(file="SPT-516.json")

# extract histories
histories <- json$changelog$histories

# get created time
created <- sapply(histories, function(x) as.POSIXct(strptime(x$created, "%Y-%m-%dT%H:%M:%S")))
last.index <- whichmax(created < date)
w <- which(created > date)

recent <- GetStatus(histories, max(w))
before <- GetStatus(histories, min(w) - 1)

d <- merge(before, recent, by=c("project", "strategy", "data_type"), all=T, suffixes=c(".old", ""))
d <- d[completed.old != completed | is.na(completed.old), -c("expected.old", "percent.old")]
d$completed.old[which(is.na(d$completed.old))] <- 0

GetStatus <- function(histories, index) {
  h <- histories[[index]]
  # get update string
  ts <- h$items[[1]]$toString

  # format string to table
  tbl <- strsplit(gsub("\\|\\|", "|", gsub("(\\(\\/\\)|\\{\\{|\\}\\}|\\(on\\))", "", ts)), "\n")[[1]]
  tbl <- data.table(t(sapply(tbl, function(x) trimws(strsplit(x, "\\|")[[1]]))))
  tbl <- tbl[, -1]
  names(tbl) <- tolower(unlist(tbl[1]))
  tbl <- tbl[-1]

  # convert character type to numeric
  tbl$completed <- as.integer(tbl$completed)
  tbl$expected <- as.integer(tbl$expected)
  tbl$percent <- as.numeric(tbl$percent)

  # return data.table
  return(tbl)
}




