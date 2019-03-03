library(plyr)
library(dplyr)
library(parallel)

cl <- makeCluster(64, type = "FORK")
setwd("~/biohack.2019/biohack2019/tau/res/")
length.file <- read.table("~/biohack.2019/biohack2019/tau/all_tr.fasta.fai", stringsAsFactors = F)

get.tab <- function(str){
  dd <- data.frame(t(unlist(strsplit(str, "@"))), stringsAsFactors = F)
  names(dd) <- c("id", "chr", "start", "stop", "mm", "gap","pident","length")
  real.length <- as.numeric(length.file[grep(gsub(">","",dd$id), length.file$V1),]$V2)
  dd$start <- as.numeric(dd$start)
  dd$stop <- as.numeric(dd$stop)
  #  dd$length <- dd$stop - dd$start
  #  dd$pid <- abs(100*(dd$length / real.length))
  return(dd)
}
fls <- grep("fasta", list.files("~/biohack.2019/biohack2019/tau/res/"), value = T)

parse.tab <- function(file){
    da <- system(paste("cat ", file, " | grep '>'", sep = ""), intern = T)
    aa <- lapply(da, get.tab)
    aa.df <- do.call(rbind, aa)
    write.table(x = aa.df, file = paste(file, ".tab.parsed", sep = ""),quote = F,sep = "\t", col.names = F, row.names = F)
}



clusterExport(cl = cl, varlist = c("parse.tab", "get.tab", "length.file"))
parLapply(cl = cl,fls, parse.tab)
stopCluster(cl = cl)

##aa.df.90 <- aa.df[aa.df$pid > 70,]


#ret.header <- function(str){
#  return(paste(str$id, str$chr, str$start, str$stop, str$mm, str$gap, sep = "@"))
#}


