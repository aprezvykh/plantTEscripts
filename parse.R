setwd("~/biohack.2019/fuck/")
strand <- function(row){
  a <- as.numeric(row[["V4"]]) - as.numeric(row[["V3"]])
  if(a<0){
    return("-")
  } else if (a > 0){
    return("+")
  }
}
flip.coords <- function(row){
  start <- row[["V3"]]
  stop <- row[["V4"]]
  if(row[["strand"]] == "-"){
    row[["V3"]] <- stop
    row[["V4"]] <- start
  }  
  return(row)
}
a.bl <- read.table("~/biohack.2019/fuck/D_all.blast", sep = ",")
a.bl$length <- 100*(a.bl$V5 / a.bl$V6)
a.bl$strand <- apply(a.bl,1,strand)


a.p <- a.bl[a.bl$length > 90 & a.bl$V7 > 90,] 
write.table(a.p, "D_all.blast.tsv", sep = "\t", quote = F)

dad <- data.frame()
for(f in 1:nrow(a.p)){
  print(f)
  da <- flip.coords(a.p[f,])
  dad <- rbind(da, dad)
}

a.p <- dad

ids <- paste(a.p$V1,"@",a.p$V2,"@",a.p$V3,
      "@",a.p$V4, "@",a.p$V5, "@",
      a.p$V6, "@",a.p$V7, "@",
      a.p$length, "@", a.p$strand, sep = "")


coords <- paste(a.p$V2, ":", a.p$V3 - 1500, "-", a.p$V4 + 1500, sep = "")

sink("D_multi_fasta.fa")
for(f in 1:length(coords)){
  sa <- system(paste("samtools faidx ../genomes/Aegilops_tauschii.Aet_v4.0.dna.toplevel.fa ", coords[f], sep = ""), intern = T)
  sa <- sub(coords[f], ids[f], sa,fixed = T)
  cat(sa, sep = "\n")
}
sink()
