library(ggplot2)
dir <- "~/biohack.2019/biohack2019/tau/res/"
setwd(dir)

ids <- grep("fasta$",grep("Sabrina", list.files(dir), value = T), value = T)
fa.length <- read.table("~/biohack.2019/biohack2019/tau/all_tr.fasta.fai")

get.len <- function(x){
  dd <- as.numeric(unlist(strsplit(x, "@"))[4]) - as.numeric(unlist(strsplit(x, "@"))[3])
  return(100*(dd / len))
}

for(f in ids){
  print(f)
  len <- fa.length[grep(strsplit(f,"-")[[1]][1], fa.length$V1),]$V2
  #paste("cat ", f, " | grep '>'", sep = "")  
  hd <- system(paste("cat ", f, " | grep '>' | sed 's/>//'", sep = ""), intern = T)  
  ss <- abs(unlist(lapply(hd, get.len)))
  ss.par <- ss[ss > 75]
  hd.par <- hd[ss.par]
  system(paste("touch ", f, ".sab.75.NEW.flt", sep = ""))
  for(i in unique(hd.par)){
    print(i)
    sta <- strsplit(i, "@")
    con <- as.numeric(sta[[1]][4]) - as.numeric(sta[[1]][3])
    if(con < 0){
      system(paste("samtools faidx ", f," ",i, " >> ", f, ".sab.NEW.75.flt", sep = ""))
      system(paste("sed -i '1 s|$|-revcomp|'",  " ", f, ".sab.new.75.flt", sep = ""))
    } else {
      system(paste("samtools faidx ", f," ",i , " >> ", f, ".sab.NEW.75.flt", sep = ""))
    }
    
  }
}

setwd(ang)

ang <- "~/biohack.2019/biohack2019/tau/res/angela/75/"
gg <- grep("flt", list.files(ang), value = T)
gg
gg[1]
ang.all <- system(paste("cat ", gg[1], " | grep '>' | sed 's/>//'"), intern = T)

ch <- unlist(lapply(strsplit(ang.all, "@"), function(x)x[[2]]))

table(t(ch))
