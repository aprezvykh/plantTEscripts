library(rtracklayer)
library(foreach)
#setwd("~/biohack.2019/biohack2019/tau/")

transp.name <- "RLC_Tura_Angela_210J24-3"
fa.length <- read.table("~/biohack.2019/biohack2019/tau/all_tr.fasta.fai")
ort <- read.table("~/biohack.2019/biohack2019/prot_prot_Q_aegilops_s_triticum_ortolog.common", header = T, stringsAsFactors = F)
ort$aegilops.real <- unlist(lapply(strsplit(ort$aegilops, "\\."), function(x)x[[1]]))

unique(gtf.ura$seqnames)

parse.genes <- function(dbg,gtf){
  sub.gtf <- gtf[gtf$seqnames == dbg[["chr"]],]
  sub.gtf <- sub.gtf[sub.gtf$type == "transcript",]
  g.idx <- nrow(sub.gtf[sub.gtf$start < dbg[["start"]],])
  prev.chr <- sub.gtf[sub.gtf$start < dbg[["start"]],]$seqnames[g.idx]
  prev.start <- sub.gtf[sub.gtf$start < dbg[["start"]],]$start[g.idx]
  prev.end <- sub.gtf[sub.gtf$start < dbg[["start"]],]$end[g.idx]
  prev.id <- sub.gtf[sub.gtf$start < dbg[["start"]],]$transcript_id[g.idx]
  g2.idx <- sub.gtf[sub.gtf$end > dbg[["start"]],][1,]
  prev.gene <- paste(prev.id, "_", prev.chr, ":", prev.start, "-", prev.end, sep = "")
  next.gene <- paste(g2.idx$transcript_id, "_", g2.idx$seqnames, ":", g2.idx$start, "-", g2.idx$end, sep = "")
  return(data.frame(prev.gene, next.gene))
}

tau <- read.table("~/biohack.2019/biohack2019/tau/res/RLC_Tura_Angela_210J24-3.alignments.fasta.tab.parsed", header = F, stringsAsFactors = F)
names(tau) <- c("id", "chr", "start", "stop", "mm", "gap","pident","length")
tau$abs.length <- abs(tau$stop - tau$start)
tau$trans.length <- fa.length[grep(transp.name, fa.length$V1),]$V2
tau$pident.real <- 100*(tau$abs.length / tau$trans.length)
tau.flt <- tau[tau$pident.real > 90,]

View(tau)

ura <- read.table("~/biohack.2019/biohack2019/ura/res/RLC_Tura_Angela_210J24-3.alignments.fasta.tab.parsed", header = F, stringsAsFactors = F)
names(ura) <- c("id", "chr", "start", "stop", "mm", "gap","pident","length")
ura$abs.length <- abs(ura$stop - ura$start)
ura$trans.length <- fa.length[grep(transp.name, fa.length$V1),]$V2
ura$pident.real <- 100*(ura$abs.length / ura$trans.length)
ura.flt <- ura[ura$pident.real > 90,]


gtf.tau <- data.frame(rtracklayer::import("~/biohack.2019/gtf/Aegilops_tauschii.Aet_v4.0.42.gtf"))
gtf.ura <- data.frame(rtracklayer::import("~/biohack.2019/gtf/Triticum_urartu.ASM34745v1.42.gtf"))

big.flt.ura <- data.frame()
for(f in 1:nrow(ura.flt)){
  print(f)
  sas <- ura.flt[f,]  
  crd <- parse.genes(sas, gtf.ura)
  big <- bind_cols(sas, crd)
  prev.parsed <- unlist(strsplit(unlist(strsplit(as.character(big$prev.gene), "_"))[1],"\\."))[1]
  next.parsed <- unlist(strsplit(unlist(strsplit(as.character(big$next.gene), "_"))[1],"\\."))[1]
  prev.par <- ort[grep(prev.parsed, ort$aegilops.real),]$triticum
  next.par <- ort[grep(next.parsed, ort$aegilops.real),]$triticum
  prev.par
  if(identical(prev.par, character(0))){
    big$prev.par <- c("No paralog")
  } else {
    big$prev.par <- prev.par
  }
  if(identical(next.par, character(0))){
    big$next.par <- c("No paralog")
  } else {
    big$next.par <- next.par
  }
  big.flt.tau <- rbind(big, big.flt.tau)
}


#write.table(x = big.flt, file =  paste(transp.name, "ind_with_par.table",sep = ""),sep = "\t")

big.flt.ura <- data.frame()
for(f in 1:nrow(ura.flt)){
  print(f)
  sas <- ura.flt[f,]  
  crd <- parse.genes(sas, gtf.ura)
  big <- bind_cols(sas, crd)
  prev.parsed <- unlist(strsplit(unlist(strsplit(as.character(big$prev.gene), "_"))[1],"\\."))[1]
  next.parsed <- unlist(strsplit(unlist(strsplit(as.character(big$next.gene), "_"))[1],"\\."))[1]
}
