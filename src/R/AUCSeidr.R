#' ---
#' title: "seidr gold standard AUC"
#' author: "Camilla Canovi"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Environment
#' Set the working dir
setwd("~/Git/lncRNAs/data/seidr")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="~/Git/lncRNAs/data/seidr")
#' ```
#' Libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(pracma))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
#' Palette
hpal <- colorRampPalette(c("blue","white","red"))(100)
cols <- c(brewer.pal(12,"Paired"), "gray")
mar <- par("mar")
# general 
inputDir <- "roc"
files <- dir(inputDir, pattern="*.roc",full.names=TRUE)
# for each file get plot and the auc
areas <- sapply(files, function(f) {  
  # read the file, skip commented lines  
  dat <- read.delim(f, header=FALSE, skip=1, comment.char="#", col.names=c("TP","FP","PR","ALGORITHM"))  
  
  if("irp"%in%dat$ALGORITHM) {
  dat.1= dat[dat$ALGORITHM=="irp",]
   # read the header
  fileHead <- scan(f, nmax=3, what="character", sep="\t", quiet=T)  
 
  #read the tail
  fileTest <- read_lines(f) %>% subset(grepl("#",.))
  
  vals <- read_lines(f) %>% subset(grepl("#",.))
  lvals <- length(vals)/3
  tabs <- tibble(ALGO=sapply(strsplit(vals[lvals + seq(1,length.out=lvals,by=2)],"\t"),"[",2),
                 #positiveEdges=as.integer(sapply(strsplit(vals[1:lvals],"\t"),"[",2)),
                 #negativeEdges=as.integer(sapply(strsplit(vals[1:lvals],"\t"),"[",3)),
                 AUC=as.double(sub(".* ","",sapply(strsplit(vals[lvals + seq(1,length.out=lvals,by=2)],"\t"),"[",1))),
                 AUPR=as.double(sub(".* ","",sapply(strsplit(vals[lvals + seq(2,length.out=lvals,by=2)],"\t"),"[",1))),
  ) 
  tabs <- tabs[tabs$ALGO=="irp",]
 #print(paste( basename(f),tabs$AUPR))
  # extract the AUC
  auc <- as.numeric(tabs$AUC) 
  # plot the data
  plot(dat.1[,2], dat.1[,1], type="l", 
       main=paste("irp (AUC = ", auc, ")"),
       xlab="False Positive Rate", ylab="True Positive Rate",
       sub=paste(fileHead[2], "Gold Standard edges out of",sum(as.integer(fileHead[2:3])), 
                 "edges for", basename(f))
  )  
  # draw random chance line
  abline(0, 1, lty=2) 
  return(auc)
 }
})
# Name correction, removing preceding folders
names(areas) <- basename(files)
irp <- areas[grepl("percent", names(areas))]
# Result message:
print(paste("The max AUC", irp[which.max(irp)], "is in:", names(which.max(irp))))