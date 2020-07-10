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
  dat <- read.delim(f, header=FALSE, skip=1, comment.char="#", col.names=c("TP","FP","PR","irp"))  
  # read the header
  fileHead <- scan(f, nmax=3, what="character", sep="\t", quiet=T)  
  #read the tail  
  fileTail <- scan(f, nmax=3, skip=(nrow(dat)+1), what="character", sep="\t", quiet=T)  
  # extract the AUC
  auc <- as.numeric(substring(fileTail[1], 7)) 
  # plot the data
  plot(dat[,2], dat[,1], type="l", 
       main=paste("irp (AUC = ", auc, ")"),
       xlab="False Positive Rate", ylab="True Positive Rate",
       sub=paste(fileHead[2], "Gold Standard edges out of",sum(as.integer(fileHead[2:3])), 
                 "edges for", basename(f))
  )  
  # draw random chance line
  abline(0, 1, lty=2)  
  return(auc)
})
# Name correction, removing preceding folders
names(areas) <- basename(files)
irp <- areas[grepl("percent", names(areas))]
# Result message:
print(paste("The max AUC", irp[which.max(irp)], "is in:", names(which.max(irp))))