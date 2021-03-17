#' ---
#' title: "lincRNA -PlantGenIE blast results"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(here)
})

#' * Helper files
suppressPackageStartupMessages({
  source(here("UPSCb-common/src/R/blastUtilities.R"))
})

#' # Data
res <- readBlast(here("data/blastn/plaza_linc_network.blt"),
          format=BM8ext)

#' There are 3911 hits from lincRNA
length(unique(res$df$query.id))

#' to 57286 subjects :-O
length(unique(res$df$subject.id))

#' Only 37 have reciprocal coverage over 70%
sel <- res$df$subject.cum.cov >= 0.7 & res$df$query.cum.cov >= 0.7

#' The ones without species are Pinus taeda, and there are a few from Ginkgo
unique(res$df[sel,"subject.id"])

#' How many species and what's their average coverage
res$df$Species <- "Pinus_taeda"
sel2 <- grepl("\\|",res$df$subject.id)
res$df$Species[sel2] <- sub(".*\\|","",res$df$subject.id[sel2])
res$df[res$df$Species=="","Species"] <- "Pinus_taeda"

mar=par("mar")
par(mar=c(10.1,4.1,0.1,0.1))
boxplot(split(res$df$query.cum.cov,res$df$Species),las=2,ylim=c(0,1))

boxplot(split(res$df$query.cum.cov,res$df$Species),notch=TRUE,las=2,ylim=c(0,0.4))

par("mar"=mar)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
