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

#' Only 3647 have a query coverage over 70%
sel <- res$df$query.cum.cov >= 0.7
sum(sel)

#' The ones without species are Pinus taeda, and there are a few from Ginkgo
res$df$Species <- "Pinus_taeda"
sel2 <- grepl("\\|",res$df$subject.id)
res$df$Species[sel2] <- sub(".*\\|","",res$df$subject.id[sel2])
res$df[res$df$Species=="","Species"] <- "Pinus_taeda"

#' How many species and what's their average coverage
res$df$Species <- "Pinus_taeda"
sel2 <- grepl("\\|",res$df$subject.id)
res$df$Species[sel2] <- sub(".*\\|","",res$df$subject.id[sel2])
res$df[res$df$Species=="","Species"] <- "Pinus_taeda"

mar=par("mar")
par(mar=c(10.1,4.1,0.1,0.1))

#' all
boxplot(split(res$df$query.cum.cov,res$df$Species),las=2,ylim=c(0,1))
boxplot(split(res$df$query.cum.cov,res$df$Species),notch=TRUE,las=2,ylim=c(0,0.4))

#' over 70% coverage
boxplot(split(res$df$query.cum.cov[sel],res$df$Species[sel]),las=2,ylim=c(0,1))

#' 707 unique lncRNAs
length(unique(res$df$query.id[sel]))

#' number of hits per lncRNAs (number on the x axis, the y axis is the number of lincRNAs with that many hits)
sp <- split(res$df$Species[sel],res$df$query.id[sel])
par("mar"=mar)
barplot(table(elementNROWS(sp)))

# number of species per lncRNA
barplot(table(elementNROWS(lapply(sp,unique))))
barplot(table(elementNROWS(lapply(sp,unique))),log="y")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
