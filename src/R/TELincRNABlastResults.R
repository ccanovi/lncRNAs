#' ---
#' title: "TE precursor - lincRNA blast results"
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
#' Very few of the lincRNA could be TEs
res <- readBlast(here("data/blastn/linc_network.fasta_TE.blt"),
          format=BM8ext)

#select TEs
sel <- res$df$query.id

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
