#' ---
#' title: "miRNA precursor - lincRNA blast results"
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
#' A lot of the miRNA precursor are found in the lncRNA set
res <- readBlast(here("data/blastn/linc_network.fasta_Pabies_SE_miRNA.precursor.blt"),
          format=BM8ext)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
