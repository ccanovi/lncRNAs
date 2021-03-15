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
res <- readBlast(here("data/blastn/species_blastn_IS_linc_network.blt"),
          format=BM8ext)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
