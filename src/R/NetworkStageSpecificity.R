#' ---
#' title: "Network stage specificity"
#' author: "Camilla Canovi & Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

#' Helper functions
source(here("UPSCb-common/src/R/expressionSpecificityUtility.R"))

#' * Data
network <- read_tsv(here("doc/network_matrix.tsv"))

#' * Sample
samples_m <- read.csv("doc/samples_B2.csv")

#' # Stage specificity
time_expression_network <- expressionSpecificity(exp.mat = network[,samples_m$ID],
                                                 tissues = as.character(samples_m$Stages),
                                                 output = "complete")

time_expression_network <- as_tibble(time_expression_network) %>% 
  rename_if(grepl("aij",colnames(.)),funs(str_replace(.,"aij\\.",""))) %>% 
  add_column(network$ID)  %>% select(S1, S2, S3, S4, S5, S6, S7, S8,peak,"network$ID")

#' # Export
write_tsv(network_tibble,path=here("doc/tissue_specificity_network.tsv"))


#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```


