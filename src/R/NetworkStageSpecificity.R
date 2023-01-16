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
  add_column(network$ID)  %>% select(score,S1, S2, S3, S4, S5, S6, S7, S8,peak,maxn,n,"network$ID")

time_expression_network_nc <- time_expression_network[! time_expression_network$`network$ID` %in% removing,]
time_expression_network_nc <- time_expression_network_nc %>% filter(grepl("TRINITY",`network$ID`))

#' # Export
write_tsv(network_tibble,path=here("doc/tissue_specificity_network.tsv"))

# From the old tibble
# Time specificity
load(here("data/analysis/DE/vst-aware.rda"))
source(here("UPSCb-common/src/R/expressionSpecificityUtility.R"))
samples_m <- read.csv("doc/samples.csv")
time_expression <- expressionSpecificity(exp.mat = vsta[,samples_m$ScilifeID],
                                         tissues = as.character(samples_m$Stages),
                                         output = "complete")


time_expression_network_nc <- read_delim(here("doc/time_expression_nc_filtered.tsv"))
time_expression_network_nc <- time_expression_network_nc[! time_expression_network_nc$Transcript.ID %in% removing,]

#bla
plot(density(time_expression_network_nc[!is.na(time_expression_network_nc[,"score"]),"score"]))
boxplot(time_expression_network_nc[time_expression_network_nc[,"score"]==1,"maxn"])
boxplot(time_expression_network_nc[time_expression_network_nc[,"score"]==1,"n"])
boxplot(time_expression_network_nc[time_expression_network_nc[,"score"]>=0.9,"n"])
time_expression_network_nc[time_expression_network_nc[,"score"]>=0.9 & time_expression_network_nc[,"score"]==8,]
time_expression_network_nc[time_expression_network_nc[,"score"]>=0.9 & time_expression_network_nc[,"n"]==8,]
#means <- rowMeans2(time_expression,cols = TRUE, na.rm = TRUE)
means <- rowMeans2(time_expression,
                   cols = is.vector(c("S1","S2","S3","S4","S5","S6","S7","S8")),
                   na.rm = TRUE)
median <- rowMedians(time_expression)
sd <- rowSds(time_expression)
zero_point <- rowMads(time_expression)



#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```


