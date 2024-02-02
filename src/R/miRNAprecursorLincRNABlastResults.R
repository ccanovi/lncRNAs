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
  library(ggpubr)
  library(tidyverse)
  library(cowplot)
})

#' * Helper files
suppressPackageStartupMessages({
  source(here("UPSCb-common/src/R/blastUtilities.R"))
})

#' # Data
#' Some miRNA precursor are found in the lncRNA set
res_prec <- readBlast(here("precursors/linc_network.fasta_Pabies_SE_miRNA.precursor.blt"),
          format=BM8ext)

#' lincRNA acting as miRNA sponge
#run the function plotEigengene.R in src/R
load(here("data/analysis/seidr/network.rda"))
rownames(dat) <- sub("\\.1","",(rownames(dat)))
dat_tibble <- as_tibble(t(dat))
samples <- read_csv(here("doc/samples_B2.csv"))

sc <- ggplot(dat_tibble, aes(x=`miRNA_7021-5p`, y=TRINITY_DN119111_c0_g1_i1)) + 
  geom_point() +
  geom_smooth(method=lm) +
  theme_classic() +
  stat_cor(method = "pearson",label.y.npc = "top",label.x.npc = "centre")
li <- plotEigengene(t(dat), "TRINITY_DN119111_c0_g1_i1", rep("stages", nrow(t(dat))),samples$Stages,title = "TRINITY_DN119111_c0_g1_i1")
mi <- plotEigengene(t(dat), "miRNA_7021-5p", rep("stages", nrow(t(dat))),samples$Stages,title = "miRNA_7021-5p")
all <- plot_grid(li,mi,sc,labels=c("A","B","C"),vjust = 1.3,hjust=-1.4,nrow=3)#hjust=-20,)
pdf("data/analysis/figures//sponge_correlation.pdf")
plot(all)
dev.off()  

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
