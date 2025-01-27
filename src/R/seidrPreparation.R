#' ---
#' title: "SE Network with lincRNAs,miRNAs and TEs"
#' author: "Camilla Canovi and Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})
#' Helper
source(here("src/R/featureSelection.R"))
#' # Data genes
samples <- read_csv(here("doc/samples_B2.csv"),col_types=cols(.default=col_character())) %>% 
  mutate(ID=sub("_S[0-9]+$","",ID))

vst_genes <- read_tsv(here("data/analysis/DE/vst-aware_genes_TEs.tsv"),
                      col_types=cols(
                        .default = col_double(),
                        ID = col_character()
                      )) %>% 
  column_to_rownames("ID")
sels_genes <- rangeFeatureSelect(counts=as.matrix(vst_genes),
                                 conditions=factor(samples$Stages),
                                 nrep=3)
vst.cutoff <- 1
sel_genes <- vst_genes[sels_genes[[vst.cutoff + 1]],]
#' # Data lincRNAs
vst_linc <- read_tsv(here("data/analysis/DE/vst-aware_linc.tsv"),
                     col_types=cols(
                       .default = col_double(),
                       ID = col_character())) %>% 
  column_to_rownames("ID")
sels_linc <- featureSelect(counts=as.matrix(vst_linc),
                           conditions=factor(samples$Stages),
                           nrep=3,exp = 0.5)

sel_linc <- vst_linc[sels_linc,]

#' # Data miRNAs (those are already filtered, just to add)
load(here("sRNA/ShortStack_genome/Pabies_SE_miRNA_filtered_exp1rep3tp1_noOutlier_zinbNormalised.rda"))
miRNA <- suppressWarnings(read_csv("sRNA/ShortStack_genome/Pabies_SE_miRNA_filtered_exp1rep3tp1_noOutlier_zinbNormalised_counts.csv",
                                   col_types=cols(
                                     .default = col_double(),
                                     ...1 = col_character()))) %>% 
  column_to_rownames("...1")
sample_info_miRNA <- read_csv("sRNA/ShortStack_genome/Pabies_SE_miRNA_sampleInfo_B2.csv",
                              col_types=cols(.default=col_character())) 

#' filter batch 1
miRNA <- miRNA[,match(sample_info_miRNA$SubmittedID,colnames(miRNA))]
#' reorder as the mRNA
#' The last 2 digits of the ID identify the matched samples
sel <- na.omit(match(substr(samples$ID,8,9),
                     substr(sample_info_miRNA$ID,8,9)))
miRNA <- miRNA[,sel]
sample_info_miRNA <- sample_info_miRNA[sel,]
colnames(miRNA) <- colnames(sel_genes)[match(substr(sample_info_miRNA$ID,8,9),substr(samples$ID,8,9))]
#' combine the mRNA, lncRNA and miRNA
dat <- bind_rows(sel_genes,sel_linc,miRNA)
dat[is.na(dat)] <- 0
save(dat,file=here("data/analysis/seidr/network.rda"))
#write the tsv to do netwok analysis on Cytoscape
#ID <- rownames(dat)
#dat <- cbind(ID,dat)
#dat_tibble <- as_tibble(dat)
#dat_tibble$ID <- sub("\\.1$","",dat_tibble$ID)
#write_tsv(dat_tibble,path=here("doc/network_matrix.tsv"))
#' # Export
dir.create(here("data/analysis/seidr"),showWarnings=FALSE)
#' * gene by column, without names matrix
write.table(t(dat),
            file=here("data/analysis/seidr/headless.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)
#' * gene names, one row
write.table(t(sub("\\.1$","",rownames(dat))),
            file=here("data/analysis/seidr/genes.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```