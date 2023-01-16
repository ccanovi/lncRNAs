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
source(here("UPSCb-common/src/R/featureSelection.R"))
#' # Data genes and TEs
vst_genes_TEs <- read_tsv(here("data/analysis/DE/vst-aware_genes+TEs.tsv"),
                          col_types=cols(
                            .default = col_double(),
                            ID = col_character()
                          )) %>% 
  column_to_rownames("ID")
samples <- read_csv(here("doc/samples_B2.csv"),col_types=cols(.default=col_character())) %>% 
  mutate(ID=sub("_S[0-9]+$","",ID))
sels_genes_TEs <- rangeFeatureSelect(counts=as.matrix(vst_genes_TEs),
                                     conditions=factor(samples$Stages),
                                     nrep=3)
vst.cutoff <- 1
sel_genes_TEs <- vst_genes_TEs[sels_genes_TEs[[vst.cutoff + 1]],]
#' # Data lincRNAs
vst_linc_new <- read_tsv(here("data/analysis/DE/vst-aware_linc_new.tsv"),
                     col_types=cols(
                       .default = col_double(),
                       ID = col_character())) %>% 
  column_to_rownames("ID")
sels_linc_new <- rangeFeatureSelect(counts=as.matrix(vst_linc_new),
                                conditions=factor(samples$Stages),
                                nrep=3)

# select 1.7 as a cutoff
c.sel <- featureSelect(counts=as.matrix(vst_linc_new),
                            conditions=factor(samples$Stages),
                            nrep=3,exp = 1.7)

sel_linc_new <- vst_linc[c.sel,]

#' # Data miRNAs (those are already filtered, just to add)
miRNA <- suppressWarnings(read_csv("sRNA/ShortStack_genome/Pabies_SE_miRNA_filtered_exp1rep3tp1_noOutlier_zinbNormalised_counts.csv",
                                   col_types=cols(
                                     .default = col_double(),
                                     X1 = col_character()))) %>% 
  column_to_rownames("X1")
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
colnames(miRNA) <- colnames(sel_genes_TEs)[match(substr(sample_info_miRNA$ID,8,9),substr(samples$ID,8,9))]
#' combine the mRNA, lncRNA and miRNA
dat <- bind_rows(sel_genes_TEs,sel_linc,miRNA)
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
            file=here("data/analysis/seidr/headless_bla.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)
#' * gene names, one row
write.table(t(sub("\\.1$","",rownames(dat))),
            file=here("data/analysis/seidr/genes_bla.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```