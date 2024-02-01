library(here)
library(tidyverse)
library(tximport)

#getting the TPM values for PlantGenie for lincRNAs

samples <- read_csv(here("doc/samples_final.csv"))
load(here("data/analysis/salmon/dds_linc.rda"))

files <- list.files(here("data/Salmon"), 
                    recursive = TRUE, 
                    pattern = "quant.sf",
                    full.names = TRUE)
# filter the files
files <- files[match(samples$ScilifeID,sub("_sortm.*","",basename(dirname(files))))]

txi <- tximport(files = files, 
               type = "salmon",
               txOut=TRUE)

tpm <-  do.call(
  cbind,
  lapply(split.data.frame(t(txi$abundance),
                          samples$ID),
         colMeans))

colnames(tpm) <- c("Stage1","Stage1","Stage1","Stage2","Stage2","Stage2","Stage3","Stage3","Stage3","Stage4","Stage4","Stage4",
                         "Stage5","Stage5","Stage5","Stage6","Stage6","Stage6","Stage7","Stage7","Stage7","Stage8","Stage8","Stage8")
colnames(dds) <- c("Stage1","Stage1","Stage1","Stage2","Stage2","Stage2","Stage3","Stage3","Stage3","Stage4","Stage4","Stage4",
  "Stage5","Stage5","Stage5","Stage6","Stage6","Stage6","Stage7","Stage7","Stage7","Stage8","Stage8","Stage8")
rownames(tpm) <- read_tsv(files[[length(files)]],show_col_types = FALSE) %>%
  select(Name) %>% unlist(use.names = FALSE)

linc_tpm <- tpm[rownames(tpm) %in% rownames(dds),]

write.csv(linc_tpm,file=here("data/analysis/salmon/tpm_linc_values.csv"),row.names = TRUE)

#getting the TPM values for PlantGenie for genes

files <- list.files(here("data/salmon"), 
                    recursive = TRUE, 
                    pattern = "quant.sf",
                    full.names = TRUE)
# filter the files
files <- files[match(samples$ScilifeID,sub("_sortm.*","",basename(dirname(files))))]

txi <- tximport(files = files, 
                type = "salmon",
                txOut=TRUE)

tpm <-  do.call(
  cbind,
  lapply(split.data.frame(t(txi$abundance),
                          samples$ID),
         colMeans))

colnames(tpm) <- c("Stage1","Stage1","Stage1","Stage2","Stage2","Stage2","Stage3","Stage3","Stage3","Stage4","Stage4","Stage4",
                   "Stage5","Stage5","Stage5","Stage6","Stage6","Stage6","Stage7","Stage7","Stage7","Stage8","Stage8","Stage8")
#colnames(dds) <- c("Stage1","Stage1","Stage1","Stage2","Stage2","Stage2","Stage3","Stage3","Stage3","Stage4","Stage4","Stage4",
#                  "Stage5","Stage5","Stage5","Stage6","Stage6","Stage6","Stage7","Stage7","Stage7","Stage8","Stage8","Stage8")
rownames(tpm) <- read_tsv(files[[length(files)]],show_col_types = FALSE) %>%
  select(Name) %>% unlist(use.names = FALSE)
rownames(tpm) <- sub("\\.[0-9]+","",rownames(tpm))
genes_tpm <- tpm[rownames(tpm) %in% ge,]

write.csv(genes_tpm,file=here("data/analysis/salmon/tpm_genes_values.csv"),row.names = TRUE)

# linc for gbrowse
linc_gbrowse <- read_tsv(here("doc/lincRNAs.tsv"))
for_gbrowse <- linc_gbrowse[linc_gbrowse$TRINITY_ID %in% li,]
for_gbrowse <- for_gbrowse %>% select(TRINITY_ID,scaffold,t_start,t_end,strand)
write.csv(for_gbrowse,file=here("doc/PlantGenIE_gbrowse.csv"),row.names = TRUE)

# We read the file as a table

network <- read_tsv(here("data/seidr/backbone/statistics.tsv"))
 
net_new <- network %>% separate("D1;D2;D3;D4;D5;D6;D7;D8;D9;NC_Score;NC_SDev",
                            c("X1","X2","X3","X4","X5","X6","X7",
                              "X8","X9","NC_Score","NC_SDev"),sep=";")
network_final <- net_new %>% select(Source,Target,Type,irp_score,NC_Score,NC_SDev)
write.csv(network_final,file=here("data/analysis/seidr/file_PlantGenie_linc.csv"),row.names = FALSE)

