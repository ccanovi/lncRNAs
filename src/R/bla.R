

library(readr)
library(dplyr)
tab <- read_tsv("/mnt/picea/storage/reference/Picea-abies/v1.0/annotation/B2G/blast2go_export_20171214_1304.txt.gz")
tab %>% select(`Sequence Name`,`Sequence Description`,`Annotation GO ID`,`Annotation GO Term`,`Annotation GO Category`,`InterPro Accession`,`InterPro Name`)
tab2 <- tab %>% select("Sequence Name","Sequence Description","Annotation GO ID","Annotation GO Term","Annotation GO Category","InterPro Accession","InterPro Name") %>% 
  mutate(`Sequence Name` = gsub("\\.\\d+$", "", `Sequence Name`))
FDN_genes12 <- read_csv("doc/FDNs/genes12.csv",col_names = "Sequence Name")
FDN_genes12_enriched <- semi_join(tab2,FDN_genes12)
write_csv2(FDN_genes12_enriched,"doc/FDNs/genes12_enriched.csv")
nASDFGH <- read_tsv("doc/FDNs/gene_table_FDN12.tsv",col_names = TRUE) %>% 
  rename("Sequence Name"="ID")
nmgtx <- left_join(FDN_genes12_enriched,nASDFGH,by="Sequence Name")
write_tsv(nmgtx,"doc/FDNs/enrichment_genes12.tsv")
yhcf <- read_tsv("doc/FDNs/enrichment_genes11.tsv",col_names = TRUE)

# First we declare where the file is
statsFile <- here("data/seidr/backbone/statistics.tsv")

# We read the file as a table
stats <- read.table(statsFile, header = T, stringsAsFactors = F, sep='\t')
bla_new <- read_tsv(here("data/seidr/backbone/statistics.tsv"),
                    col_select = c("Source","Target","Type","irp_score",
                                   "D1;D2;D3;D4;D5;D6;D7;D8;D9;D10;D11;NC_Score;NC_SDev"))
 
meh <- bla_new %>% separate("D1;D2;D3;D4;D5;D6;D7;D8;D9;D10;D11;NC_Score;NC_SDev",
                            c("X1","X2","X3","X4","X5","X6","X7",
                              "X8","X9","X10","X11","NC_Score","NC_SDev"),sep=";")
                
bla_final <- meh %>% select(Source,Target,Type,irp_score,NC_Score,NC_SDev)
write.csv(bla_final,file=here("data/analysis/seidr/file_PlantGenie.csv"),row.names = FALSE)
laST_Bv <- read.csv(here("data/analysis/seidr/file_PlantGenie.csv"))
