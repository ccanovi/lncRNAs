

library(readr)
library(dplyr)
tab <- read_tsv("/mnt/picea/storage/reference/Picea-abies/v1.0/annotation/B2G/blast2go_export_20171214_1304.txt.gz")
tab %>% select(`Sequence Name`,`Sequence Description`,`Annotation GO ID`,`Annotation GO Term`,`Annotation GO Category`,`InterPro Accession`,`InterPro Name`)
tab2 <- tab %>% select("Sequence Name","Sequence Description","Annotation GO ID","Annotation GO Term","Annotation GO Category","InterPro Accession","InterPro Name") %>% 
  mutate(`Sequence Name` = gsub("\\.\\d+$", "", `Sequence Name`))
FDN_genes5 <- read_csv("doc/FDNs/genes5.csv",col_names = "Sequence Name")
FDN_genes5_enriched <- semi_join(tab2,FDN_genes5)
write_csv2(FDN_genes5_enriched,"doc/FDNs/genes5_enriched.csv")
nASDFGH <- read_tsv("doc/FDNs/gene_table_FDN5.tsv",col_names = TRUE) %>% 
  rename("Sequence Name"="ID")
nmgtx <- left_join(FDN_genes5_enriched,nASDFGH,by="Sequence Name")
write_tsv(nmgtx,"doc/FDNs/enrichment_genes5.tsv")
yhcf <- read_tsv("doc/FDNs/enrichment_genes11.tsv",col_names = TRUE)

