
library(tidyverse)
library(tximport)
library(here)
dir(here("data/"))

#Salmon

salmon <- list.files("data/Salmon", 
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)

salmon.g <- suppressMessages(tximport(files = salmon, 
                                  type = "salmon",txOut=TRUE))


counts <- round(salmon.g$counts)               
colnames(counts) <- sub(".*/","",sub("_L00[1,2]_sort.*","",salmon))

# Time specificity
source(here("UPSCb-common/src/R/expressionSpecificityUtility.R"))
expressionSpecificity(counts,samples$Stages,output = "complete")

rowMeans2()
rowMedians()
rowSds()
rowMads()


#GMAP

gmap <- read_table2("data/GMAP/GMAP_all.gff3",
                    skip = 1,
                    col_names = c("Gene","X2","X3",
                                  "start","end","X6","strand","X8","path"),
                    col_types = cols_only("Gene" = col_character(),
                                          "start" = col_double(),
                                          "end" = col_double(),
                                          "strand" = col_character(),
                                          "path" = col_character())) %>% 
                        mutate(path,Path=factor(gsub(";.*","",path))) %>% 
                        mutate(path,TRINITY_ID=factor(gsub(".*;|.*=","",path))) 

gmap_mult <- read_table2("data/GMAP/mult_gene",
                         col_names = c("Gene","X2","X3",
                                       "start","end","X6","strand","X8","path"),
                         col_types = cols_only("Gene" = col_character(),
                                               "start" = col_double(),
                                               "end" = col_double(),
                                               "strand" = col_character(),
                                               "path" = col_character()))


gmap_uniq <- read_table2("data/GMAP/uniq_gene",
                         col_names = c("Gene","X2","X3",
                                       "start","end","X6","strand","X8","path"),
                         col_types = cols_only("Gene" = col_character(),
                                               "start" = col_double(),
                                               "end" = col_double(),
                                               "strand" = col_character(),
                                               "path" = col_character()))


gmap_transloc <- read_table2("data/GMAP/transloc_gene",
                             col_names = c("Gene","X2","X3",
                                           "start","end","X6","strand","X8","path"),
                             col_types = cols_only("Gene" = col_character(),
                                                   "start" = col_double(),
                                                   "end" = col_double(),
                                                   "strand" = col_character(),
                                                   "path" = col_character()))


#BedToolsIntersect

intersect <- read_table2("data/GMAP/BedToolsIntersect/GMAP_all-Eugene-gene-only.tsv",
                      col_names = c("Gene","X2","X3",
                                    "start","end","X6","strand","X8","path","gene_intersection","X11","X12",
                                    "start_overlap","end_overlap","X15","strand_overlap","X17","GENE_ID","total_overlap"),
                      col_types = cols_only("Gene" = col_character(),
                                            "start" = col_double(),
                                            "end" = col_double(),
                                            "strand" = col_character(),
                                            "path" = col_character(),
                                            "gene_intersection" = col_character(),
                                            "start_overlap" = col_double(),
                                            "end_overlap" = col_double(),
                                            "strand_overlap" = col_character(),
                                            "GENE_ID" = col_character(),
                                            "total_overlap" = col_character()))

#BedToolsSubtract

subtract <- read_table2("data/GMAP/BedToolsSubtract/GMAP_all-Eugene-gene-only_no-self_no-full-overlap.gff3",
                        col_names = c("Gene","X2","X3",
                                      "start","end","X6","strand","X8","path",
                                      "gene_subtraction","start_subtract","end_subtract",
                                      "X13","X14",
                                      "strand_subtract","total_subtract"),
                        col_types = cols_only("Gene" = col_character(),
                                              "start" = col_double(),
                                              "end" = col_double(),
                                              "strand" = col_character(),
                                              "path" = col_character(),
                                              "gene_subtraction" = col_character(),
                                              "start_subtract" = col_double(),
                                              "end_subtract" = col_double(),
                                              "strand_subtract" = col_character(),
                                              "total_subtract" = col_character()))
# Diamond

diamond <- read_table2("data/DIAMOND/uniref90.dmnd_Trinity.blt.gz",
                      col_names = c("TRINITY_ID", "reference_ID", "matches", "alignment_length", "mismatch",
                                    "gap_open", "start_trinity", "end_trinity", "start_ref", "end_ref",
                                    "evalue", "bitscore", "trinity_length", "ref_length", "taxonomy"))
