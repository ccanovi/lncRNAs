
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
colnames(counts) <- sub(".*/","",sub("_L00[1,2]_sort.*","",sub("[1,2]_[1].*_.*DXX_","", sub("_dual.*","",salmon))))


# Time specificity
load(here("data/analysis/DE/vst-aware.rda"))
source(here("UPSCb-common/src/R/expressionSpecificityUtility.R"))
samples_m <- read.csv("doc/samples.csv")
time_expression <- expressionSpecificity(exp.mat = vsta[,samples_m$ScilifeID],
                                         tissues = as.character(samples_m$Stages),
                                         output = "complete")

#bla
plot(density(time_expression[!is.na(time_expression[,"score"]),"score"]))
boxplot(time_expression[time_expression[,"score"]==1,"maxn"])
boxplot(time_expression[time_expression[,"score"]==1,"n"])
boxplot(time_expression[time_expression[,"score"]>=0.9,"n"])
time_expression[time_expression[,"score"]>=0.9 & time_expression[,"score"]==8,]
time_expression[time_expression[,"score"]>=0.9 & time_expression[,"n"]==8,]
#means <- rowMeans2(time_expression,cols = TRUE, na.rm = TRUE)
means <- rowMeans2(time_expression,
                   cols = is.vector(c("S1","S2","S3","S4","S5","S6","S7","S8")),
                   na.rm = TRUE)
median <- rowMedians(time_expression)
sd <- rowSds(time_expression)
zero_point <- rowMads(time_expression)


#GMAP

gmap_old <- read_table2("data/GMAP/GMAP_all.gff3",
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

gmap <- read_table2("data/GMAP/GMAP_all.gff3",
                                skip = 1,
                                col_names = c("Gene","X2","X3",
                                              "start","end","X6","strand","X8","path"),
                                col_types = cols_only("Gene" = col_character(),
                                                      "start" = col_double(),
                                                      "end" = col_double(),
                                                      "strand" = col_character(),
                                                      "path" = col_character())) %>% 
  separate(path,
           sep = ";",
           into = c("path_total","Transcript.ID")) %>% 
  mutate(Transcript.ID,Transcript.ID=gsub(".*=","",Transcript.ID)) %>% 
  mutate(path_total,path_total=gsub(".*=","",path_total))

thingy <- read_table2("data/GMAP/GMAP_all.gff3")
  
gmap_tib <- gmap %>% select("Gene","path_total","Transcript.ID")

gmap_mult <- read_table2("data/GMAP/mult_gene",
                         col_names = c("Gene","X2","X3",
                                       "start","end","X6","strand","X8","path"),
                         col_types = cols_only("Gene" = col_character(),
                                               "start" = col_double(),
                                               "end" = col_double(),
                                               "strand" = col_character(),
                                               "path" = col_character())) %>% 
  separate(path,
           sep = ";",
           into = c("path_mult","Transcript.ID")) %>% 
  mutate(Transcript.ID,Transcript.ID=gsub(".*=","",Transcript.ID)) %>% 
  mutate(path_mult,path_mult=gsub(".*=","",path_mult))

gmap_mult_tib <- gmap_mult %>% select("Gene","path_mult","Transcript.ID")

partial <- left_join(gmap_tib, gmap_mult_tib, by = NULL, copy=FALSE)
gmap_uniq <- read_table2("data/GMAP/uniq_gene",
                         col_names = c("Gene","X2","X3",
                                       "start","end","X6","strand","X8","path"),
                         col_types = cols_only("Gene" = col_character(),
                                               "start" = col_double(),
                                               "end" = col_double(),
                                               "strand" = col_character(),
                                               "path" = col_character())) %>% 
  separate(path,
           sep = ";",
           into = c("path_uniq","Transcript.ID")) %>% 
  mutate(Transcript.ID,Transcript.ID=gsub(".*=","",Transcript.ID)) %>% 
  mutate(path_uniq,path_uniq=gsub(".*=","",path_uniq))

gmap_uniq_tib <- gmap_uniq %>% select("Gene","path_uniq","Transcript.ID")

partial2 <- left_join(partial, gmap_uniq_tib, by = NULL, copy=FALSE)

gmap_transloc <- read_table2("data/GMAP/transloc_gene",
                             col_names = c("Gene","X2","X3",
                                           "start","end","X6","strand","X8","path"),
                             col_types = cols_only("Gene" = col_character(),
                                                   "start" = col_double(),
                                                   "end" = col_double(),
                                                   "strand" = col_character(),
                                                   "path" = col_character())) %>% 
  separate(path,
           sep = ";",
           into = c("path_transloc","Transcript.ID")) %>% 
  mutate(Transcript.ID,Transcript.ID=gsub(".*=","",Transcript.ID)) %>% 
  mutate(path_transloc,path_transloc=gsub(".*=","",path_transloc))

gmap_transloc_tib <- gmap_transloc %>% select("Gene","path_transloc","Transcript.ID")

gmap_total <- left_join(partial2, gmap_transloc_tib, by = NULL, copy=FALSE)

#BedToolsIntersect

intersect <- read_table2("data/GMAP/BedToolsIntersect/GMAP_all-Eugene-gene-only.tsv",
                      col_names = c("Gene","X2","X3",
                                    "start","end","X6","strand","X8","path","gene_intersection","X11","X12",
                                    "start_overlap","end_overlap","X15","strand_intersect","X17","GENE_ID","total_intersect"),
                      col_types = cols_only("Gene" = col_character(),
                                            "start" = col_double(),
                                            "end" = col_double(),
                                            "strand" = col_character(),
                                            "path" = col_character(),
                                            "gene_intersection" = col_character(),
                                            "start_overlap" = col_double(),
                                            "end_overlap" = col_double(),
                                            "strand_intersect" = col_character(),
                                            "GENE_ID" = col_character(),
                                            "total_intersect" = col_character())) %>% 
  separate(path,
           sep = ";",
           into = c("path_intersect","Transcript.ID")) %>% 
  mutate(Transcript.ID,Transcript.ID=gsub(".*=","",Transcript.ID)) %>% 
  mutate(path_intersect,path_intersect=gsub(".*=","",path_intersect)) %>% 
  mutate(GENE_ID,GENE_ID=gsub(".*=","",GENE_ID))

bla_intersect <- intersect %>% select("Gene","path_intersect","Transcript.ID","strand_intersect","total_intersect")

intersect_tib <- left_join(gmap_total, bla_intersect, by = NULL, copy=FALSE)


intersect2 <- read_table2("data/GMAP/BedToolsIntersect2/GMAP_all-Eugene-gene-only.tsv",
                          col_names = c("Gene","X2","X3",
                                        "start","end","X6","strand","X8","path","gene_intersection","X11","X12",
                                        "start_overlap","end_overlap","X15","strand_intersect","X17","GENE_ID","total_intersect"),
                          col_types = cols_only("Gene" = col_character(),
                                                "start" = col_double(),
                                                "end" = col_double(),
                                                "strand" = col_character(),
                                                "path" = col_character(),
                                                "gene_intersection" = col_character(),
                                                "start_overlap" = col_double(),
                                                "end_overlap" = col_double(),
                                                "strand_intersect" = col_character(),
                                                "GENE_ID" = col_character(),
                                                "total_intersect" = col_character())) %>% 
  separate(path,
           sep = ";",
           into = c("path_intersect","Transcript.ID")) %>% 
  mutate(Transcript.ID,Transcript.ID=gsub(".*=","",Transcript.ID)) %>% 
  mutate(path_intersect,path_intersect=gsub(".*=","",path_intersect)) %>% 
  mutate(GENE_ID,GENE_ID=gsub(".*=","",GENE_ID))

bla_intersect <- intersect2 %>% select("Gene","path_intersect","Transcript.ID","strand_intersect","total_intersect")

intersect_tib <- left_join(gmap_total, bla_intersect, by = NULL, copy=FALSE)

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
                                              "total_subtract" = col_character())) %>% 
  separate(path,
           sep = ";",
           into = c("path_subtract","Transcript.ID")) %>% 
  mutate(Transcript.ID,Transcript.ID=gsub(".*=","",Transcript.ID)) %>% 
  mutate(path_subtract,path_subtract=gsub(".*=","",path_subtract)) 

bla_sub <- subtract %>% select("Gene","path_subtract","Transcript.ID","strand_subtract","total_subtract")

gmap_final <- left_join(intersect_tib, bla_sub, by = NULL, copy=FALSE)

#trmap Cuffmerge
trmap <- read_table2("data/trmap/CuffMerge/trmap.fasta_filt",
                     col_names = c("type","Gene","strand",
                                   "start","end",
                                   "exon","position"),
                     col_types = cols_only("type"= col_character(),
                                           "Gene" = col_character(),
                                           "strand" = col_character(),
                                           "start" = col_double(),
                                           "end" = col_double(),
                                           "exon" = col_character(),
                                           "position" = col_character()))

trmap_u <- trmap %>% filter(type == "r")

trmap_uniq <- read_delim("data/trmap/trinity/uniq/trmapIDs.txt",
                         delim = "/",
                         col_names = "Transcript.ID")

trmap_mult <- read_delim("data/trmap/trinity/mult/trmapIDs.txt",
                         delim = "/",
                         col_names = "Transcript.ID")
                          
trmap_transloc <- read_delim("data/trmap/trinity/transloc/trmapIDs.txt",
                             delim = "/",
                             col_names = "Transcript.ID")


# Diamond

#diamond <- read_table2("data/DIAMOND/uniref90.dmnd_Trinity.blt.gz",
#                      col_names = c("TRINITY_ID", "reference_ID", "matches", "alignment_length", "mismatch",
#                                    "gap_open", "start_trinity", "end_trinity", "start_ref", "end_ref",
#                                    "evalue", "bitscore", "trinity_length", "ref_length", "taxonomy"))

diamond <- read_table2("data/DIAMOND/uniref90.dmnd_Trinity.blt.gz",
                       col_names = c("TRINITY_ID", "reference_ID", 
                                     "identity", "alignment_length", "mismatch",
                                     "gap_open", "start_trinity", "end_trinity", 
                                     "start_ref", "end_ref",
                                     "evalue", "bitscore", "trinity_length", 
                                     "ref_length", "taxonomy")) %>% 
  mutate(trinity_cov=alignment_length/trinity_length,
         ref_cov=alignment_length/ref_length) %>% 
  group_by(TRINITY_ID) %>% 
  arrange(desc(identity),desc(trinity_cov),desc(ref_cov)) %>% 
  slice(1) %>% ungroup()
# TODO we need to pre-process the results to keep the best hit

uniref_id <- unique(diamond$reference_ID)

# Taxonomy
#tax_map <- readRDS(here("uniref/annotation/uniref90.id.rds"))

f<-function(lines,pos){
    return(lines[lines$V1 %in% uniref_id,])
}

#lines <- read_delim(here("uniref/annotation/uniref90_id-table.txt"),delim=" ",n_max=6)

tax_map <- read_delim_chunked(
    here("uniref/annotation/uniref90_id-table.txt"),
    chunk_size=1e6,delim=" ",callback=DataFrameCallback$new(f))

# Does not look too bad
mar <- par("mar")
par(mar=c(10.1,4.1,0.1,0.1))
barplot(sort(table(tax_map$V2),decreasing=TRUE)[1:20],las=2)
par(mar=mar)


