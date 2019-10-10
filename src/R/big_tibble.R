library(Biostrings)
library(tidyverse)
library(here)
dir(here("data/"))
ref <- readDNAStringSet("data/trinity/Trinity.fasta")
bla <- sub(" .*","",names(ref))
bla_new <- alphabetFrequency(ref)
bla_GC <- rowSums(alphabetFrequency(ref)[,c("C","G")]) / width(ref)
length_ref <- width(ref)
#' create the big tibble
CNCI_read <- read_delim(file = "data/CNCI/CNCI.index",
                        delim="\t",
                        skip = 1,
                        col_names = c("Transcript.ID", "index", "score", "start", "end", "length")) %>%
    mutate(Transcript.ID,Transcript.ID=gsub(">| .*","",Transcript.ID))
CNCI_tib <- CNCI_read$index

PLEK_read <- read_delim(file="data/PLEK/Trinity.txt",
                        delim="\t",
                        skip = 1,
                        col_names = c("Type","Score","Transcript.ID")) %>% 
    mutate(Transcript.ID,Transcript.ID=gsub(">| .*","",Transcript.ID))
PLEK_index <- left_join(CNCI_read, PLEK_read, by = NULL, copy=FALSE)
PLEK_tib <- PLEK_index$index

CPC2_read <- read_delim(file = "data/CPC2/results.txt",
                        delim="\t")
CPC2_tib <- CPC2_read$label

transdecoder <- read_table2("data/Transdecoder/Trinity.fasta.transdecoder.IDs.txt",
                            col_names = c("Transcript.ID","X2","X3",
                                          "Type","AALength","Score","Coord"),
                            col_types = cols_only("Transcript.ID"=col_character(),
                                                  "Type"=col_character(),
                                                  "AALength"=col_number(),
                                                  "Score"=col_character(),
                                                  "Coord"=col_character())) %>% 
  mutate(Transcript.ID,Transcript.ID=gsub(">| .*","",Transcript.ID)) %>%  
  mutate(Type,Type=factor(sub("type:","",Type))) %>% 
  mutate(Score,Score=parse_double(sub(".*=","",Score),locale=locale(decimal_mark = "."))) %>% 
  mutate(Coord,Strand=factor(gsub(".*\\(|\\)","",Coord))) %>% 
  mutate(Coord,Start=parse_integer(gsub(".*:|-.*","",Coord))) %>% 
  mutate(Coord,End=parse_integer(gsub(".*[0-9]+-|\\(.*","",Coord)))
transdecoder_index <- left_join(CNCI_read, transdecoder, by = NULL, copy=FALSE)
transdecoder_tib <- transdecoder_index$Type

final_tibble <- tibble(TRINITY_ID = bla, GC_content = bla_GC, length = length_ref, CNCI = CNCI_tib, PLEK = PLEK_tib, CPC2 = CPC2_tib)
#Transdecoder = transdecoder_tib)

sbam <- final_tibble %>% filter(length >= 200, CNCI == "noncoding", PLEK == "noncoding", CPC2 == "noncoding")
#Transdecoder == "complete") 
#Transdecoder == "internal", Transdecoder == "3prime_partial", Transdecoder == "5prime_partial")


                  
