#' ---
#' title: "lincRNAs fasta file chunking"
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
  library(Biostrings)
  library(here)
})

#' * Data
fasta <- readDNAStringSet(here("data/analysis/DE/lincRNAs_backbone2.fasta"))

#' # Chunk
chunks <- breakInChunks(length(fasta),nchunk=1000)

#' # Export
dev.null <- sapply(1:length(chunks),
                   function(i,s,e,fa){
                     writeXStringSet(fa[s[i]:e[i]],
                                     file=paste0(here("data/lnctar/"),"linc.",sprintf("%04d",i)))
                   },start(chunks),end(chunks),fasta)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

