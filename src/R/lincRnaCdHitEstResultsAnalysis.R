#' ---
#' title: "cd-hit-est clusters"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Libraries
suppressPackageStartupMessages({
  library(here)
  library(IRanges)
  library(readr)
  library(S4Vectors)
})

#' # Data
#' ## 90% percent identity
dat90 <- scan(file=here("data/cdhit_new/linc_network_id90.fasta.clstr"),
            what="character",sep="\n")
pos90 <- grepl(">Cluster",dat90)

rngs90 <- as(Rle(!pos90),"IRanges")

tab90 <- table(width(rngs90))

#' The vast majority is unique
barplot(tab90,log="y")

#' ## 80% percent identity
dat80 <- scan(file=here("data/cdhit_new/linc_network_id80.fasta.clstr"),
              what="character",sep="\n")
pos80 <- grepl(">Cluster",dat80)

rngs80 <- as(Rle(!pos80),"IRanges")

tab80 <- table(width(rngs80))

#' The vast majority is unique
barplot(tab80,log="y")

#' Comparison
n <- as.character(sort(as.integer(union(names(tab90),names(tab80)))))
m <- rbind(id80=tab80[n],id90=tab90[n])
colnames(m) <- n           
barplot(m,beside=TRUE,log="y",
        legend.text=c("id80","id90"),las=2)

#' # Export
dir.create(here("data/analysis/cdhit_new"),recursive=TRUE,showWarnings=FALSE)
write_tsv(data.frame(ID=gsub(".*>|\\.\\.\\..*","",dat90[!pos90]),
           cluster=rep(sub(".* ","",dat90[pos90]),width(rngs90))),
          file=here("data/analysis/cdhit_new/cd-hit-est_id90_cluster-membership.tsv"))

write_tsv(data.frame(ID=gsub(".*>|\\.\\.\\..*","",dat80[!pos80]),
                     cluster=rep(sub(".* ","",dat80[pos80]),width(rngs80))),
          file=here("data/analysis/cdhit_new/cd-hit-est_id80_cluster-membership.tsv"))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
