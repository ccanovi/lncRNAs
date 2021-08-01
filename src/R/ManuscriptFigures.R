#' ---
#' title: "Manuscript figures"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(here)
  library(matrixStats)
  library(parallel)
  library(tidyverse)
})

#' * Data
metadata <- readRDS(file=here("data/metadata.rds"))

#' * Directory
dir.create(here("data/analysis/figures"),showWarnings=FALSE)

#' # Stats
#' * How many are lincRNA
sum(metadata$non_coding)

#' * How many are expressed and under regulation (filtered for seidr)
sum(metadata$seidr) 

#' * Do they form gene families?
cl <- metadata %>% select(cdhit_cluster) %>% filter(!is.na(cdhit_cluster))
nrow(unique(cl))
table(table(cl))

#' * Are any of these conserved in other species
# Check lincRNAPlantGenIEBlastResults.R

#' # Figures
#' ## Expression level
#' Note that the mRNA expression here is probably very biased by the fact that we considered almost 400k transcripts for the vst.
#' Without the limit on the x axis, the mRNA are the ones most highly expressed.
p <- ggplot(metadata %>% select(TRINITY_ID,seidr,non_coding,starts_with("vst")) %>% 
  mutate(x=rowMeans(across(where(is.numeric))),
         type=ifelse(seidr,"exp. lincRNA",ifelse(non_coding,"lincRNA","mRNA"))),
  aes(x=x,group=type,col=type)) +
  geom_density() + ggtitle("average expression distribution") +
  scale_x_continuous(name="averaged vst",limits=c(0,1.5))

plot(p)

png(filename=here("data/analysis/figures/Figure1-vst-expression-inset.png"))
plot(p)
dev.off()

#' The same for the raw data
#' This is more representative of the distribution (we might want to integrate a comparison to the expression data of the "real" transcripts)
p <- ggplot(metadata %>% select(TRINITY_ID,seidr,non_coding,starts_with("raw")) %>% 
         mutate(x=rowMeans(across(where(is.numeric))),
                type=ifelse(seidr,"exp. lincRNA",ifelse(non_coding,"lincRNA","mRNA"))),
       aes(x=x,group=type,col=type)) +
  geom_density() + ggtitle("average expression distribution") +
  scale_x_continuous(name="log10(counts)",trans="log10")

plot(p)

png(filename=here("data/analysis/figures/Figure1-raw-expression-inset.png"))
plot(p)
dev.off()

#' ## Stage  specificity
p <- ggplot(metadata %>% select(TRINITY_ID,seidr,non_coding,score) %>% 
         mutate(type=ifelse(seidr,"exp. lincRNA",ifelse(non_coding,"lincRNA","mRNA"))),
       aes(x=score,group=type,col=type)) +
  geom_density() + ggtitle("stage specificity distribution") +
  scale_x_continuous(name="stage specificity")

plot(p)

png(filename=here("data/analysis/figures/Figure1-stage-specificity-inset.png"))
plot(p)
dev.off()

#' ## Number of exons
p <- ggplot(metadata %>% select(TRINITY_ID,seidr,non_coding,Exon) %>% 
  mutate(type=ifelse(seidr,"exp. lincRNA",ifelse(non_coding,"lincRNA","mRNA")),
         x=unlist(mclapply(strsplit(Exon,"\\|"),function(e){sort(as.integer(e),decreasing=TRUE)[1]},mc.cores=16L))),
  aes(x=x,group=type,col=type)) +
  geom_density() + ggtitle("number of exon distribution") +
  scale_x_continuous(name="number of exons",limits=c(1,6)) +
  scale_y_continuous(name="density",limits=c(0,2))

plot(p)

png(filename=here("data/analysis/figures/Figure1-number-of-exons-inset.png"))
plot(p)
dev.off()

#' ## GC content 
p <- ggplot(metadata %>% select(TRINITY_ID,seidr,non_coding,GC) %>% 
         mutate(type=ifelse(seidr,"exp. lincRNA",ifelse(non_coding,"lincRNA","mRNA"))),
       aes(x=GC,group=type,col=type)) +
  geom_density() + ggtitle("gc content distribution") +
  scale_x_continuous(name="GC content")

plot(p)

png(filename=here("data/analysis/figures/Figure1-gc-content-inset.png"))
plot(p)
dev.off()

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```


