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
  library(GenomicRanges)
})

#' * Data
metadata <- readRDS(file=here("data/metadata.rds"))

sel <- select(metadata,TRINITY_ID,coding,lincRNAs)
cod <- sel %>% filter(coding==TRUE)
cod <- add_column(cod,Type="coding")
li <- sel %>% filter(lincRNAs==TRUE)
li <- add_column(li,Type="lincRNA")
total <- rbind(cod,li)
total <- select(total,TRINITY_ID,Type)
meta <- left_join(metadata,total, by= "TRINITY_ID")
tot <- meta %>% filter(is.na(Type))
tot <- tot %>% select(TRINITY_ID,Type)
tot <- tot %>%  add_column(type= "other non-coding")
tot <- tot %>%  select(-Type)
total <- rename(total,Type="type")
all <- rbind(total,tot)
meta <- left_join(meta,all,by= "TRINITY_ID")

other_nc <- meta %>% select(TRINITY_ID, non_coding)
other_nc <- other_nc %>% filter(non_coding == TRUE)
 
#' * Directory
dir.create(here("data/analysis/figures"),showWarnings=FALSE)

#' # Stats
#' * How many are lincRNA
sum(metadata$non_coding)

#' * How many are expressed and under regulation (filtered for seidr)
sum(metadata$seidr) 

#' * Do they form gene families?
cl <- metadata %>% select(cdhit_cluster) %>% filter(!is.na(cdhit_cluster))
cl_f <- metadata %>% select(cdhit_cluster,TRINITY_ID,gmap_type,reference_ID,taxon) %>%
                    filter(!is.na(cdhit_cluster))
transloc <- cl_f %>% filter(gmap_type == "transloc")
no_gmap <-  cl_f %>% filter(is.na(gmap_type))

#' How many sequences
nrow(cl)

#' How many unique sequences
nrow(unique(cl))
tab <- table(table(cl))

#' How many clusters with >1 member sequence
sum(tab[-1])

#' How many sequences in such clusters in total
sum(tab[-1]*as.integer(names(tab[-1])))

#' Unique loci multi cluster
mccl <- unique(unlist(cl,use.names = FALSE)[duplicated(unlist(cl))])
df <- metadata %>% filter(cdhit_cluster %in% mccl & !is.na(gmap_loc)) %>% select(gmap_loc,cdhit_cluster) %>% as.data.frame()

#' Clusters split by exact location (definitely an under estimate as some clusters have overlapping loci that could represent isoforms too)
uq <- lapply(split(df$gmap_loc,df$cdhit_cluster),unique)

#' Among the clusters, which have a unique loci
uloci <- uq[elementNROWS(uq)==1]

#' Some of the unique loci are actually multiple / transloc. We ignore that fact for now and simply filter such loci
#' to count only those loci that have multiple isoforms within the exact same coordinates
sum(!grepl("\\|",uloci))

#' * Are any of these conserved in other species
# Check lincRNAPlantGenIEBlastResults.R

#' # Figures
#' ## Expression level
#' Note that the mRNA expression here is probably very biased by the fact that we considered almost 400k transcripts for the vst.
#' Without the limit on the x axis, the mRNA are the ones most highly expressed.
p <- ggplot(meta %>% select(TRINITY_ID,type,starts_with("vst")) %>% 
   filter(type != "other non-coding") %>% 
   mutate(x=rowMeans(across(where(is.numeric))),
         type=type)) +
  aes(x=x,group=type,col=type) +
  geom_density(size=1.5) + ggtitle("average expression distribution") +
  scale_x_continuous(name="averaged vst",limits=c(0,1.5)) +
  theme_classic() +
  theme(text=element_text(size=25)) +
  theme(plot.title = element_text(face = "bold"))

plot(p)

ggsave(filename=here("data/analysis/figures//Figure1-vst-expression-inset_last.png"),device ="png",dpi = 600)
dev.off()
#' The same for the raw data
#' This is more representative of the distribution (we might want to integrate a comparison to the expression data of the "real" transcripts)
p <- ggplot(meta %>% select(TRINITY_ID,type,starts_with("raw")) %>% 
              filter(type != "other non-coding") %>%
              mutate(x=rowMeans(across(where(is.numeric))),
              type=type)) +
       aes(x=x,group=type,col=type) +
  geom_density(size=1.5) + ggtitle("average expression distribution") +
  scale_x_continuous(name="log10(counts)",trans="log10") +
  theme_classic() +
  theme(text=element_text(size=25)) +
  theme(plot.title = element_text(face = "bold"))

plot(p)

ggsave(filename=here("data/analysis/figures//Figure1-raw-expression-inset_last.png"),device ="png",dpi = 600)
dev.off()

#' ## Stage  specificity
p <- ggplot(meta %>% select(TRINITY_ID,type,score) %>% 
              filter(type != "other non-coding") %>%
              mutate(type=type)) +
       aes(x=score,group=type,col=type) +
  geom_density(size=1.5) + ggtitle("stage specificity distribution") +
  scale_x_continuous(name="stage specificity") +
  theme_classic() +
  theme(text=element_text(size=25)) +
  theme(plot.title = element_text(face = "bold"))

plot(p)

ggsave(filename=here("data/analysis/figures//Figure1-stage-specificity-inset_last.png"),device ="png",dpi = 600)
dev.off()

#' ## Number of exons
p <- ggplot(meta %>% select(TRINITY_ID,type,Exon) %>% 
              filter(type != "other non-coding") %>%
              mutate(type=type,
         x=unlist(mclapply(strsplit(Exon,"\\|"),function(e){sort(as.integer(e),decreasing=TRUE)[1]},mc.cores=16L)))) +
  aes(x=x,group=type,col=type) +
  geom_density(size=1.5) + ggtitle("number of exon distribution") +
  scale_x_continuous(name="number of exons",limits=c(1,6)) +
  scale_y_continuous(name="density",limits=c(0,2)) +
  theme_classic() +
  theme(text=element_text(size=25)) +
  theme(plot.title = element_text(face = "bold"))

plot(p)

ggsave(filename=here("data/analysis/figures//Figure1-number-of-exons-inset_last.png"),device ="png",dpi = 600)
dev.off()

#' ## GC content 
p <- ggplot(meta %>% select(TRINITY_ID,type,GC) %>% 
              filter(type != "other non-coding") %>%
              mutate(type=type)) +
       aes(x=GC,group=type,col=type) +
  geom_density(size=1.5) + ggtitle("gc content distribution") +
  scale_x_continuous(name="GC content") +
  theme_classic() +
  theme(text=element_text(size=25)) +
  theme(plot.title = element_text(face = "bold")) 

plot(p)

ggsave(filename=here("data/analysis/figures//Figure1-gc-content-inset_last.png"),device ="png",dpi = 600)
dev.off()

#' boxplot coding - non-coding - score
nc_score <- meta %>% select(score,maxn,coding,lincRNAs) %>% 
  filter(coding|lincRNAs) %>% 
  mutate(type=factor(ifelse(coding,"coding","lincRNA")),
         score=factor(ifelse(score >= 0.9,"specific","aspecific"))) %>% 
  filter(!is.na(score)) %>%
  ggplot(aes(y=maxn,x=type,fill=type)) + 
  geom_boxplot() + 
  theme_classic()  +
  scale_y_continuous(trans="log10",name="expression (log10)") +
  facet_wrap(~ score) +
  theme(text=element_text(size=20)) +
  theme(plot.title = element_text(face = "bold")) +
  theme(legend.position="none")
plot(nc_score)

ggsave(filename=here("data/analysis/figures//nc_score.png"),device ="png",dpi = 600)
dev.off()

# redo, maybe as proportions

perc_spec <- meta %>% select(score,peak,coding,lincRNAs) %>% 
  filter(coding|lincRNAs) %>% 
  mutate(type=factor(ifelse(coding,"coding","lincRNA")),
         score=factor(ifelse(score >= 0.9,"specific","aspecific"))) %>% 
  filter(score=="specific") %>% 
  group_by(score,peak) %>% dplyr::count(type) %>% group_by(score) %>% mutate(n=n/sum(n))  %>% 
  ggplot(aes(x=peak,y=n,fill=type)) + 
  geom_col() +
  theme_classic()  +
  facet_wrap(~ type) +
  theme(text=element_text(size=20)) +
  theme(plot.title = element_text(face = "bold")) +
  theme(legend.position="none") +
  xlab("maximum expression") + ylab("percentage")

plot(perc_spec)

ggsave(filename=here("data/analysis/figures//specificity_percentage.png"),device ="png",dpi = 600)
dev.off()

library(cowplot)
all <- plot_grid(nc_score,perc_spec,labels=c("A","B"),vjust = 1.3) #,ncol=3)#hjust=-20,)
plot(all)
ggsave(filename=here("data/analysis/figures//supplementary_figure.png"),device ="png",dpi = 600)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```


