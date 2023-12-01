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

linc <- read_delim("~/Git/lncRNAs/doc/time_expression_nc_filtered.tsv",
                        delim = " ")
lincRNA <- linc$Transcript.ID 
lincRNA <- as_tibble(lincRNA)
lincRNA <-add_column(lincRNA,linc= TRUE)
metadata <- left_join(metadata, lincRNA,by =c("TRINITY_ID"="value"))
hee <- select(metadata,TRINITY_ID,coding,linc)
co <- hee %>% filter(coding==TRUE)
co <- add_column(co,Type="coding")
li <- hee %>% filter(linc==TRUE)
li <- add_column(li,Type="lincRNA")
total <- rbind(co,li)
total <- select(total,TRINITY_ID,Type)
metadata <- left_join(metadata,total, by= "TRINITY_ID")
ble <- metadata %>% filter(is.na(Type))
nle <- ble %>% select(TRINITY_ID,Type)
nle <- nle %>%  add_column(type= "other non-coding")
nle <- nle %>%  select(-Type)
total <- rename(total,type=Type)
all <- rbind(total,nle)
metadata <- left_join(metadata,all,by= "TRINITY_ID")

TEs <- metadata %>% select(TRINITY_ID,TE_query.id,TE_query.cum.cov,TE_subject.cum.cov,seidr,non_coding,coding,linc)
TEs_filtered <- TEs %>% filter(!is.na(TE_query.id) &
                               !is.na(TE_query.cum.cov) & 
                               !is.na(TE_subject.cum.cov))

only_TEs <- TEs_filtered %>% filter(coding == FALSE &
                                    non_coding == FALSE)

only_FALSE_TEs <- metadata %>% filter(coding == FALSE &
                                  non_coding == FALSE &
                                  linc == TRUE &
                                  seidr == TRUE &
                                  !is.na(TE_query.id) &
                                  !is.na(TE_query.cum.cov) & 
                                  !is.na(TE_subject.cum.cov))

false_linc_network <- metadata %>% filter(non_coding == FALSE &
                                  seidr == TRUE &
                                  linc == TRUE)


other_nc <- metadata %>% select(TRINITY_ID, non_coding)
other_nc <- other_nc %>% filter(non_coding == TRUE)


gmap_bla <- metadata %>% filter(gmap_type == "uniq" |
                                gmap_type == "mult" |
                                gmap_type == "transloc")

gmap_tot <- rbind(gmap_u,gmap_m,gmap_t) 

ti_ex_filt <- metadata %>% filter(!is.na(score)) 
#' * Directory
dir.create(here("data/analysis/figures_new"),showWarnings=FALSE)

#' # Stats
#' * How many are lincRNA
sum(metadata$non_coding)

#' * How many are expressed and under regulation (filtered for seidr)
sum(metadata$seidr) 

#' * Do they form gene families?
cl <- metadata %>% select(cdhit_cluster) %>% filter(!is.na(cdhit_cluster))
bla <- metadata %>% select(cdhit_cluster,TRINITY_ID,gmap_type,reference_ID,taxon) %>%
                    filter(!is.na(cdhit_cluster))
transloc <- bla %>% filter(gmap_type == "transloc")
no_gmap <-  bla %>% filter(is.na(gmap_type))

#' How many sequences
nrow(cl)

#' How many unique sequences
nrow(unique(cl))
tab <- table(table(cl))

#' How many unique sequences in clusters with >1 member sequence
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

linc <- metadata %>% filter(non_coding == TRUE &
                            closest_length > 1000 &
                            seidr == TRUE)
#' * Are any of these conserved in other species
# Check lincRNAPlantGenIEBlastResults.R

#' # Figures
#' ## Expression level
#' Note that the mRNA expression here is probably very biased by the fact that we considered almost 400k transcripts for the vst.
#' Without the limit on the x axis, the mRNA are the ones most highly expressed.
p <- ggplot(metadata %>% select(TRINITY_ID,type,starts_with("vst")) %>% 
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

ggsave(filename=here("data/analysis/figures_new//Figure1-vst-expression-inset_last.png"),device ="png",dpi = 600)
dev.off()
#' The same for the raw data
#' This is more representative of the distribution (we might want to integrate a comparison to the expression data of the "real" transcripts)
p <- ggplot(metadata %>% select(TRINITY_ID,type,starts_with("raw")) %>% 
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

ggsave(filename=here("data/analysis/figures_new//Figure1-raw-expression-inset_last.png"),device ="png",dpi = 600)
dev.off()

#' ## Stage  specificity
p <- ggplot(metadata %>% select(TRINITY_ID,type,score) %>% 
              filter(type != "other non-coding") %>%
              mutate(type=type)) +
       aes(x=score,group=type,col=type) +
  geom_density(size=1.5) + ggtitle("stage specificity distribution") +
  scale_x_continuous(name="stage specificity") +
  theme_classic() +
  theme(text=element_text(size=25)) +
  theme(plot.title = element_text(face = "bold"))

plot(p)

ggsave(filename=here("data/analysis/figures_new//Figure1-stage-specificity-inset_last.png"),device ="png",dpi = 600)
dev.off()

#' ## Number of exons
p <- ggplot(metadata %>% select(TRINITY_ID,type,Exon) %>% 
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

ggsave(filename=here("data/analysis/figures_new//Figure1-number-of-exons-inset_last.png"),device ="png",dpi = 600)
dev.off()

#' ## GC content 
p <- ggplot(metadata %>% select(TRINITY_ID,type,GC) %>% 
              filter(type != "other non-coding") %>%
              mutate(type=type)) +
       aes(x=GC,group=type,col=type) +
  geom_density(size=1.5) + ggtitle("gc content distribution") +
  scale_x_continuous(name="GC content") +
  theme_classic() +
  theme(text=element_text(size=25)) +
  theme(plot.title = element_text(face = "bold")) 

plot(p)

ggsave(filename=here("data/analysis/figures_new//Figure1-gc-content-inset_last.png"),device ="png",dpi = 600)
dev.off()

#' boxplot coding - non-coding - score
nc_score <- metadata %>% select(score,maxn,coding,linc) %>% 
  filter(coding|linc) %>% 
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

ggsave(filename=here("data/analysis/figures_new//nc_score.png"),device ="png",dpi = 600)
dev.off()

  
#+ facet_wrap(linc,type)

# redo, maybe as proportions

perc_spec <- metadata %>% select(score,peak,coding,linc,) %>% 
  filter(coding|linc) %>% 
  mutate(type=factor(ifelse(coding,"coding","lincRNA")),
         score=factor(ifelse(score >= 0.9,"specific","aspecific"))) %>% 
  filter(score=="specific") %>% 
  group_by(score,peak) %>% count(type) %>% group_by(score) %>% mutate(n=n/sum(n))  %>% 
  ggplot(aes(x=peak,y=n,fill=type)) + 
  geom_col() +
  theme_classic()  +
  facet_wrap(~ type) +
  theme(text=element_text(size=20)) +
  theme(plot.title = element_text(face = "bold")) +
  theme(legend.position="none") +
  xlab("maximum expression") + ylab("percentage")

plot(perc_spec)

ggsave(filename=here("data/analysis/figures_new//specificity_percentage.png"),device ="png",dpi = 600)
dev.off()

library(cowplot)
all <- plot_grid(nc_score,perc_spec,labels=c("A","B"),vjust = 1.3) #,ncol=3)#hjust=-20,)
plot(all)
ggsave(filename=here("data/analysis/figures_new//supplementary_figure.png"),device ="png",dpi = 600)


dat <- metadata %>% select(score,peak,coding,linc,) %>%
  filter((coding|linc) & !is.na(score)) %>%
  mutate(type=factor(ifelse(coding,"coding","lincRNA")),
         score=factor(ifelse(score >= 0.9,"specific","aspecific"))) %>%
  group_by(score,peak) %>% count(type) %>% group_by(score) %>% mutate(n=n/sum(n))

ggplot(dat %>% ungroup(),aes(x=peak,y=n,fill=type)) + geom_col() + facet_wrap(vars(score,type))


lalalalala <- metadata %>% select(score,maxn,coding,linc) %>% 
  filter(coding|linc) %>% 
  mutate(type=factor(ifelse(coding,"coding","lincRNA")),
         score=factor(ifelse(score >= 0.9,"specific","aspecific"))) %>% 
  filter(!is.na(score))
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```


