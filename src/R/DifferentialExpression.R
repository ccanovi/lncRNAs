#' ---
#' title: "Differential Expression lincRNAs and coding"
#' author: "Camilla Canovi"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup

#' * Libraries
suppressPackageStartupMessages({
(library(data.table))
(library(DESeq2))
(library(gplots))
(library(here))
(library(hyperSpec))
(library(RColorBrewer))
(library(tidyverse))
(library(VennDiagram))
})
#' * Helper files
suppressMessages(source(here("UPSCb-common/src/R/featureSelection.R")))
suppressMessages(source(here("UPSCb-common/src/R/plotMA.R")))
suppressMessages(source(here("UPSCb-common/src/R/volcanoPlot.R")))
suppressMessages(source(here("UPSCb-common/src/R/gopher.R")))

#' * Graphics
pal=brewer.pal(8,"Dark2")
pal12 <- brewer.pal(12,"Paired")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Data lincRNAs
#' ```{r load, echo=FALSE,eval=FALSE}
#' ```
load(here("data/analysis/salmon/dds_linc.rda"))
load(here("data/analysis/DE/vst-aware_linc.rda"))

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Check the different contrasts
resultsNames(dds)

res_2vs1 <- results(dds,c("Stages","S2","S1"), filter = rowMedians(counts(dds)))
res_3vs2 <- results(dds,c("Stages","S3","S2"), filter = rowMedians(counts(dds)))
res_4vs3 <- results(dds,c("Stages","S4","S3"), filter = rowMedians(counts(dds)))
res_5vs4 <- results(dds,c("Stages","S5","S4"), filter = rowMedians(counts(dds)))
res_6vs5 <- results(dds,c("Stages","S6","S5"), filter = rowMedians(counts(dds)))
res_7vs6 <- results(dds,c("Stages","S7","S6"), filter = rowMedians(counts(dds)))
res_8vs7 <- results(dds,c("Stages","S8","S7"), filter = rowMedians(counts(dds)))

#'  # Number of DE genes in different stages of SE
#'  Extract names of all significantly DE genes in the experiment
#' ##all genes
sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
    if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
    if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
    if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
} 

#all
res_list <- list(res_2vs1, res_3vs2, res_4vs3, res_5vs4, res_6vs5, res_7vs6, res_8vs7)
names(res_list) <- c("res_2vs1", "res_3vs2", "res_4vs3", "res_5vs4", "res_6vs5", "res_7vs6", "res_8vs7")
res_sig_list <- lapply(res_list, sigDeg)

all <- lapply(res_sig_list, function(x){
    rownames(x)
})

#up
up <- lapply(res_sig_list, function(x){
    ab <- sigDeg(x, genes="up")
    rownames(ab)
})

#down
down <- lapply(res_sig_list, function(x){
    ab <- sigDeg(x, genes="down")
    rownames(ab)
})

nr_DElincRNAs <- cbind(elementNROWS(up), elementNROWS(down))
colnames(nr_DElincRNAs) <- c("up_regulated", "down_regulated")
rownames(nr_DElincRNAs) <- c("1~2", "2~3", "3~4", "4~5", "5~6", "6~7", "7~8")

barplot(t(nr_DElincRNAs),
        beside = TRUE,
        legend.text = TRUE,
        xlab = "Stage comparison",
        ylab = "Number of DE lincRNAs",
        ylim = c(0,500),
        col = pal12[c(2,3)],
        args.legend = list(bty = "n", x = "top")
)

#' with bigger font size
barplot(t(nr_DElincRNAs), 
        beside = TRUE, 
        col = pal12[c(2,3)], 
        ylim = c(0, 500),
        cex.axis = 1.4, 
        cex.names = 1.4)
mtext(side=1, line=3, "Stage comparison", cex=1.4)
mtext(side=2, line=2.5, "Number of DE lincRNAs", cex=1.4)
mtext(side=3, line=2, "Number of up- and down-regulated lincRNAs", font=2, cex=1.6)

legend("top", bty = "n",
       fill = pal12[c(2,3)],
       legend=c("up-regulated", "down-regulated"), cex = 1.2)

linc <- data.frame(numbers=c(0,247,80,70,65,127,90,
                           0,76,74,59,127,142,132),
                    Stage_comparison = c("1~2", "2~3", "3~4", "4~5", "5~6", "6~7", "7~8",
                                         "1~2", "2~3", "3~4", "4~5", "5~6", "6~7", "7~8"),
                    regulation = c("up-regulated","up-regulated","up-regulated","up-regulated","up-regulated","up-regulated","up-regulated",
                                   "down-regulated","down-regulated","down-regulated","down-regulated","down-regulated","down-regulated","down-regulated"))
l <- ggplot(linc,aes(Stage_comparison,numbers,fill=regulation)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  labs(fill = NULL) +
  labs(x = "Stage comparison") +
  labs(y = "Number of DE lincRNAs") +
  labs(title = "Number of up- and down-regulated lincRNAs") +
  theme_classic() +
  scale_fill_manual(values=pal12[c(2,3)]) +
  theme(text=element_text(size=12)) +
  theme(plot.title = element_text(face = "bold"))
l

#' * Data coding 
#' ```{r load, echo=FALSE,eval=FALSE}
#' ```
load(here("data/analysis/salmon/dds_genes_TEs.rda"))
load(here("data/analysis/DE/vst-aware_genes_TEs.rda"))

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Check the different contrasts
resultsNames(dds)

res_2vs1 <- results(dds,c("Stages","S2","S1"), filter = rowMedians(counts(dds)))
res_3vs2 <- results(dds,c("Stages","S3","S2"), filter = rowMedians(counts(dds)))
res_4vs3 <- results(dds,c("Stages","S4","S3"), filter = rowMedians(counts(dds)))
res_5vs4 <- results(dds,c("Stages","S5","S4"), filter = rowMedians(counts(dds)))
res_6vs5 <- results(dds,c("Stages","S6","S5"), filter = rowMedians(counts(dds)))
res_7vs6 <- results(dds,c("Stages","S7","S6"), filter = rowMedians(counts(dds)))
res_8vs7 <- results(dds,c("Stages","S8","S7"), filter = rowMedians(counts(dds)))

#'  # Number of DE genes in different stages of SE
#'  Extract names of all significantly DE genes in the experiment
#' ##all genes
sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
  if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
  if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
  if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
} 

#all
res_list <- list(res_2vs1, res_3vs2, res_4vs3, res_5vs4, res_6vs5, res_7vs6, res_8vs7)
names(res_list) <- c("res_2vs1", "res_3vs2", "res_4vs3", "res_5vs4", "res_6vs5", "res_7vs6", "res_8vs7")
res_sig_list <- lapply(res_list, sigDeg)

all <- lapply(res_sig_list, function(x){
  rownames(x)
})

#up
up <- lapply(res_sig_list, function(x){
  ab <- sigDeg(x, genes="up")
  rownames(ab)
})

#down
down <- lapply(res_sig_list, function(x){
  ab <- sigDeg(x, genes="down")
  rownames(ab)
})

nr_DEcoding <- cbind(elementNROWS(up), elementNROWS(down))
colnames(nr_DEcoding) <- c("up_regulated", "down_regulated")
rownames(nr_DEcoding) <- c("1~2", "2~3", "3~4", "4~5", "5~6", "6~7", "7~8")

barplot(t(nr_DEcoding),
        beside = TRUE,
        legend.text = TRUE,
        xlab = "Stage comparison",
        ylab = "Number of DE coding",
        ylim = c(0,12000),
        col = pal12[c(2,3)],
        args.legend = list(bty = "n", x = "top")
)

#' with bigger font size
barplot(t(nr_DEcoding), 
        beside = TRUE, 
        col = pal12[c(2,3)], 
        ylim = c(0,12000),
        cex.axis = 1.4, 
        cex.names = 1.4)
mtext(side=1, line=3, "Stage comparison", cex=1.4)
mtext(side=2, line=2.5, "Number of DE coding", cex=1.4)
mtext(side=3, line=2, "Number of up- and down-regulated coding", font=2, cex=1.6)

legend("top", bty = "n",
       fill = pal12[c(2,3)],
       legend=c("up-regulated", "down-regulated"), cex = 1.2)

ga <- data.frame(numbers=c(343,6897,3994,620,5383,11709,3218,
                           541,5282,3015,1939,4901,7848,2751),
                 Stage_comparison = c("1~2", "2~3", "3~4", "4~5", "5~6", "6~7", "7~8",
                                      "1~2", "2~3", "3~4", "4~5", "5~6", "6~7", "7~8"),
                 regulation = c("up-regulated","up-regulated","up-regulated","up-regulated","up-regulated","up-regulated","up-regulated",
                                "down-regulated","down-regulated","down-regulated","down-regulated","down-regulated","down-regulated","down-regulated"))


g <- ggplot(ga,aes(Stage_comparison,numbers,fill=regulation)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  labs(fill = NULL) +
  labs(x = "Stage comparison") +
  labs(y = "Number of DE coding") +
  labs(title = "Number of up- and down-regulated coding") +
  theme_classic() +
  scale_fill_manual(values=pal12[c(2,3)]) +
  theme(text=element_text(size=12)) +
  theme(plot.title = element_text(face = "bold"))
g

#for figure3 after running the 2 BiologicalQA
library(cowplot)
fig2 <- plot_grid(m,g,p,l,labels=c("A","B","C","D"),ncol=2)#vjust = 1.3,ncol=2)#hjust=-20,)
plot(fig2)
ggsave(filename=here("data/analysis/figures//PCA_DE_figure2.png"),device ="png",dpi = 600)

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


