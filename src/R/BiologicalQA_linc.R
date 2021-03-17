#' ---
#' title: "LincRNAs Biological QA"
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
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(parallel)
  library(pander)
  library(plotly)
  library(tidyverse)
  library(tximport)
  library(vsn)
  library(zinbwave)
  library(Rtsne)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Metadata
#' Sample information

samples <- read_csv(here("doc/samples_final.csv"))

#' # Raw data
filelist <- list.files(here("data/Salmon"), 
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)

#' Sanity check to ensure that the data is sorted according to the sample info

filelist <- filelist[match(samples$ScilifeID,sub("_sortmerna.*","",basename(dirname(filelist))))]

stopifnot(all(match(sub("_sortmerna.*","",basename(dirname(filelist))),
                    samples$ScilifeID) == 1:nrow(samples)))

#' name the file list vector
names(filelist) <- samples$ID

#' Read the expression at the gene level

counts <- suppressMessages(round(tximport(files = filelist, 
                                  type = "salmon",
                                  txOut=TRUE)$counts))

#' combine technical replicates

samples$ID <- sub("_L00[1,2]", "",
                  samples$ScilifeID)
counts <- do.call(
  cbind,
  lapply(split.data.frame(t(counts),
                          samples$ID),
         colSums))

csamples <- samples[,-1]
csamples <- csamples[match(colnames(counts),csamples$ID),]

#' read the expression for the pool of lincRNAs we found 
linc_read <- read_delim("~/Git/lncRNAs/doc/time_expression_nc_filtered.tsv",
                        delim = " ")

linc <- linc_read$Transcript.ID

counts <- counts[linc, ]


#' ## Quality Control
#' * Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' * Let us take a look at the sequencing depth, colouring by Stages

dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(csamples)

ggplot(dat,aes(x,y,fill=csamples$Stages)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90,size=4),axis.title.x=element_blank())

#' * Display the per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected, considering we have lincRNAs, caracterised by a really low signal.
ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual samples colored by Stages. 

dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Stages=csamples$Stages[match(ind,csamples$ID)])

ggplot(dat,aes(x=values,group=ind,col=Stages)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(here("data/analysis/salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("data/analysis/salmon/raw-unormalised-gene-expression_data_linc.csv"))

#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate. 

csamples$Stages <- factor(csamples$Stages)

#' there are a lot of zeros, so we use zinbwave

se <- SummarizedExperiment(assays=list(counts=as.matrix(counts)),
                           colData=as.data.frame(csamples))

zinb <- zinbwave(se,K=0,epsilon=1e12,
                 X="~Stages",
                 observationalWeights=TRUE)

save(zinb,file=here("data/analysis/salmon/zinb.rda"))

#' Check the size factors (_i.e._ the sequencing library size effect)
dds <- DESeqDataSet(zinb,design=~Stages)
dds <- DESeq(dds, 
             sfType = "poscounts", 
             useT = TRUE, 
             minmu = 1e-6)

save(dds,file=here("data/analysis/salmon/dds_linc.rda"))

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)
save(vst,file=here("data/analysis/DE/vst-blind_linc.rda"))

#' ## Variance Stabilising Transformation
vsda <- varianceStabilizingTransformation(dds, blind=FALSE)
vsta <- assay(vsda)
vsta <- vsta - min(vsta)
save(vsta,file=here("data/analysis/DE/vst-aware_linc.rda"))

# prepare the data to build the network
#ID <- rownames(vsta)
#vsta <- cbind(ID,vsta)
#vsta_tibble <- as_tibble(vsta)
#write_tsv(vsta_tibble,path=here("data/analysis/DE/vst-aware_linc.tsv"))
#thing <- read_tsv(here("data/analysis/DE/vst-aware_linc.tsv"))
#library(Biostrings)
#seq <- readDNAStringSet("data/trinity/Trinity.fasta")
#names(seq) <- sub(" .*","",names(seq))
#IDs <- thing$ID
#writeXStringSet(seq[IDs],file="~/Git/lncRNAs/data/analysis/DE/linc_all.fasta")

#' * Validation

#' Check the variance stabilisation. It could be worse, considering we have a pool of lincRNAs.
meanSdPlot(log2(counts(dds)[!mcols(dds)$allZero,]+1))
meanSdPlot(log2(assay(zinb)+1))
meanSdPlot(vst[rowSums(vst)>0,])
meanSdPlot(vsta[rowSums(vsta)>0,]) 

#' ## QC on the normalised data
#' ### PCA
load("data/analysis/DE/vst-aware_linc.rda")
load("data/analysis/salmon/dds_linc.rda")
pc <- prcomp(t(vsta))
percent <- round(summary(pc)$importance[2,]*100)

#' * Cumulative components effect
#' 
#' We define the number of variable of the model
nvar=1

#' And the number of possible combinations

nlevel=nlevels(dds$Stages)

#' We plot the percentage explained by the different components, the
#' red line represent the number of variable in the model, the orange line
#' the number of variable combinations.
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)
  
#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    csamples)

p <- ggplot(pc.dat,aes(x=PC2,y=PC1,col=Stages,text=dds$ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis")
          #,subtitle="variance stabilized counts"

plot(p + labs(x=paste("PC2 (",percent[2],"%)",sep=""),
              y=paste("PC1 (",percent[1],"%)",sep="")))

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")),
         yaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")))

#' ### Heatmap
#' 
#' Filter for noise
#' 
conds <- factor(csamples$Stages)
sels <- rangeFeatureSelect(counts=vsta,
                           conditions=conds,
                           nrep=3)
vst.cutoff <- 1

#' * Heatmap of "all" genes

mar <- par("mar")

par(mar=c(0.05,0.05,0.05,0.05)) 
hm <- heatmap.2(t(scale(t(vsta[sels[[vst.cutoff+1]],]))),
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = conds,
                col=hpal)

plot(as.hclust(hm$colDendrogram),xlab="",sub="",labels=conds)



#' ## Conclusion
# The Biological QA is good, considering it's based on lincRNAs. We have no outliers. 
# The sequencing depth decreased comparing to the previous analysis.
# Looking at the PCA, it could be interesting to do DE analysis to see in particular what's going on
# in S3 and S6. I consider those two stages relevant, because I think things are changes here. 
# 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
