#' ---
#' title: "Differential Expression lincRNAs"
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

#' TODO remember to add a function to check for correlation between logFC and 
#' transcript length - check for the effective length

#' * Functions
#' 1. plot specific gene expression
#' ```{r edit1, echo=FALSE,eval=FALSE}
#' CHANGEME - here you need to change the variables in the 
#' plot to display the expression values accross your samples
#' The example below has 2 variables MGenotype and MDay. These 
#' need replacing by the variable(s) of interest in your project
#' ```
"line_plot" <- function(dds=dds,vst=vst,gene_id=gene_id){
    message(paste("Plotting",gene_id))
    sel <- grepl(gene_id,rownames(vsta))
    stopifnot(sum(sel)==2)

    p <- ggplot(bind_cols(as.data.frame(colData(dds)),
                          data.frame(value=vsta[sel,])),
                aes(x=Stages,y=value,group=ID)) +
        geom_point() + geom_smooth() +
        scale_y_continuous(name="VST expression") + 
        ggtitle(label=paste("Expression for: ",gene_id))
    
    suppressMessages(suppressWarnings(plot(p)))
    return(NULL)
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vsta,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("data/analysis/DE"),
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds),
                              expression_cutoff=0,
                              debug=FALSE){
    
    if(length(contrast)==1){
        res <- results(dds,name=contrast)
    } else {
        res <- results(dds,contrast=contrast)
    }
    
    if(plot){
        par(mar=c(5,5,5,5))
        volcanoPlot(res)
        par(mar=mar)
    }
    
    # a look at independent filtering
    if(plot){
        plot(metadata(res)$filterNumRej,
             type="b", ylab="number of rejections",
             xlab="quantiles of filter")
        lines(metadata(res)$lo.fit, col="red")
        abline(v=metadata(res)$filterTheta)
    }
    
    if(verbose){
        message(sprintf("The independent filtering cutoff is %s, removing %s of the data",
                        round(metadata(res)$filterThreshold,digits=5),
                        names(metadata(res)$filterThreshold)))
        
        max.theta <- metadata(res)$filterNumRej[which.max(metadata(res)$filterNumRej$numRej),"theta"]
        message(sprintf("The independent filtering maximises for %s %% of the data, corresponding to a base mean expression of %s (library-size normalised read)",
                        round(max.theta*100,digits=5),
                        round(quantile(counts(dds,normalized=TRUE),probs=max.theta),digits=5)))
    }
    
    if(plot){
        qtl.exp=quantile(counts(dds,normalized=TRUE),probs=metadata(res)$filterNumRej$theta)
        dat <- data.frame(thetas=metadata(res)$filterNumRej$theta,
                          qtl.exp=qtl.exp,
                          number.degs=sapply(lapply(qtl.exp,function(qe){
                              res$padj <= padj & abs(res$log2FoldChange) >= lfc & 
                                  ! is.na(res$padj) & res$baseMean >= qe
                          }),sum))
        
        plot(ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                 geom_line() + geom_point() +
                 scale_x_continuous("quantiles of expression") + 
                 scale_y_continuous("base mean expression") +
                 geom_hline(yintercept=expression_cutoff,
                            linetype="dotted",col="red"))
        
        if(debug){
            p <- ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                geom_line() + geom_point() +
                scale_x_continuous("quantiles of expression") + 
                scale_y_log10("base mean expression") + 
                geom_hline(yintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs)) + 
                     geom_line() + geom_point() +
                     geom_hline(yintercept=dat$number.degs[1],linetype="dashed") +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes"))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs[1] - number.degs),aes()) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Cumulative number of DE genes"))
            
            plot(ggplot(data.frame(x=dat$thetas[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
            
            plot(ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("base mean of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
        }
        p <- ggplot(data.frame(x=dat$qtl.exp[-1],
                               y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
            geom_line() + geom_point() +
            scale_x_log10("base mean of expression") + 
            scale_y_continuous("Number of DE genes per interval") + 
            geom_vline(xintercept=expression_cutoff,
                       linetype="dotted",col="red")
        suppressMessages(suppressWarnings(plot(p)))
    }
    
    sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & 
        res$baseMean >= expression_cutoff
    
    if(verbose){
        message(sprintf("There are %s genes that are DE with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s",
                        sum(sel),
                        lfc,padj,expression_cutoff))
    }
    
    val <- rowSums(vsta[sel,sample_sel])==0
    if (sum(val) >0){
        warning(sprintf("There are %s DE genes that have no vst expression in the selected samples",sum(val)))
        sel[sel][val] <- FALSE
    }    
    
    if(export){
        if(!dir.exists(default_dir)){
            dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
        }
        write.csv(res,file=file.path(default_dir,paste0(default_prefix,"results.csv")))
        write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
    }
    if(plot){
        heatmap.2(t(scale(t(vsta[sel,sample_sel]))),
                  distfun = pearson.dist,
                  hclustfun = function(X){hclust(X,method="ward.D2")},
                  trace="none",col=hpal,labRow = FALSE,
                  labCol=labels[sample_sel]
        )
    }
    return(list(all=rownames(res[sel,]),
                up=rownames(res[sel & res$log2FoldChange > 0,]),
                dn=rownames(res[sel & res$log2FoldChange < 0,])))
}

#' * Data
#' ```{r load, echo=FALSE,eval=FALSE}
#' CHANGEME - here you are meant to load an RData object
#' that contains a DESeqDataSet object. If you ran the 
#' biological QA template, you need not change anything
#' ```
load(here("data/analysis/salmon/dds_linc.rda"))

#' ## Normalisation for visualisation
load(here("data/analysis/DE/vst-aware_linc.rda"))

#' ## Gene of interests
#' ```{r goi, echo=FALSE,eval=FALSE}
#' CHANGEME - Here, you can plot the expression pattern of your gene of
#' interest. You need to have a list of genes in a text file, one geneID per line
#' The ID should exist in your vst data.
#' Note that for the plot to work, you also need to edit the first function (`line_plot`)
#' at the top of this file
#' ```
goi <- read_lines(here("doc/goi.txt"))
stopifnot(all(goi %in% rownames(vsta)))
dev.null <- lapply(goi,line_plot,dds=dds,vst=vsta)

#' ## Differential Expression
#dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Check the different contrasts
resultsNames(dds)

res_2vs1 <- results(dds,c("Stages","S2","S1"), filter = rowMedians(counts(dds)))
#res_3 <- results(dds,c("Stages","S3","res_2vs1"), filter = rowMedians(counts(dds)))
res_3vs2 <- results(dds,c("Stages","S3","S2"), filter = rowMedians(counts(dds)))
res_4vs3 <- results(dds,c("Stages","S4","S3"), filter = rowMedians(counts(dds)))
res_5vs4 <- results(dds,c("Stages","S5","S4"), filter = rowMedians(counts(dds)))
res_6vs5 <- results(dds,c("Stages","S6","S5"), filter = rowMedians(counts(dds)))
res_7vs6 <- results(dds,c("Stages","S7","S6"), filter = rowMedians(counts(dds)))
res_8vs7 <- results(dds,c("Stages","S8","S7"), filter = rowMedians(counts(dds)))

#' for the network

blah <- rownames(res_2vs1)
blahh <- as_tibble(res_2vs1)
blahh <- add_column(blahh,blah)
blahh <- blahh %>% select(log2FoldChange,padj,blah)
colnames(blahh) <- c("log2FoldChange_2vs1","padj_2vs1","blah")
blahh2 <- as_tibble(res_3vs2)
blahh2 <- add_column(blahh2,blah)
blahh2 <- blahh2 %>% select(log2FoldChange,padj,blah)
colnames(blahh2) <- c("log2FoldChange_3vs2","padj_3vs2","blah")
blahh3 <- as_tibble(res_4vs3)
blahh3 <- add_column(blahh3,blah)
blahh3 <- blahh3 %>% select(log2FoldChange,padj,blah)
colnames(blahh3) <- c("log2FoldChange_4vs3","padj_4vs3","blah")
blahh4 <- as_tibble(res_5vs4)
blahh4 <- add_column(blahh4,blah)
blahh4 <- blahh4 %>% select(log2FoldChange,padj,blah)
colnames(blahh4) <- c("log2FoldChange_5vs4","padj_5vs4","blah")
blahh5 <- as_tibble(res_6vs5)
blahh5 <- add_column(blahh5,blah)
blahh5 <- blahh5 %>% select(log2FoldChange,padj,blah)
colnames(blahh5) <- c("log2FoldChange_6vs5","padj_6vs5","blah")
blahh6 <- as_tibble(res_7vs6)
blahh6 <- add_column(blahh6,blah)
blahh6 <- blahh6 %>% select(log2FoldChange,padj,blah)
colnames(blahh6) <- c("log2FoldChange_7vs6","padj_7vs6","blah")
blahh7 <- as_tibble(res_8vs7)
blahh7 <- add_column(blahh7,blah)
blahh7 <- blahh7 %>% select(log2FoldChange,padj,blah)
colnames(blahh7) <- c("log2FoldChange_8vs7","padj_8vs7","blah")
tibble <- left_join(blahh, blahh2, by = NULL, copy=FALSE)
tibble <- left_join(tibble, blahh3, by = NULL, copy=FALSE)
tibble <- left_join(tibble, blahh4, by = NULL, copy=FALSE)
tibble <- left_join(tibble, blahh5, by = NULL, copy=FALSE)
tibble <- left_join(tibble, blahh6, by = NULL, copy=FALSE)
tibble <- left_join(tibble, blahh7, by = NULL, copy=FALSE)
tibble <- column_to_rownames(tibble, var ="blah")
tibble_DE_linc <- tibble
write.table(tibble_DE_linc,sep= "\t",col.names = NA, quote = FALSE, file=here("doc/tibble_DE_linc.tsv"))

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

# export R object with all significantly DEGs
#save(res_sig_list, file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DEGs_DESeqResultsian_p001_lfc05.Rda")

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

nr_DElinc <- cbind(elementNROWS(up), elementNROWS(down))
colnames(nr_DElinc) <- c("up-regulated", "down-regulated")
rownames(nr_DElinc) <- c("1~2", "2~3", "3~4", "4~5", "5~6", "6~7", "7~8")

barplot(t(nr_DElinc),
        beside = TRUE,
        legend.text = TRUE,
        xlab = "Stage comparison",
        ylab = "Number of DE lincRNAs",
        ylim = c(0,2000),
        col = pal12[c(2,3)],
        args.legend = list(bty = "n", x = "top")
)

barplot2(nr_DElinc, 
         beside = TRUE, 
         legend.text = TRUE,
         xlab = "Stage",
         ylab = "Number of DE genes",
         ylim = c(0,2000))

#' with bigger font size
barplot(t(nr_DElinc), 
        beside = TRUE, 
        col = pal12[c(2,3)], 
        ylim = c(0, 2000),
        cex.axis = 1.4, 
        cex.names = 1.4)
mtext(side=1, line=3, "Stage comparison", cex=1.4)
mtext(side=2, line=2.5, "Number of DE linc", cex=1.4)
mtext(side=3, line=2, "Number of up- and down-regulated linc", font=2, cex=1.6)

legend("top", bty = "n",
       fill = pal12[c(2,3)],
       legend=c("up-regulated", "down-regulated"), cex = 1.2)



# working on genes+TEs
load(here("data/analysis/salmon/dds_genes+TEs.rda"))
load(here("data/analysis/DE/vst-aware_genes+TEs.rda"))




#' ## Results
#' #' ```{r res, echo=FALSE,eval=FALSE}
#' CHANGEME - here you need to define the contrast you want to 
#' study - see the example in the next block. 
#' 
#' The `contrast` can be given
#' by name, as a list (numerator/denominator) or as a vector of weight (e.g. c(0,1));
#' read the DESeq2 vignette for more info
#' 
#' The `label` argument is typically one (or a combination) of the metadata stored in colData
#' 
#' The function allows for looking at the independent filtering results using `debug=TRUE`
#' 
#' If you are not satisfied with the default from DESeq2, you can set your own cutoff using `expression_cutoff`
#' 
#' You can change the default output file prefix using `default_prefix`
#' 
#' You can select the set of samples to be added to the `heatmap`, using the `sample_sel` argument. It takes a logical vector.
#' 
#' ```

#' ```{r contrast, echo=FALSE,eval=FALSE}
#'contrast1 <- extract_results(dds=dds,vst=vst,contrast="CHANGEME")
#' ```

#' ### Venn Diagram
#' ```{r venn, echo=FALSE,eval=FALSE}
#' CHANGEME - Here, you typically would have run several contrasts and you want
#' to assess their overlap plotting VennDiagrams.
#' 
#' In the examples below, we assume that these resutls have been saved in a list
#' called `res.list`
#' ```

#' #### All DE genes
#' ```{r venn2, echo=FALSE,eval=FALSE}
#'grid.newpage()
#'grid.draw(venn.diagram(lapply(res.list,"[[","all"),
#'                       NULL,
#'                       fill=pal[1:3]))
#' ```

#' #### DE genes (up in mutant)
#' ```{r venn3, echo=FALSE,eval=FALSE}
#'grid.newpage()
#'grid.draw(venn.diagram(lapply(res.list,"[[","up"),
#'                       NULL,
#'                       fill=pal[1:3]))
#' ```

#' #### DE genes (up in control)
#' ```{r venn4, echo=FALSE,eval=FALSE}
#'grid.newpage()
#'grid.draw(venn.diagram(lapply(res.list,"[[","dn"),
#'                       NULL,
#'                       fill=pal[1:3]))
#' ```

#' ### Gene Ontology enrichment
#' ```{r go, echo=FALSE,eval=FALSE}
#' Once you have obtained a list of candidate genes, you most probably want
#' to annotate them.
#' 
#' In the following example, we first identify the background; _i.e._ the
#' population of expressed genes. We select the genes expressed in a least
#' 2 replicate of one condition at a cutoff of `exp`.
#' 
#' Next we run the enrichment, in the example against `athaliana` using 
#' the gofer3 REST API (interfaced through the gopher.R script loaded at the
#' beginning of this fil).
#' 
#' Finally we export the go enrichment as a complete table and as a table consisting
#' of only the `id` and `padj` columns. The latter can be used as input for _e.g._
#' REVIGO.
#' ```
background <- rownames(vst)[featureSelect(vst,dds$MGenotype,exp=CHANGEME)]

enr.list <- lapply(res.list,function(r){
    lapply(r,gopher,background=background,task="go",url="athaliana")
})

dev.null <- lapply(names(enr.list),function(n){
    r <- enr.list[[n]]
    write_delim(r$all$go,path=file.path(file.path(here("data/analysis/DE",
                                              paste0(n,"-all-DE-genes_GO-enrichment.txt")))))
    write_delim(r$all$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/DE",
                                                       paste0(n,"-all-DE-genes_GO-enrichment_for-REVIGO.txt")))))
    write_csv(r$up$go,path=file.path(file.path(here("data/analysis/DE",
                                                 paste0(n,"-up-DE-genes_GO-enrichment.txt")))))
    write_delim(r$up$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/DE",
                                                                        paste0(n,"-up-DE-genes_GO-enrichment_for-REVIGO.txt")))))    
    write_csv(r$dn$go,path=file.path(file.path(here("data/analysis/DE",
                                                 paste0(n,"-down-DE-genes_GO-enrichment.txt")))))
    write_delim(r$dn$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/DE",
                                                                        paste0(n,"-down-DE-genes_GO-enrichment_for-REVIGO.txt")))))    
})

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


