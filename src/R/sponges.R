#' ---
#' title: "lincRNAs acting as sponges"
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
  (library(here))
  (library(tidyverse))
  (library(RColorBrewer))
  (library(reshape2))
})

#' * Helper files
suppressMessages(source(here("src/R/plotEigengene.R")))

#' # Load dataset and samples
load(here("data/analysis/seidr/network.rda"))
samples <- read_csv(here("doc/samples_B2.csv"))

#' # miRNA, sponge 1 and target genes
plotEigengene(dat, "miRNA_12146-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_12146-5p")
plotEigengene(dat, "TRINITY_DN49527_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN49527_c0_g1_i1")
plotEigengene(dat, "MA_10428924g0020",rep("bla", nrow(dat)),samples$Stages,title="MA_10428924g0020")

#' # miRNA, sponge 2 and target genes
plotEigengene(dat, "miRNA_12146-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_12146-5p")
plotEigengene(dat, "TRINITY_DN10084_c0_g5_i1",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN10084_c0_g5_i1")
plotEigengene(dat, "MA_10428924g0020",rep("bla", nrow(dat)),samples$Stages,title="MA_10428924g0020")

#' # miRNA, sponge 3 and target genes
plotEigengene(dat, "miRNA_12146-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_12146-5p")
plotEigengene(dat, "TRINITY_DN64598_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN64598_c0_g1_i3")
plotEigengene(dat, "MA_10428924g0020",rep("bla", nrow(dat)),samples$Stages,title="MA_10428924g0020")

#' # miRNA, sponge 4 and target genes
plotEigengene(dat, "miRNA_12146-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_12146-5p")
plotEigengene(dat, "TRINITY_DN14989_c0_g3_i1",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN14989_c0_g3_i1")
plotEigengene(dat, "MA_10428924g0020",rep("bla", nrow(dat)),samples$Stages,title="MA_10428924g0020")

#' # miRNA, sponge 5 and target genes
plotEigengene(dat, "miRNA_1613-3p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_1613-3p")
plotEigengene(dat, "TRINITY_DN11832_c0_g1_i5",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN11832_c0_g1_i5")
plotEigengene(dat, "MA_945672g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_945672g0010")
plotEigengene(dat, "MA_488753g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_488753g0010")
plotEigengene(dat, "MA_127386g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_127386g0010")
plotEigengene(dat, "MA_505857g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_505857g0010")
plotEigengene(dat, "MA_10428475g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10428475g0010")

#' # miRNA, sponge 6 and target genes
plotEigengene(dat, "miRNA_1613-3p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_1613-3p")
plotEigengene(dat, "TRINITY_DN11832_c0_g1_i7",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN11832_c0_g1_i7")
plotEigengene(dat, "MA_945672g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_945672g0010")
plotEigengene(dat, "MA_488753g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_488753g0010")
plotEigengene(dat, "MA_127386g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_127386g0010")
plotEigengene(dat, "MA_505857g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_505857g0010")
plotEigengene(dat, "MA_10428475g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10428475g0010")

#' # miRNA, sponge 7 and target genes
plotEigengene(dat, "miRNA_1613-3p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_1613-3p")
plotEigengene(dat, "TRINITY_DN11832_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN11832_c0_g1_i2")
plotEigengene(dat, "MA_945672g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_945672g0010")
plotEigengene(dat, "MA_488753g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_488753g0010")
plotEigengene(dat, "MA_127386g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_127386g0010")
plotEigengene(dat, "MA_505857g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_505857g0010")
plotEigengene(dat, "MA_10428475g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10428475g0010")

#' # miRNA, sponge 8 and target genes
plotEigengene(dat, "miRNA_1613-3p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_1613-3p")
plotEigengene(dat, "TRINITY_DN7435_c0_g1_i7",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN7435_c0_g1_i7")
plotEigengene(dat, "MA_945672g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_945672g0010")
plotEigengene(dat, "MA_488753g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_488753g0010")
plotEigengene(dat, "MA_127386g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_127386g0010")
plotEigengene(dat, "MA_505857g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_505857g0010")
plotEigengene(dat, "MA_10428475g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10428475g0010")

#' # miRNA, sponge 9 and target genes
plotEigengene(dat, "miRNA_203-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_203-5p")
plotEigengene(dat, "TRINITY_DN56799_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN56799_c0_g1_i1")
plotEigengene(dat, "MA_10432357g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10432357g0010")
plotEigengene(dat, "MA_10433485g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10433485g0010")
plotEigengene(dat, "MA_10702g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10702g0010")
plotEigengene(dat, "MA_10426674g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10426674g0010")
plotEigengene(dat, "MA_11005g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_11005g0010")

#' # miRNA, sponge 10 and target genes
plotEigengene(dat, "miRNA_203-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_203-5p")
plotEigengene(dat, "TRINITY_DN64598_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN64598_c0_g1_i3")
plotEigengene(dat, "MA_10432357g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10432357g0010")
plotEigengene(dat, "MA_10433485g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10433485g0010")
plotEigengene(dat, "MA_10702g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10702g0010")
plotEigengene(dat, "MA_10426674g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10426674g0010")
plotEigengene(dat, "MA_11005g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_11005g0010")

#' # miRNA, sponge 11 and target genes
plotEigengene(dat, "miRNA_25425-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_25425-5p")
plotEigengene(dat, "TRINITY_DN22441_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN22441_c0_g1_i1")
plotEigengene(dat, "MA_10437032g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10437032g0010")
plotEigengene(dat, "MA_985g0020",rep("bla", nrow(dat)),samples$Stages,title="MA_985g0020")
plotEigengene(dat, "MA_913482g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_913482g0010")
plotEigengene(dat, "MA_10433629g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10433629g0010")
plotEigengene(dat, "MA_9198405g0020",rep("bla", nrow(dat)),samples$Stages,title="MA_9198405g0020")

#' # miRNA, sponge 12 and target genes
plotEigengene(dat, "miRNA_28675-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_28675-5p")
plotEigengene(dat, "TRINITY_DN22441_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN22441_c0_g1_i1")
plotEigengene(dat, "MA_10437032g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10437032g0010")
plotEigengene(dat, "MA_985g0020",rep("bla", nrow(dat)),samples$Stages,title="MA_985g0020")
plotEigengene(dat, "MA_913482g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_913482g0010")
plotEigengene(dat, "MA_10433629g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10433629g0010")
plotEigengene(dat, "MA_9198405g0020",rep("bla", nrow(dat)),samples$Stages,title="MA_9198405g0020")

#' # miRNA, sponge 13 and target genes
plotEigengene(dat, "miRNA_29263-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_29263-5p")
plotEigengene(dat, "TRINITY_DN22441_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN22441_c0_g1_i1")
plotEigengene(dat, "MA_10437032g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10437032g0010")
plotEigengene(dat, "MA_985g0020",rep("bla", nrow(dat)),samples$Stages,title="MA_985g0020")
plotEigengene(dat, "MA_913482g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_913482g0010")
plotEigengene(dat, "MA_10433629g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10433629g0010")
plotEigengene(dat, "MA_9198405g0020",rep("bla", nrow(dat)),samples$Stages,title="MA_9198405g0020")

#' # miRNA, sponge 14 and target genes
plotEigengene(dat, "miRNA_2964-3p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_2964-3p")
plotEigengene(dat, "TRINITY_DN52354_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN52354_c0_g1_i4")
plotEigengene(dat, "MA_11619g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_11619g0010")
plotEigengene(dat, "MA_10429662g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10429662g0010")

#' # miRNA, sponge 15 and target genes
plotEigengene(dat, "miRNA_31183-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_31183-5p")
plotEigengene(dat, "TRINITY_DN19947_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN19947_c0_g1_i2")
plotEigengene(dat, "MA_205171g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_205171g0010")
plotEigengene(dat, "MA_10434752g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10434752g0010")

#' # miRNA, sponge 16 and target genes
plotEigengene(dat, "miRNA_31183-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_31183-5p")
plotEigengene(dat, "TRINITY_DN11757_c0_g2_i2",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN11757_c0_g2_i2")
plotEigengene(dat, "MA_205171g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_205171g0010")
plotEigengene(dat, "MA_10434752g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10434752g0010")

#' # miRNA, sponge 17 and target genes
plotEigengene(dat, "miRNA_31183-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_31183-5p")
plotEigengene(dat, "TRINITY_DN11757_c0_g2_i3",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN11757_c0_g2_i3")
plotEigengene(dat, "MA_205171g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_205171g0010")
plotEigengene(dat, "MA_10434752g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10434752g0010")

#' # miRNA, sponge 18 and target genes
plotEigengene(dat, "miRNA_32209-3p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_32209-3p")
plotEigengene(dat, "TRINITY_DN9567_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN9567_c0_g1_i4")
plotEigengene(dat, "MA_9061214g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_9061214g0010")
plotEigengene(dat, "MA_10431711g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10431711g0010")
plotEigengene(dat, "MA_820g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_820g0010")
plotEigengene(dat, "MA_39627g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_39627g0010")

#' # miRNA, sponge 19 and no target genes
plotEigengene(dat, "miRNA_3489-3p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_3489-3p")
plotEigengene(dat, "TRINITY_DN20880_c0_g1_i1",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN20880_c0_g1_i1")

#' # miRNA, sponge 20 and target genes
plotEigengene(dat, "miRNA_36092-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_36092-5p")
plotEigengene(dat, "TRINITY_DN13269_c0_g1_i4",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN13269_c0_g1_i4")
plotEigengene(dat, "MA_323184g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_323184g0010")
plotEigengene(dat, "MA_264506g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_264506g0010")
plotEigengene(dat, "MA_13307g0020",rep("bla", nrow(dat)),samples$Stages,title="MA_13307g0020")
plotEigengene(dat, "MA_80318g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_80318g0010")
plotEigengene(dat, "MA_10436921g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10436921g0010")

#' # miRNA, sponge 21 and target genes
plotEigengene(dat, "miRNA_36730-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_36730-5p")
plotEigengene(dat, "TRINITY_DN14121_c0_g1_i2",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN14121_c0_g1_i2")
plotEigengene(dat, "MA_35900g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_35900g0010")

#' # miRNA, sponge 22 and target genes
plotEigengene(dat, "miRNA_5305-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_5305-5p")
plotEigengene(dat, "TRINITY_DN1187_c0_g1_i7",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN1187_c0_g1_i7")
plotEigengene(dat, "MA_10435426g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10435426g0010")
plotEigengene(dat, "MA_10427079g0030",rep("bla", nrow(dat)),samples$Stages,title="MA_10427079g0030")
plotEigengene(dat, "MA_10434973g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10434973g0010")
plotEigengene(dat, "MA_9486006g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_9486006g0010")
plotEigengene(dat, "MA_81888g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_81888g0010")
plotEigengene(dat, "MA_10428192g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10428192g0010")
plotEigengene(dat, "MA_10090939g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10090939g0010")
plotEigengene(dat, "MA_639470g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_639470g0010")
plotEigengene(dat, "MA_251653g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_251653g0010")
plotEigengene(dat, "MA_205171g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_205171g0010")
plotEigengene(dat, "MA_9881339g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_9881339g0010")
plotEigengene(dat, "MA_10430890g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10430890g0010")
plotEigengene(dat, "MA_637235g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_637235g0010")
plotEigengene(dat, "MA_8729g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_8729g0010")
plotEigengene(dat, "MA_10426931g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10426931g0010")

#' # miRNA, sponge 23 and target genes
plotEigengene(dat, "miRNA_5305-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_5305-5p")
plotEigengene(dat, "TRINITY_DN1187_c0_g1_i12",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN1187_c0_g1_i12")
plotEigengene(dat, "MA_10435426g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10435426g0010")
plotEigengene(dat, "MA_10427079g0030",rep("bla", nrow(dat)),samples$Stages,title="MA_10427079g0030")
plotEigengene(dat, "MA_10434973g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10434973g0010")
plotEigengene(dat, "MA_9486006g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_9486006g0010")
plotEigengene(dat, "MA_81888g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_81888g0010")
plotEigengene(dat, "MA_10428192g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10428192g0010")
plotEigengene(dat, "MA_10090939g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10090939g0010")
plotEigengene(dat, "MA_639470g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_639470g0010")
plotEigengene(dat, "MA_251653g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_251653g0010")
plotEigengene(dat, "MA_205171g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_205171g0010")
plotEigengene(dat, "MA_9881339g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_9881339g0010")
plotEigengene(dat, "MA_10430890g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10430890g0010")
plotEigengene(dat, "MA_637235g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_637235g0010")
plotEigengene(dat, "MA_8729g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_8729g0010")
plotEigengene(dat, "MA_10426931g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10426931g0010")

#' # miRNA, sponge 24 and target genes
plotEigengene(dat, "miRNA_5305-5p",rep("bla", nrow(dat)),samples$Stages,title="miRNA_5305-5p")
plotEigengene(dat, "TRINITY_DN12482_c0_g1_i3",rep("bla", nrow(dat)),samples$Stages,title="TRINITY_DN12482_c0_g1_i3")
plotEigengene(dat, "MA_10435426g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10435426g0010")
plotEigengene(dat, "MA_10427079g0030",rep("bla", nrow(dat)),samples$Stages,title="MA_10427079g0030")
plotEigengene(dat, "MA_10434973g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10434973g0010")
plotEigengene(dat, "MA_9486006g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_9486006g0010")
plotEigengene(dat, "MA_81888g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_81888g0010")
plotEigengene(dat, "MA_10428192g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10428192g0010")
plotEigengene(dat, "MA_10090939g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10090939g0010")
plotEigengene(dat, "MA_639470g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_639470g0010")
plotEigengene(dat, "MA_251653g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_251653g0010")
plotEigengene(dat, "MA_205171g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_205171g0010")
plotEigengene(dat, "MA_9881339g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_9881339g0010")
plotEigengene(dat, "MA_10430890g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10430890g0010")
plotEigengene(dat, "MA_637235g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_637235g0010")
plotEigengene(dat, "MA_8729g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_8729g0010")
plotEigengene(dat, "MA_10426931g0010",rep("bla", nrow(dat)),samples$Stages,title="MA_10426931g0010")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
