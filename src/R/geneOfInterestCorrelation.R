# libraries
library(dplyr)
library(here)
library(readr)

# data
dat <- read_tsv(here("doc/network_matrix.tsv"),show_col_types = FALSE) %>% 
  tibble::column_to_rownames("ID") %>% as.data.frame()

# sample
smpl <- read_csv(here("doc/samples_final.csv"),show_col_types = FALSE) %>% 
  filter(!duplicated(ID)) %>% select(ID,Stages)

# sanity
stopifnot(all(colnames(dat) == smpl$ID))

# calc medians
mat <- sapply(split.data.frame(t(dat),smpl$Stages),matrixStats::colMedians)
rownames(mat) <- rownames(dat)

# goi
goi <- c("miRNA_8938-3p","TRINITY_DN44006_c0_g1_i1", "MA_8886174g0010",
"MA_10430713g0010","MA_709163g0010")#,"TRINITY_DN37193_c0_g1_i1","TRINITY_DN4541_c1_g1_i1",
#"MA_10437010g0020","MA_1537g0010","MA_81888g0010","MA_10433607g0010",
#"MA_17048g0010","MA_42702g0010")#,"MA_168804g0010","MA_10427413g0010","MA_8075895g0010")

# permutations
res <- t(apply(combn(1:length(goi),2),2,function(inx){
  sp <- cor.test(mat[goi,][inx[1],],mat[goi,][inx[2],],method="spearman")
  pe <- cor.test(mat[goi,][inx[1],],mat[goi,][inx[2],],method="pearson")
  
  c(from=goi[inx[1]],
    to=goi[inx[2]],
    spearman=sp$estimate,
    spearman_pval=sp$p.value,
    pearson=pe$estimate,
    pearson_pval=pe$p.value)
}))


