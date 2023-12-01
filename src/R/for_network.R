load(here("data/analysis/seidr/network.rda"))
dat <- t(dat)
dat <- as.data.frame(dat)
dat <- rownames_to_column(dat, var = "ID")
lincRNAs <- dat %>% filter(grepl("TRINITY",ID))
load(here("doc/extra_filtering.rda"))
lincRNAs <- lincRNAs[! lincRNAs$ID %in% removing,]
write_tsv(lincRNAs,file=here("doc/lincRNAs_old_network.tsv"))
in_both <- calculate.overlap(list(both=both_final,network=lincRNAs$ID,filename=NULL))
BC_lin <- calculate.overlap(list(BI=only_nc_BC$TRINITY_ID,network=lincRNAs$ID,filename=NULL))
t_map <- calculate.overlap(list(trm=only_nc_trmap$TRINITY_ID,network=lincRNAs$ID,filename=NULL))

BC_1000_lin <- calculate.overlap(list(BI=only_nc_BC_1000$TRINITY_ID,network=lincRNAs$ID,filename=NULL))
lin_1000 <- BC_1000_lin$a2

