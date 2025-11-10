#??װmaftools??
#BiocManager::install("maftools")
library(maftools)

#??ȡ?ϲ???maf?ļ?
laml <- read.maf(maf = "combined_maf_value.txt")

?tmb # maf ==> Estimates Tumor Mutation Burden in terms of per megabases
tmb_table_wt_log = tmb(maf = laml)
write.table(tmb_table_wt_log,file="TMB_log.txt",sep="\t",row.names=F)
