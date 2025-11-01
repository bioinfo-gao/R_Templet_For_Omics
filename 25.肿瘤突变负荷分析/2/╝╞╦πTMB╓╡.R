#安装maftools包
#BiocManager::install("maftools")
library(maftools)

#读取合并的maf文件
laml <- read.maf(maf = "combined_maf_value.txt")

#计算tmb值
tmb_table_wt_log = tmb(maf = laml)
write.table(tmb_table_wt_log,file="TMB_log.txt",sep="\t",row.names=F)
