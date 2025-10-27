#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("survival")
#install.packages("survminer")


library(limma)
library(survival)
library(survminer)

expFile="geneExp.txt"     #
cliFile="Survival_SupplementalTable_S1_xena_sp"      #?ٴ??????ļ?

getwd()
setwd("C:\\Users\\zhen-\\Code\\R_code\\R_For_DS_Omics\\05.PFS")      #???ù???Ŀ¼

#
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[3]
data=rt
e=ifelse(data[,gene]>median(data[,gene]), "High", "Low")
data=cbind(as.data.frame(data), Type)

#??ȡ?ٴ??????ļ?
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[,c("PFI.time", "PFI")]
cli=na.omit(cli)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365
cli=as.matrix(cli)
row.names(cli)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(cli))

#???ݺϲ???????????
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(as.data.frame(cli), data)

#?Ƚϸߵͱ?????֮???????????죬?õ???????pֵ
diff=survdiff(Surv(futime, fustat) ~ Type, data=rt)
pValue=1-pchisq(diff$chisq, df=1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=", sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Type, data = rt)

		
#????????????
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           surv.median.line = "hv",
		           legend.title=gene,
		           legend.labs=c("High level", "Low level"),
		           xlab="Time(years)",
		           ylab="Progression free survival",
		           break.time.by = 1,
		           palette=c("red", "blue"),
		           risk.table=F,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)

#????????????
pdf(file=paste0(gene, ".PFS.pdf"), width=5.5, height=5, onefile=FALSE)
print(surPlot)
dev.off()


