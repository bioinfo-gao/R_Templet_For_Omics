
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("survival")
#install.packages("survminer")


#引用包
library(limma)
library(survival)
library(survminer)

expFile="geneExp.txt"     #表达数据文件
cliFile="Survival_SupplementalTable_S1_xena_sp"      #临床数据文件
#setwd("C:\\biowolf\\Gene\\11.PFS")      #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[3]
data=rt
#根据目标基因表达量对样品进行分组
Type=ifelse(data[,gene]>median(data[,gene]), "High", "Low")
data=cbind(as.data.frame(data), Type)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[,c("PFI.time", "PFI")]
cli=na.omit(cli)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365
cli=as.matrix(cli)
row.names(cli)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(cli))

#数据合并并输出结果
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(as.data.frame(cli), data)

#比较高低表达组之间的生存差异，得到显著性p值
diff=survdiff(Surv(futime, fustat) ~ Type, data=rt)
pValue=1-pchisq(diff$chisq, df=1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=", sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Type, data = rt)

		
#绘制生存曲线
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

#输出生存曲线
pdf(file=paste0(gene, ".PFS.pdf"), width=5.5, height=5, onefile=FALSE)
print(surPlot)
dev.off()


