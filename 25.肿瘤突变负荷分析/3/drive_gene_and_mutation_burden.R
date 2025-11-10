setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/25.肿瘤突变负荷分析/3")


#???ð?
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

expFile="geneExp.txt"     #?????????ļ?
tmbFile="TMB.txt"         #????ͻ???????ļ?

#??ȡ?????????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[1]

#ɾ????????Ʒ
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\", rownames(tumorData))
data=avereps(tumorData)

#??ȡ????ͻ?为?ɵ??ļ?
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(tmb)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\", rownames(tmb))
tmb=avereps(tmb)
#???ݺϲ???????????
sameSample=intersect(row.names(data), row.names(tmb))
data=data[sameSample,,drop=F]
tmb=tmb[sameSample,,drop=F]
rt=cbind(data, tmb)

#?????Է???
x=as.numeric(rt[,gene])
y=log2(as.numeric(rt[,"total_perMB_log"])+1)
df1=as.data.frame(cbind(x,y))
corT=cor.test(x, y, method="spearman")
p1=ggplot(df1, aes(x, y)) + 
			xlab(paste0(gene, " expression"))+ylab("Tumor mutation burden")+
			geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
			stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))

#??????????ͼ??
pdf(file="cor.pdf",width=5,height=5)
print(p2)
dev.off()



