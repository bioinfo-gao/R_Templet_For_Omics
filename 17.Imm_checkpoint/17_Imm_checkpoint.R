setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/17.Imm_checkpoint")
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(corrplot)

pFilter=0.05             
geneName="TSPAN6"           
expFile="combined_RNAseq_counts.txt"     
#expFile="../01.New_TCGA/combined_RNAseq_counts.txt"      
geneFile="gene.txt"       


rt=read.table(expFile, header=T, sep="\t", check.names=F)
#rt=as.matrix(rt)
rt[1:2, 1:5]

rownames(rt)=make.unique(rt[,1]) # 1
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
	
e=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data), as.vector(gene[,1]))
data=t(data[c(geneName, sameGene),])
data=log2(data+1)

#ɾ?up=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=t(avereps(data))

#???s.numeric(data[geneName,])
outTab=data.frame()
for(i in sameGene){
	if(i==geneName){next}
    y=as.numeric(data[i,])
	corT=cor.test(x, y, method = 'pearson')
	cor=corT$estimate
	pvalue=corT$p.value
	if(pvalue<pFilter){
		outTab=rbind(outTab, cbind(Query=geneName, Gene=i, cor, pvalue))
	}
}
#???te.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)

#???a=t(data[c(geneName, as.vector(outTab[,2])),])
M=cor(data)

#???f(file="corpot.pdf",width=15,height=15)
corrplot(M,
         order="original",
         method = "color",
         number.cex = 0.7,
         addCoef.col = "black",
         diag = TRUE,
         tl.col="black",
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()


