
library(limma)            #引用包
expFile="combined_RNAseq_FPKM.txt"      #表达数据文件
geneFile="gene.txt"       #基因列表文件


#读取输入文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#读取基因列表文件,提取铜死亡相关基因的表达量
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#输出结果
outTab=rbind(ID=colnames(geneExp),geneExp)
write.table(outTab, file="铜死亡基因表达量.txt", sep="\t", quote=F, col.names=F)

#####共表达分析
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
RNA=data[,group==0]
conNum=length(group[group==1])       #正常组样品数目
treatNum=length(group[group==0])     #肿瘤组样品数目
sampleType=c(rep(1,conNum), rep(2,treatNum))


rt1=rt1=read.table("铜死亡基因表达量.txt", header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
cuproptosis=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
cuproptosis=avereps(cuproptosis)
cuproptosis=cuproptosis[rowMeans(cuproptosis)>0,]

#删掉正常样品
group=sapply(strsplit(colnames(cuproptosis),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
cuproptosis=cuproptosis[,group==0]

#相关性检验
corFilter=0.3           
pvalueFilter=0.05  
outTab=data.frame()
for(i in row.names(RNA)){
  if(sd(RNA[i,])>0.1){
    test=wilcox.test(data[i,] ~ sampleType)
    if(test$p.value<0.05){
      for(j in row.names(cuproptosis)){
        x=as.numeric(RNA[i,])
        y=as.numeric(cuproptosis[j,])
        corT=cor.test(x,y)
        cor=corT$estimate
        pvalue=corT$p.value
        if((cor>corFilter) & (pvalue<pvalueFilter)){
          outTab=rbind(outTab,cbind(Cuproptosis=j,RNA=i,cor,pvalue,Regulation="postive"))
        }
        if((cor< -corFilter) & (pvalue<pvalueFilter)){
          outTab=rbind(outTab,cbind(Cuproptosis=j,RNA=i,cor,pvalue,Regulation="negative"))
        }
      }
    }
  }
}

#输出相关性的结果
write.table(file="Result.txt",outTab,sep="\t",quote=F,row.names=F)

#提取铜死亡相关RNA的表达量
cuproptosisRNA=unique(as.vector(outTab[,"RNA"]))
cuproptosisRNAexp=data[cuproptosisRNA,]
cuproptosisRNAexp=rbind(ID=colnames(cuproptosisRNAexp), cuproptosisRNAexp)
write.table(cuproptosisRNAexp,file="cuproptosisExp.txt",sep="\t",quote=F,col.names=F)
