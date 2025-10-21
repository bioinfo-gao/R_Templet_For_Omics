

#install.packages("survivalROC")

library(survivalROC)
library(survival)
library(survminer)
library(timeROC)
#setwd("D:\\biowolf\\80geneFilter\\12.RocFilter")                      #工作目录（需修改）
rocFilter=0                                                                     #ROC过滤值
rt=read.table("surSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)      #读取输入文件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
	   roc=timeROC(T=rt$futime, 
	                   delta=rt$fustat, 
	                   marker = rt[,i], 
	                   cause=1, 
	                   weighting='aalen',time=5,ROC=TRUE)
	   if(roc$AUC[2]>rocFilter){
	       sigGenes=c(sigGenes,i)
	       outTab=rbind(outTab,cbind(gene=i,AUC=roc$AUC[2]))
	   }
}
write.table(outTab,file="ROC.xls",sep="\t",row.names=F,quote=F)    #输出基因和p值表格文件
rocSigExp=rt[,sigGenes]
rocSigExp=cbind(id=row.names(rocSigExp),rocSigExp)
write.table(rocSigExp,file="rocSigExp.txt",sep="\t",row.names=F,quote=F)

#可视化
gene=colnames(rt)[3]
ROC_rt=timeROC(T=rt$time, delta=rt$event,
               marker=rt[,gene], cause=1,
               weighting='aalen',
               times=c(1,3,5,10), ROC=TRUE)
pdf(file=paste0(gene, ".ROC.pdf"), width=5, height=5)
plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=10,col='yellow',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3])),
         paste0('AUC at 10 years: ',sprintf("%.03f",ROC_rt$AUC[4]))),
       col=c("green",'blue','red','yellow'),lwd=2,bty = 'n')
dev.off()