#install.packages("survival")
# if (!require("BiocManager", quietly = TRUE))
#         install.packages("BiocManager")
# BiocManager::install("limma")
# install.packages("tidyverse") # very slow, need 30 sec to finish 
# conda install conda-forge::r-tidyverse

library(tidyverse)
library(survival)
library(limma)
#setwd() 

expFile="diffGeneExp.txt"

cliFile="time.txt"  

rt2=read.table(expFile, header=T, check.names=F, row.names=1)

exp_data_T = rt2%>% dplyr::select(str_which(colnames(.), "-01A")) #

nT = ncol(exp_data_T) 

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$futime=cli$futime/365

group1=sapply(strsplit(colnames(rt2),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)
conNum=length(group1[group1==1])       
treatNum=length(group1[group1==0])     
Type=c(rep(1,conNum), rep(2,treatNum))


#ɾ��������Ʒ
tumorData=exp_data_T
tumorData=as.matrix(tumorData)
tumorData=t(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)
#���ݺϲ���������
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt2=cbind(cli, data)
#����ϲ��������
outTab=cbind(ID=row.names(rt2), rt2)
outTab=outTab[,-ncol(outTab)]


#�������
pFilter=0.05                                                                 #�����Թ�������
rt=outTab    #��ȡ�����ļ�                                                   #�������Ϊ��λ������30������Ϊ��λ������365
outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,4:ncol(rt)])){
	   if(sd(rt[,gene])<0.1){
	     next}
	   a=rt[,gene]<=median(rt[,gene])
	  
		 rt1=rt[a,]
	   b=setdiff(rownames(rt),rownames(rt1))
	   rt2=rt[b,]
	   surTab1=summary(survfit(Surv(futime, fustat) ~ 1, data = rt1))
	   surTab2=summary(survfit(Surv(futime, fustat) ~ 1, data = rt2))
	   survivalTab1=cbind(time=surTab1$time, surv=surTab1$surv,lower=surTab1$lower,upper=surTab1$upper)
	   survivalTab1=survivalTab1[survivalTab1[,"time"]<5,]
	   if(class(survivalTab1)=="matrix"){
	     survivalTab1=survivalTab1[nrow(survivalTab1),]
	   }
	   survivalTab2=cbind(time=surTab2$time, surv=surTab2$surv,lower=surTab2$lower,upper=surTab2$upper)
	   survivalTab2=survivalTab2[survivalTab2[,"time"]<5,]
	   if(class(survivalTab2)=="matrix"){
	     survivalTab2=survivalTab2[nrow(survivalTab2),]
	   }
	   fiveYearsDiff=abs(survivalTab1["surv"]-survivalTab2["surv"])

     #km����
	   diff=survdiff(Surv(futime, fustat) ~a,data = rt)
	   pValue=1-pchisq(diff$chisq,df=1)
	   fit=survfit(Surv(futime, fustat) ~ a, data = rt)
	   #cox����
	   cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
	   coxSummary = summary(cox)
	   coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	
	   if((pValue<pFilter) & (coxP<pFilter) & (fiveYearsDiff>0.15)){
	       sigGenes=c(sigGenes,gene)
	       outTab=rbind(outTab, 
	                    cbind(gene=gene,
	                          KM=pValue,
	                          HR=coxSummary$conf.int[,"exp(coef)"],
	                          HR.95L=coxSummary$conf.int[,"lower .95"],
	                          HR.95H=coxSummary$conf.int[,"upper .95"],
			                      coxPvalue=coxP) )
		 }
}
write.table(outTab,file="survival.xls",sep="\t",row.names=F,quote=F)    #��������pֵ�����ļ�
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="surSigExp.txt",sep="\t",row.names=F,quote=F)

