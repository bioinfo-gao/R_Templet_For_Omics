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
getwd()
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



class(exp_data_T)
tumorData=exp_data_T
tumorData[1:2, 1:3]
dim(tumorData)

# tumorData=as.matrix(tumorData)
tumorData=t(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt2=cbind(cli, data)
head(rt2)
rt2[1:2, 1:5]
outTab=cbind(ID=row.names(rt2), rt2)

#outTab[1:10,ncol(outTab), drop = F]
outTab[,-ncol(outTab)] # why do this step?  <<=== ??
# outTab[1:2, 1:5]
# outTab[1:2, 1:5]


#
pFilter=0.05                                                                 
rt=outTab    #
outTab=data.frame()
sigGenes=c("futime","fustat")
rt[1:5, (ncol(rt)-4):ncol(rt)]

#class(survivalTab1)
for(gene in colnames(rt[,4:ncol(rt)])){
    #gene = "PSMB3"
    
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
    
    if(is.matrix(survivalTab1)){ #  BAD code   class(survivalTab1) [1] "matrix" "array"   if(class(survivalTab1)=="matrix"){
        survivalTab1=survivalTab1[nrow(survivalTab1),]
    }
    
    survivalTab2=cbind(time=surTab2$time, surv=surTab2$surv,lower=surTab2$lower,upper=surTab2$upper)
    survivalTab2=survivalTab2[survivalTab2[,"time"]<5,]
    
    if(is.matrix(survivalTab2)){
        survivalTab2=survivalTab2[nrow(survivalTab2),]
    }
    
    fiveYearsDiff=abs(survivalTab1["surv"]-survivalTab2["surv"]) # <<====
    
    
    diff=survdiff(Surv(futime, fustat) ~a,data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    fit=survfit(Surv(futime, fustat) ~ a, data = rt)
    #cox
    cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt) # <<============== SURVIVAL DIFFERENCE 
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

write.table(outTab,file="survival.xls",sep="\t",row.names=F,quote=F)    #

# KM	HR	HR.95L	HR.95H	coxPvalue
# 1. KM (Kaplan-Meier)
# 解释： KM < 0.05 通常表示高表达组和低表达组的生存曲线分离显著，即基因表达水平分组对生存时间有显著影响。
# 2. HR (Hazard Ratio) - 风险比
# HR > 1： 基因表达水平越高，风险越大（生存时间越短）。通常是不良预后指标。
# HR < 1： 基因表达水平越高，风险越小（生存时间越长）。通常是良好预后指标。
# HR = 1： 基因表达水平与生存风险无关。
# 3. HR.95L 和 HR.95H
# HR.95L	Hazard Ratio 95% Confidence Interval Lower Bound	风险比 95% 置信区间的下限。
# HR.95H	Hazard Ratio 95% Confidence Interval Upper Bound	风险比 95% 置信区间的上限。
# 4. coxPvalue
# coxPvalue < 0.05（两者是同一件事情的不同表现形式）。
# KM P-value 是单因素分组检验（高 vs 低），coxPvalue 是回归模型检验（连续变量或分组变量的风险贡献），两者通常互相印证，但 Cox P 值是更严谨的量化指标。
# ZG 注： 基因表达是连续的， cox 更靠谱
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="surSigExp.txt",sep="\t",row.names=F,quote=F)

