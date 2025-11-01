
library(readr)
clinical=read_tsv("clinical.tsv")

#提取样品ID
ID=clinical$case_submitter_id
#提取年龄
age=clinical$age_at_index
#提取性别
gender=clinical$gender
#提取生存时间
time=clinical$days_to_death
#提取生存状态
status=clinical$vital_status
#提取TMN分期
pathologicT=clinical$ajcc_pathologic_t
pathologicM=clinical$ajcc_pathologic_m
pathologicN=clinical$ajcc_pathologic_n
#提取stage分期
pathologicStage=clinical$ajcc_pathologic_stage

#合并信息
TCGA_merge=cbind(ID,
                    age,
                    gender,
                    time,
                    status,
                    pathologicT,
                    pathologicM,
                    pathologicN,
                    pathologicStage)
#删除缺失值
TCGA_merge[which(TCGA_merge=="'--")]=NA
TCGA_clinical=na.omit(TCGA_merge)
TCGA_clinical=as.data.frame(TCGA_clinical)
#删除重复ID
duplicated(TCGA_clinical$ID)
TCGA_clinical<-TCGA_clinical[!duplicated(TCGA_clinical$ID),]

#导出文件
rownames(TCGA_clinical)=TCGA_clinical$ID
TCGA_clinical=TCGA_clinical[,2:ncol(TCGA_clinical)]
write.csv(TCGA_clinical,file = "TCGA_clinical.csv",quote = F)
