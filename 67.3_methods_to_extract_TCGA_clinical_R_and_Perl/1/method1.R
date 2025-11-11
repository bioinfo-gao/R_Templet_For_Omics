
library(readr)

setwd("C:/Users/zhen-/Code/R_code/R_Omics_DS/67.3_methods_to_extract_TCGA_clinical_R_and_Perl")

clinical=read_tsv("./1/clinical.cart.2023-02-03/clinical.tsv")

ID=clinical$case_submitter_id

age=clinical$age_at_index

gender=clinical$gender

time=clinical$days_to_death

status=clinical$vital_status

pathologicT=clinical$ajcc_pathologic_t
pathologicM=clinical$ajcc_pathologic_m
pathologicN=clinical$ajcc_pathologic_n

pathologicStage=clinical$ajcc_pathologic_stage


TCGA_merge=cbind(ID,
                    age,
                    gender,
                    time,
                    status,
                    pathologicT,
                    pathologicM,
                    pathologicN,
                    pathologicStage)

TCGA_merge[which(TCGA_merge=="'--")]=NA
TCGA_clinical=na.omit(TCGA_merge)
TCGA_clinical=as.data.frame(TCGA_clinical)

duplicated(TCGA_clinical$ID)
TCGA_clinical<-TCGA_clinical[!duplicated(TCGA_clinical$ID),]


rownames(TCGA_clinical)=TCGA_clinical$ID
TCGA_clinical=TCGA_clinical[,2:ncol(TCGA_clinical)]
write.csv(TCGA_clinical,file = "TCGA_clinical.csv",quote = F)
