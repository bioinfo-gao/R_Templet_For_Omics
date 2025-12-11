library(methylSig)
library(gplots);
library(pheatmap);
library(lattice);

file_non=list.files(".",".*cov10_non.fmt")
file_res=list.files(".",".*cov10.fmt")
n1=length(file_non)
n2=length(file_res)
fileList=c(file_non,file_res)
smp1=paste0("non",c(1:n1))
smp2=paste0("res",c(1:n2))
sample.id = c(smp1,smp2)
treatment = c(rep(0,n1),rep(1,n2))
meth <- methylSigReadData(fileList, sample.ids = sample.id, assembly = "hg19",
                          treatment = treatment, context = "CpG", destranded=TRUE, num.cores=4, 
                          minCount=10,  quiet=TRUE)  ######## minCount to define the minimum coverage

### Call DMR: Try different window size to get the corresponoding No.DMR and check which size give the most DMR ###
for(wd in seq(25,100,25))
{
    methTile = methylSigTile(meth,win.size = wd)
    myDiffSigbothTile = methylSigCalc(methTile, groups=c(1,0), 
                                      min.per.group=2,local.disp=TRUE, winsize.disp=200, num.cores=4)
    save(methTile,myDiffSigbothTile,file=paste0("Pre_non_vs_res_Diff-cov10_mingroup2_wdsz",wd,".rda"))
}


########### generate plots ######################################

sample.id = c(paste0("N",c(1,2,3)),paste0("R",c(1:7)))

mtdf=c()
wd=25
flnm=paste0("Pre_non_vs_res_Diff-cov10_mingroup2_wdsz",wd)
dt=paste0(flnm,".rda")
load(dt)
idall=methTile@data.ids
covs=methTile@data.coverage
nmcs=methTile@data.numCs
mi=apply(methTile[, "coverage"],1,min)
tt=(mi>0)
x = methTile[,"numCs"]/methTile[, "coverage"]
idall2=idall[tt]
colnames(x) = meth@sample.ids
colnames(x)=sample.id
x2=x
x_sd=apply(x2,1,function(x) sd(x,na.rm=T))
# save(x,x_sd,file="CpG_cov10_wdsz25_SD.rda")

# h=hist(x_sd)
# xx=h$mids
# yy=h$counts

################ generate a denstiy polt #############################
d=density(x_sd)
pdf("CpG-sd_densityPlot_naRowRM.pdf")
plot(d$x,d$y,type="l",col="blue",xlab="SD",ylab="Density")
dev.off()

######## Clustering samples method 1: hierarchical clulstering ######
ct=0.2 ########## choose Methylation regions standard variation higher than a cut off to cluster samples
vv=(x_sd>ct)
x2=x[vv,]

rownames(x2) =rep(NA, NROW(x2))
x2[is.na(x2)]=0

corrmat= cor(x2, use="pairwise.complete.obs")

pdf(file = paste("Carraway_sampleCorheatmap-Dens",ct,".pdf",sep=""))
heatmap.2(corrmat, na.rm=T, breaks=100,#labRow=F,
          hclustfun = function(x) hclust(x,method="single"),
          col="bluered",trace="none", symm=T, keysize=1,density.info="none")

dev.off() 

################ Clustering method 2: PCA analysis ##################################
library("devtools")
# install_github("kassambara/factoextra")
library("factoextra")
mm=t(x2)
res.pca <- prcomp(mm,  scale = FALSE)
smps=res.pca$x[,c(1,2)]
res.kms=kmeans(smps,2)
cls=c(rep("non",3),rep("res",7))

labls=methTile@sample.ids
pdf(file = paste("PCA_analysis-Dens",ct,".pdf",sep=""))
p <- fviz_pca_ind(res.pca, label="none",  geom="point", habillage=cls,pointsize = 2,
                  addEllipses=FALSE)  ##, ellipse.level=0.9
p
dev.off()
