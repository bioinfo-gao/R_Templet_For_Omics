# if not runnning using bsub:
#
# R -e 'library(DoGs);library(AtacSeq);AtacSeq:::installAtacSeq("hg19","/scratch/projects/bbc/aiminy_project/AtacSeq")'
#
#
installOtherPackage <- function(){
    
    library("DESeq2")
    library("RColorBrewer")
    library("gplots")
    library("dplyr")
    library("BiocParallel")
    library("ggplot2")
    library(data.table)
    library(stringr)
    
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    
    # source("https://bioconductor.org/biocLite.R")
    # biocLite("ChIPseeker")
    # biocLite("GenomeInfoDbData")
    # biocLite("DO.db")
    # biocLite("GO.db")
    # biocLite("GenomicAlignments")
    # biocLite("DelayedArray")
    # 
    # install.packages("bit64")
    # install.packages("blob")
    # install.packages("colorspace")
    # 
    # install.packages("roxygen2")
    # install.packages("devtools")
    # install_github("kassambara/easyGgplot2",force = TRUE)
    
    library(ChIPseeker)
    library(devtools)
    library(easyGgplot2)
}

installAtacSeq <- function(genome, DATA.DIR){
    
    if (!dir.exists(DATA.DIR)) {
        dir.create(DATA.DIR, recursive = TRUE)
        }
    
    #export PYTHONPATH=~/miniconda3/lib/python3.6;
    
    cmd0 = "bash ~/atac_dnase_pipelines/install_genome_data.sh"
    cmd1 = paste(cmd0,genome,DATA.DIR,collapse = " ")
    print(cmd1)
    
    system(cmd1)
    
}

#R -e 'library(DoGs);library(AtacSeq);AtacSeq:::submitJob("hg19","/scratch/projects/bbc/aiminy_project/AtacSeq")'

submitJob <- function(genome,DATA.DIR){
    
    if (!dir.exists(DATA.DIR))
    {
        dir.create(DATA.DIR, recursive = TRUE)
    }
    
    job.name <- "AtacSeq"
    
    Rfun1 <- 'library(AtacSeq);re <- AtacSeq:::installAtacSeq('
    
    Rinput <- paste0('\\"',genome,'\\",',
                     '\\"',DATA.DIR,'\\"')
    Rfun2 <- ')'
    
    Rfun <- paste0(Rfun1,Rinput,Rfun2)
    
    cmd.gff <- DoGs:::createBsubJobArrayRfun(Rfun,job.name,wait.job.name=NULL)
    
    system(cmd.gff)
    
}

testAtacSeq <- function(input.file.dir,output.file.dir){
    
    file.1 <- list.files(input.bam.dir,pattern=".bam$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    r.lib <- Sys.getenv("R_LIBS_USER")
    
    if (!dir.exists(dirname(output.file.dir)))
    {
        dir.create(dirname(output.file.dir), recursive = TRUE)
    }
    
    input <- paste(file.1,collapse = " ")
    output <- output.file.dir
    
    cmd1 <- paste("sh",file.path(r.lib,"AtacSeq/bin/bash/testAtacSeq.sh"),input,output,sep = " ")
    
    print(cmd1)
    
    system(cmd1)
    
}

#R -e 'library(DoGs);library(AtacSeq);AtacSeq:::submitJob("hg19","/scratch/projects/bbc/aiminy_project/AtacSeq")'

submitJob4testAtacSeq <- function(input.file.dir,output.file.dir){
    
    job.name <- "RunAtac"
    
    Rfun1 <- 'library(AtacSeq);re <- AtacSeq:::testAtacSeq('
    
    Rinput <- paste0('\\"',input.file.dir,'\\",',
                     '\\"',output.file.dir,'\\"')
    Rfun2 <- ')'
    
    Rfun <-paste0(Rfun1,Rinput,Rfun2)
    
    cmd.gff <- DoGs:::createBsubJobArrayRfun(Rfun,job.name,wait.job.name=NULL)
    
    system(cmd.gff)
    
}


testAtacSeq0 <- function() {
    
    #file.1 <- list.files(input.bam.dir,pattern=".bam$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    r.lib <- Sys.getenv("R_LIBS_USER")
    
    #if (!dir.exists(dirname(output.file.dir)))
    #{
    #  dir.create(dirname(output.file.dir), recursive = TRUE)
    #
    #}
    
    #input <- paste(file.1,collapse = " ")
    #output <- output.file.dir
    
    #cmd0="export PATH=/usr/bin:$PATH"
    
    cmd1 <- paste("sh",file.path(r.lib,"AtacSeq/bin/bash/testAtacSeq.sh"),sep = " ")
    
    print(cmd1)
    
    system(cmd1)
    
}

#R -e 'library(DoGs);library(AtacSeq);AtacSeq:::submitJob4testAtacSeq0()'

submitJob4testAtacSeq0 <- function(){
    
    job.name <- "RunAtac"
    
    Rfun1 <- 'library(AtacSeq);re <- AtacSeq:::testAtacSeq0()'
    
    #Rinput <- paste0('\\"',input.file.dir,'\\",',
    #          '\\"',output.file.dir,'\\"')
    #Rfun2 <- ')'
    
    Rfun <- Rfun1
    
    cmd.gff <- DoGs:::createBsubJobArrayRfun(Rfun,job.name,wait.job.name=NULL)
    
    #cmd0="module load java/1.8.0_60;export _JAVA_OPTIONS=-Xms256M -Xmx728M -XX:ParallelGCThreads=1;module unload python/2.7.3;
    #unset PYTHONPATH;
    #export PYTHONPATH=/nethome/axy148/anaconda3/envs/nothing/lib/python2.7/site-packages"
    
    cmd0 = "module load java/1.8.0_60;export _JAVA_OPTIONS=-Xms256M -Xmx728M -XX:ParallelGCThreads=1;module unload python/2.7.3"
    # unset PYTHONPATH;
    #  export PYTHONPATH=/nethome/axy148/anaconda3/envs/nothing/lib/python2.7/site-packages"
    
    #cmd0="module load java/1.8.0_60;export _JAVA_OPTIONS=-Xms256M -Xmx728M -XX:ParallelGCThreads=1;
    #module unload python/2.7.3;
    #unset PYTHONPATH;
    #source activate python2;
    #export PYTHONPATH=/nethome/axy148/anaconda3/envs/python2/lib/python2.7/site-packages"
    
    cmd.gff2 = paste(cmd0,cmd.gff,sep=";")
    
    system(cmd.gff2)
    
}

# R -e 'library(DoGs);library(AtacSeq);AtacSeq:::testAtacSeqNonCluster("/media/aiminyan/DATA/AtacSeq_Input","/media/aiminyan/DATA/AtacSeq_Output")'

testAtacSeqNonCluster <- function(input.fastq.dir,output) {
    
    file.1 <- list.files(input.fastq.dir,pattern="*.fq.gz$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    r.lib <- Sys.getenv("R_LIBS_USER")
    
    if (!dir.exists(dirname(output)))
    {
        dir.create(dirname(output), recursive = TRUE)
    }
    
    #input <- paste(file.1,collapse = " ")
    #output <- output.file.dir
    
    #cmd0="export PATH=/usr/bin:$PATH"
    
    cmd1 <- paste("sh",file.path(r.lib,"AtacSeq/bin/bash/testAtacSeq1.sh"),file.1[1],file.1[2],output,sep = " ")
    
    print(cmd1)
    
    system(cmd1)
    
}

CutStringByNFromEnd <- function(SequenceSampleID,n){
    
    SequenceSampleID2 <- lapply(SequenceSampleID,function(u,n){
        m <- nchar(u)
        ff <- strtrim(u,m-n)
        ff
    },n)
    
    SequenceSampleID3 <- unlist(SequenceSampleID2)
    
    return(SequenceSampleID3)
}

generateFq <- function(f,tm,pattern,TorC){
    
    f1 <- f[grep(pattern = pattern,f$Sample.Name),]
    f2 <- f1[which(f1$Hours==tm),]
    s <- unique(f2$Sample)
    
    x <- lapply(1:length(s),function(u,s,f2){
        
        v <- paste(paste0(TorC,u,"_1"),f2[which(f2$Sample==s[u]),]$fq.files[1],paste0(TorC,u,"_2"),f2[which(f2$Sample==s[u]),]$fq.files[2],sep=" ")
        v
    },s,f2)
    
    xx <- paste(do.call(c,x),collapse = " ")
    
    xx
    
}

# sample.info = "~/pegasus/Project/Alejandro_AtacSeq/ATACSeq_sample_mapping.csv"
# title = "IL-2vsPBS"
# species = "mm10"
# input.fq.dir = "~/pegasus/Project/Alejandro_atac/DATA/Formatted"
# output = "~/pegasus/Project/Alejandro_AtacSeq"
#
# sample.info = "~/AiminDropbox/Alejandro_AtacSeq_uploaded_2/ATACSeq_sample_mapping.csv"
# title = "IL-2vsPBS"
# species = "mm10"
# input.fq.dir = "~/AiminDropbox/Alejandro_atac_fq_input/DATA/Formatted"
# output = "~/pegasus/Project/Alejandro_AtacSeq"
#
#
# AtacSeq:::testAtacSeqNonCluster2(sample.info,input.fq.dir,title,species,output)
#
testAtacSeqNonCluster2 <- function(sample.info,input.fq.dir,title,species,output)
{
    
    f <- read.csv(sample.info)
    
    fq.files <- list.files(path = input.fq.dir, pattern= "*.fastq$",full.names = TRUE, recursive = TRUE)
    
    fq.files.2 <- cbind.data.frame(fq.files,CutStringByNFromEnd(as.character(basename(fq.files)),12),stringsAsFactors = F)
    
    colnames(fq.files.2)[2] =  "BaseSpace.Sample.ID"
    
    print(f)
    
    print(fq.files.2)
    
    ff <- merge(fq.files.2,f,by="BaseSpace.Sample.ID")
    
    print(ff)
    
    #fq1 <- ff[which((ff$Hours==16.00)&&(grep("PBS",ff$Sample.Name))),]$fq.files
    
    #fq2 <- ff[which((ff$Hours==16.00)&&(grep("IL2",ff$Sample.Name))),]$fq.files
    
    #print(fq1)
    
    #print(fq2)
    
    
    
    if (!dir.exists(dirname(output)))
    {
        dir.create(dirname(output), recursive = TRUE)
    }
    
    cmd0 = "unset PYTHONPATH"
    cmd1 = paste0("TITLE=",title)
    cmd2 = paste0("SPECIES=",species)
    cmd3 = paste0("WORK=",file.path(output,"$TITLE"))
    cmd4 = "mkdir -p $WORK;cd $WORK"
    cmd5 = "$HOME/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8"
    
    cmd6.t = generateFq(ff,16.00,"IL2","-fastq")
    cmd6.ck = generateFq(ff,16.00,"PBS","-ctl_fastq")
    
    cmd6 = paste(cmd6.t,cmd6.ck,sep = " ")
    
    cmd7 = "-enable_idr -idr_suffix -gensz mm -out_dir $WORK"
    cmd8 = paste(cmd5,cmd6,cmd7)
    
    cmd = paste(cmd0,cmd1,cmd2,cmd3,cmd4,cmd8,sep=";")
    
    cat(cmd,"\n")
    
    system(cmd)
    
    #return(ff)
    
}

# sample.info = "~/pegasus/Project/Alejandro_AtacSeq/ATACSeq_sample_mapping.csv"
# species = "mm10"
# input.fq.dir = "~/pegasus/Project/Alejandro_atac/DATA/Formatted"
# output = "~/pegasus/Project/Alejandro_AtacSeq"
#
# AtacSeq:::testAtacSeqNonCluster3(sample.info,input.fq.dir,species,output)
#
testAtacSeqNonCluster3 <- function(sample.info,input.fq.dir,species,output) {
    
    f <- read.csv(sample.info)
    
    fq.files <- list.files(path = input.fq.dir, pattern= "*.fastq$",full.names = TRUE, recursive = TRUE)
    
    fq.files.2 <- cbind.data.frame(fq.files,CutStringByNFromEnd(as.character(basename(fq.files)),12),stringsAsFactors = F)
    
    colnames(fq.files.2)[2] =  "BaseSpace.Sample.ID"
    
    #print(f)
    
    #print(fq.files.2)
    
    ff <- merge(fq.files.2,f,by="BaseSpace.Sample.ID")
    
    print(ff)
    
    t <- unique(ff$Hours)
    s <- unique(as.character(ff$Sample.Name))
    
    ss <- unique(as.character(sapply(strsplit(s,"_"), `[`, 1)))
    #ss <- strsplit(s,"_")[1]
    
    #print(t)
    #print(ss)
    
    #fq1 <- ff[which((ff$Hours==16.00)&&(grep("PBS",ff$Sample.Name))),]$fq.files
    
    #fq2 <- ff[which((ff$Hours==16.00)&&(grep("IL2",ff$Sample.Name))),]$fq.files
    
    #print(fq1)
    
    #print(fq2)
    
    if (!dir.exists(dirname(output)))
    {
        dir.create(dirname(output), recursive = TRUE)
    }
    
    
    cmd0 = "unset PYTHONPATH"
    
    cmd.l <- lapply(t, function(u,ss,ff,species,output){
        
        lapply(ss,function(y,ff,species,output,u){
            
            
            title= paste0(y,"-at-",u)
            
            if(title!="IL2-at-16")
            {
                cmd1 = paste0("TITLE=",title)
                cmd2 = paste0("SPECIES=",species)
                cmd3 = paste0("WORK=",file.path(output,"$TITLE"))
                cmd4 = "mkdir -p $WORK;cd $WORK"
                cmd5 = "$HOME/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8"
                
                cmd6 = generateFq(ff,u,y,"-fastq")
                #cmd6.ck = generateFq(ff,16.00,"PBS","-ctl_fastq")
                
                #cmd6 = paste(cmd6.t,cmd6.ck,sep = " ")
                
                cmd7 = "-enable_idr -idr_suffix -gensz mm -out_dir $WORK"
                cmd8 = paste(cmd5,cmd6,cmd7)
                
                cmd = paste(cmd0,cmd1,cmd2,cmd3,cmd4,paste0(cmd8,"&"),sep=";")
                
                
                cat(title,"\n")
                cat(cmd,"\n")
                cat("\n")
                
                system(cmd)
                
            }
            
            
        },ff,species,output,u)
        
    },ss,ff,species,output)
    #return(ff)
    
}

CutStringByNFromEnd <- function(SequenceSampleID,n){
    
    SequenceSampleID2 <- lapply(SequenceSampleID,function(u,n){
        m <- nchar(u)
        ff <- strtrim(u, m-n)
        ff
    }, n)
    
    SequenceSampleID3 <- unlist(SequenceSampleID2)
    
    return(SequenceSampleID3)
}

# extractAtacSeq("/Users/axy148/.CMVolumes/AiminDropbox/Alejandro_AtacSeq_uploaded_2","~/Dropbox/Alejandro_AtacSeq_out")
#
extractAtacSeq <- function(input.atac.dir,output.file.dir){
    
    if (!dir.exists(dirname(output.file.dir)))
    {
        dir.create(dirname(output.file.dir), recursive = TRUE)
    }
    #rep1-pr.IDR0.1.filt.narrowPeak.gz
    #peak/macs2/idr/pseudo_reps/rep
    file.1 <- list.files(input.atac.dir,pattern="*pr.IDR0.1.filt.narrowPeak.gz$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    file.2 <- file.1[grep("pseudo_reps",file.1)]
    
    file.3 <- file.2[-grep("pooled_pseudo_reps_0.1",file.2)]
    
    
    file.4 <- cbind.data.frame(ID = CutStringByNFromEnd(basename(file.3),34),file.3,stringsAsFactors=FALSE)
    
    file.id <- unique(file.4$ID)
    
    cmd.l <- lapply(file.id, function(u,file.4){
        
        cat(u,"\n")
        
        if(u == "IL-2vsPBS-2"){x = "IL2-at-16"
        cmd0 <- paste("gunzip -c",paste(file.4[file.4$ID == u,]$file.3,collapse = " "),">",file.path(output.file.dir,paste0(x,".idr.concatenated.bed")),sep = " ")
        cmd1 <- paste("sort -k1,1 -k2,2n",file.path(output.file.dir,paste0(x,".idr.concatenated.bed")),">",file.path(output.file.dir,paste0(x,".idr.concatenated.sorted.bed")),sep =" ")
        cmd2 <- paste("mergeBed -i", file.path(output.file.dir,paste0(x,".idr.concatenated.sorted.bed")),">",file.path(output.file.dir,paste0(x,".ppr.merged.bed")),sep=" ")
        cmd <- paste(cmd0,"wait",cmd1,"wait",cmd2,sep = ";")
        }else{
            cmd0 <- paste("gunzip -c",paste(file.4[file.4$ID == u,]$file.3,collapse = " "),">",file.path(output.file.dir,paste0(u,".idr.concatenated.bed")),sep = " ")
            cmd1 <- paste("sort -k1,1 -k2,2n",file.path(output.file.dir,paste0(u,".idr.concatenated.bed")),">",file.path(output.file.dir,paste0(u,".idr.concatenated.sorted.bed")),sep =" ")
            cmd2 <- paste("mergeBed -i", file.path(output.file.dir,paste0(u,".idr.concatenated.sorted.bed")),">",file.path(output.file.dir,paste0(u,".ppr.merged.bed")),sep=" ")
            cmd <- paste(cmd0,"wait",cmd1,"wait",cmd2,sep = ";")
        }
        
        print(cmd)
        
        system(cmd)
        
    },file.4)
    
    # print(file.4)
    
}

generate.ppr.merged.bed <- function() {
    ppr.merged.bed <- system("ls -lrth /Users/axy148/Dropbox/Alejandro_AtacSeq_out/*.merged.bed | cut -d' ' -f12",intern = T)
    output.dir.pegasus <- "~/pegasus/aiminy_project/Alejandro_AtacSeq_out"
    cmd0 <- paste("cat",paste(ppr.merged.bed,collapse = " "),">",file.path(output.dir.pegasus,"cat.ppr.merged.bed"))
    print(cmd0)
    system(cmd0)
    
    cmd1 <- paste("mergeBed -i",file.path(output.dir.pegasus,"cat.ppr.merged.bed"),">",file.path(output.dir.pegasus,"ppr.merged.bed"))
    
    system(cmd1)
}


generate.union.ppr.merged.bed <- function(){
    
    ppr.merged.bed <- system("ls -lrth /Users/axy148/Dropbox/Alejandro_AtacSeq_out/*.merged.bed | cut -d' ' -f12",intern = T)
    
    output.dir.pegasus <- "~/Dropbox/Alejandro_AtacSeq_out"
    
    cmd1 <- paste("bedops -u",paste(ppr.merged.bed[1:6],collapse = " "),"| uniq - >",file.path(output.dir.pegasus,"union.ppr.merged.bed"))
    
    system(cmd1)
    
}

# input.atac.dir <- "/Users/axy148/.CMVolumes/AiminDropbox/Alejandro_AtacSeq_uploaded_2"
# output.dir <- "~/Dropbox/Alejandro_AtacSeq_out_2"
# ppr.merged.bed.dir <- "/Users/axy148/Dropbox/Alejandro_AtacSeq_out"
# getCount(input.atac.dir,ppr.merged.bed.dir,output.dir)

getCount <- function(input.atac.dir,ppr.merged.bed.dir,output.dir) {
    
    if (!dir.exists(dirname(output.dir)))
    {
        dir.create(dirname(output.dir), recursive = TRUE)
    }
    
    #file.1 <- list.files(input.atac.dir,pattern="*tn5_pooled.tagAlign.gz$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    file.1 <- list.files(input.atac.dir,pattern="*tn5.tagAlign.gz$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    file.2 <- file.1[grep("pseudo_reps",file.1)]
    
    #file.2 <- file.1[grep("pooled_pseudo_reps",file.1)]
    
    file.3 <- cbind.data.frame(ID = basename(dirname(dirname(dirname(dirname(dirname(file.2)))))),rep=basename(dirname(dirname(file.2))),ppr= basename(dirname(file.2)),fileName=basename(file.2),file.2,stringsAsFactors=FALSE)
    
    file.4 <- file.3[-which(file.3$ID %in% c("IL-2vsPBS")),]
    
    file.5 <- cbind.data.frame(ID2 = paste0(file.4$ID,"_",file.4$rep,"_",file.4$ppr),file.4)
    
    file.id <- unique(file.5$ID2)
    
    cmd.l <- lapply(file.id, function(u,file.5){
        
        
        x <- file.5[which(file.5$ID2 == u),]
        
        cat(as.character(x$ID2),"\n")
        
        if(as.character(x$ID) == "IL-2-at-16"){
            
            cmd0 <- paste("gzip -d",x$file.2,"-c >",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),sep = " ")
            
            cmd1 <- paste("python ~/Dropbox/Aimin_project/AtacSeq/inst/deseq2_wrappers/workflow_for_raw_counts/shift_atac_tagalign.py",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),file.path(output.dir,paste0(x$ID2,"slopped")),sep = " ")
            
            cmd2 <- paste("bedtools coverage -counts -a",file.path(ppr.merged.bed.dir,"ppr.merged.bed"),"-b",file.path(output.dir,paste0(x$ID2,"slopped")),"| cut -f4 >",file.path(output.dir,x$ID2),sep = " ")
            
            cmd <- paste(cmd0,"wait",cmd1,"wait",cmd2,sep = ";")
            
        }else
        {
            cmd0 <- paste("gzip -d",x$file.2,"-c >",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),sep = " ")
            cmd1 <- paste("python ~/Dropbox/Aimin_project/AtacSeq/inst/deseq2_wrappers/workflow_for_raw_counts/shift_atac_tagalign.py",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),file.path(output.dir,paste0(x$ID2,"slopped")),sep = " ")
            cmd2 <- paste("bedtools coverage -counts -a",file.path(ppr.merged.bed.dir,"ppr.merged.bed"),"-b",file.path(output.dir,paste0(x$ID2,"slopped")),"| cut -f4 >",file.path(output.dir,x$ID2),sep = " ")
            cmd <- paste(cmd0,"wait",cmd1,"wait",cmd2,sep = ";")
        }
        print(cmd)
        system(cmd)
    },file.5)
}

makeCountTable <- function(output.dir) {
    file.count <- system("ls -lrth ~/Dropbox/Alejandro_AtacSeq_out/*ppr? | cut -d' ' -f13",intern = T)
    system(paste("echo",paste(basename(file.count),collapse = "\\\t"),">",file.path(output.dir,"count.txt")))
    cmd <- paste("paste",paste(file.count,collapse = " "), ">>",file.path(output.dir,"count.txt"))
    system(cmd)
}

# makeCountTable2(output.dir)
#
makeCountTable2 <- function(output.dir) {
    file.count <- system("ls -lrth ~/Dropbox/Alejandro_AtacSeq_out_2/*pr?  | awk '{print $9}'",intern = T)
    system(paste("echo",paste(basename(file.count),collapse = "\\\t"),">",file.path(output.dir,"count.txt")))
    cmd <- paste("paste",paste(file.count,collapse = " "), ">>",file.path(output.dir,"count.txt"))
    system(cmd)
}

trimGene <- function(g){
    
    y <- lapply(g, function(u){
        
        if(!is.na(str_locate(u,"-")[1])){
            x <- str_sub(u,1,str_locate(u,"-")[1]-1)
        }else
        {
            x <- u
        }
        x
    })
    
    gg <- unlist(y)
    
    ggg <- gg[!is.na(gg)]
    ggg
    
}

trimGene2 <- function(g){
    
    y <- lapply(g, function(u){
        
        if(!is.na(str_locate(u,"-")[1])){
            x <- str_sub(u,str_locate(u,"at")[2]+2,str_locate(u,"ppr")[1]-2)
        }else
        {
            x <- u
        }
        x
    })
    
    gg <- unlist(y)
    
    ggg <- gg[!is.na(gg)]
    ggg
    
}

trimGene3 <- function(g){
    
    y <- lapply(g, function(u){
        
        if(!is.na(str_locate(u,".")[1])){
            x <- str_sub(u,str_locate(u,"at")[2]+2,str_locate(u,"rep")[1]-2)
        }else
        {
            x <- u
        }
        x
    })
    
    gg <- unlist(y)
    
    ggg <- gg[!is.na(gg)]
    ggg
    
}

trimGene4 <- function(g){
    
    y <- lapply(g, function(u){
        
        if(!is.na(str_locate(u,"\\.")[1])){
            x <- str_sub(u,1,str_locate(u,"\\.")[1]-1)
        }else
        {
            x <- u
        }
        x
    })
    
    gg <- unlist(y)
    
    ggg <- gg[!is.na(gg)]
    ggg
    
}


Analysis.1 <- function() {
    counts.metadata <- cbind.data.frame(Sample = basename(file.count),Condition=trimGene(basename(file.count)),CellCyclePoint=trimGene2(basename(file.count)))
    counts.metadata[counts.metadata$Condition=="IL",2] <- "IL2"
    
    write.table(counts.metadata, file = file.path("~/Dropbox/Alejandro_AtacSeq_out","counts.metadata"), append = FALSE, quote = F, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = F,
                col.names = TRUE, qmethod = c("escape", "double"),
                fileEncoding = "")
    
    
    library("DESeq2")
    library("RColorBrewer")
    library("gplots")
    library("dplyr")
    library("BiocParallel")
    library("ggplot2")
    library(data.table)
    
    # Metadata
    
    padjCutoff <- 0.01/20
    foldChangeCutoff <- 0.10
    register(MulticoreParam(32))
    parallelFlag <- TRUE
    
    countData=data.frame(read.table("~/Dropbox/Alejandro_AtacSeq_out/count.txt",header=T,sep='\t'))
    
    colData=data.frame(read.table("~/Dropbox/Alejandro_AtacSeq_out/counts.metadata",header=T,sep='\t'))
    countData.1 <- na.omit(countData)
    
    any(is.na(countData.1))
    sum(is.na(countData.1))
    
    #rse_gene <- rse_gene[rowSums(assay()) != 0, ]
    
    colData.1 <- colData[colData$CellCyclePoint == 16.00,]
    
    g.0.75 <- as.character(colData[colData$CellCyclePoint == 0.75,]$Sample)
    g.0.75.1 <- gsub("-",".",g.0.75)
    
    countData.2 <- countData.1[,-which(colnames(countData.1) %in% g.0.75.1)]
    
    colData.1$Sample <- gsub("-",".",colData.1$Sample)
    colData.1
    coldata <- data.frame(row.names=colnames(countData.2), condition=colData.1$Condition)
    #effect of high vs low density
    #we use density & condition as factors and include the interaction term
    ddsFullCountTable<-DESeqDataSetFromMatrix(
        countData=countData.2,
        colData=coldata,
        design= ~condition)
    
    dds<-DESeq(ddsFullCountTable,parallel=parallelFlag,betaPrior=FALSE)
    
    res_density=results(dds, name="Density_HD_vs_LD")
    png("maplot.density.png")
    plotMA(res_density,ylim=c(-1,1))
    dev.off()
    
    #Treatment effect in low density samples
    res_treatment_ld= results(dds, contrast=c("Condition","IL2","PBS"))
    png("maplot.treatment.ld.png")
    plotMA(res_treatment_ld,ylim=c(-1,1))
    dev.off()
    
    #Treatment effect for high density samples
    res_treatment_hd= results(dds, list( c("Condition_treated_vs_control","DensityHD.Conditiontreated") ))
    png("maplot.treatment.hd.png")
    plotMA(res_treatment_hd,ylim=c(-1,1))
    dev.off()
    
    #THIS IS AN EXAMPLE OF A SINGLE FACTOR ANALYSIS -- SAMPLE IS THE SINGLE INDEPENDENT VARIABLE.
    
    (condition <- factor(c(rep("FP", 2), rep("PBS", 2),rep("IL2", 2))))
    
    coldata <- data.frame(row.names=colnames(countData.2), condition=colData.1$Condition)
    ddsIndividualFullCountTable<-DESeqDataSetFromMatrix(
        countData=countData.2,
        colData=coldata,
        design= ~condition)
    
    dds <- DESeq(ddsIndividualFullCountTable)
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    
    install.packages("locfit")
    rld <- rlog(dds,fitType='local')
    head(assay(rld))
    hist(assay(rld))
    library(RColorBrewer)
    (mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])
    
    sampleDists <- as.matrix(dist(t(assay(rld))))
    library(gplots)
    png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
    heatmap.2(as.matrix(sampleDists), key=F, trace="none",
              col=colorpanel(100, "black", "white"),
              ColSideColors=mycols[condition], RowSideColors=mycols[condition],
              margin=c(10, 10), main="At 16 hours")
    dev.off()
    
    rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
        require(genefilter)
        require(calibrate)
        require(RColorBrewer)
        rv = rowVars(assay(rld))
        select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        pca = prcomp(t(assay(rld)[select, ]))
        fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
        if (is.null(colors)) {
            if (nlevels(fac) >= 3) {
                colors = brewer.pal(nlevels(fac), "Paired")
            }   else {
                colors = c("black", "red")
            }
        }
        pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
        pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
        pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
        pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
        plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
        with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
        legend(legendpos, legend=levels(fac), col=colors, pch=20)
        #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
        #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
        #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
    }
    png("qc-pca.png", 1000, 1000, pointsize=20)
    rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
    dev.off()
    
    ddsIndividual<-DESeq(ddsIndividualFullCountTable,parallel=parallelFlag,betaPrior=FALSE)
    #create matrices to store the results from the PER-SAMPLE univariate comparison
    #get the individual comparisons
    conditionsToCompare <- data.frame(
        matrix(
            c('earlyG1.HD.treated', 'earlyG1.HD.controls',
              'earlyG1.LD.treated', 'earlyG1.LD.controls',
              'lateG1.HD.treated','lateG1.HD.controls',
              'lateG1.LD.treated','lateG1.LD.controls',
              'SG2M.HD.treated','SG2M.HD.controls',
              'SG2M.LD.treated','SG2M.LD.controls',
              'lateG1.HD.treated','earlyG1.HD.treated',
              'lateG1.LD.treated','earlyG1.LD.treated',
              'SG2M.HD.treated','lateG1.HD.treated',
              'SG2M.LD.treated','lateG1.LD.treated',
              'SG2M.HD.treated','earlyG1.HD.treated',
              'SG2M.LD.treated','earlyG1.LD.treated'),
            ncol = 2,
            byrow = TRUE),
        stringsAsFactors = FALSE)
    
    numCols <- 15
    numRows <- nrow(countData)
    diffMat <- matrix(, ncol = numCols, nrow = numRows)
    confidenceMatt <- matrix(,ncol=numCols, nrow=numRows)
    diffMat <- matrix(,ncol=numCols,nrow=numRows)
    diffMat[,1] = (res_density$padj<=padjCutoff)*(abs(res_density$log2FoldChange) >= foldChangeCutoff)*(res_density$log2FoldChange)
    diffMat[,2] = (res_treatment_ld$padj <= padjCutoff)*(abs(res_treatment_ld$log2FoldChange)>=foldChangeCutoff)*(res_treatment_ld$log2FoldChange)
    diffMat[,3] = (res_treatment_hd$padj <= padjCutoff)*(abs(res_treatment_hd$log2FoldChange)>=foldChangeCutoff)*(res_treatment_hd$log2FoldChange)
    
    colNameEntries=c("HD_vs_LD","LD_treated_vs_control","HD_treated_vs_control")
    for (i in 1:12){
        res <- results(ddsIndividual,
                       contrast = c("Sample",
                                    conditionsToCompare[i, 1],
                                    conditionsToCompare[i, 2]),
                       parallel = parallelFlag)
        colNameEntries=append(colNameEntries,paste(conditionsToCompare[i,1],"_vs_",conditionsToCompare[i,2],sep=""))
        #plot the MAplot
        #png(paste(conditionsToCompare[i,1],"_vs_",conditionsToCompare[i,2],".png",sep=""))
        #plotMA(res,ylim=c(-1,1))
        #dev.off()
        diffMat[ , 3 + i] <- (res$padj <= padjCutoff)*(abs(res$log2FoldChange) >=foldChangeCutoff)*(res$log2FoldChange)
    }
    # Replace NAs
    # Remember to take abs so that +1s and -1s don't cancel each other out
    diffMat[is.na(diffMat)] <- 0
    colnames(diffMat) <- colNameEntries
    
    peakNames <- read.table("ppr.merged.bed",header=TRUE,colClasses = c("character"))
    alldata=data.frame(peakNames,as.data.frame(diffMat))
    # Write diffMat to disk
    write.table(alldata,file = "diffMat_multifactor.STRINGENT.csv",sep='\t')
    
    
    #IF YOU WANT TO QC THE REPLICATES, YOU CAN PLOT RLOG-TRANSFORMED REPLICATE PAIRS AND CHECK CORRELATION
    #log-transform and plot pairs of replicates
    rld <- rlog(ddsIndividual)
    for(i in seq(1,24,2)){
        png(paste(as.character(i),"_",as.character(i+1),".png",sep = ""))
        plot( assay(rld)[, i:i+1], col="#00000020", pch=20, cex=0.3, xlab = colData[i,"Sample"],ylab=colData[i+1,"Sample"])
        dev.off()
    }
}

Analysis.2 <- function() {
    
    library("DESeq2")
    library("RColorBrewer")
    library("gplots")
    library("dplyr")
    library("BiocParallel")
    library("ggplot2")
    library(data.table)
    
    # Metadata
    
    padjCutoff <- 0.01/20
    foldChangeCutoff <- 0.10
    register(MulticoreParam(32))
    parallelFlag <- TRUE
    
    countData=data.frame(read.table("~/Dropbox/Alejandro_AtacSeq_out_2/count.txt",header=T,sep='\t'))
    
    counts.metadata <- cbind.data.frame(Sample = colnames(countData),Condition=CutStringByNFromEnd(colnames(countData),15),Stage=trimGene3(colnames(countData)),stringsAsFactors=F)
    
    counts.metadata$Condition <- trimGene4(counts.metadata$Condition)
    
    counts.metadata[which(as.character(counts.metadata$Condition) == "IL"),]$Condition <- "IL2"
    
    counts.metadata.1 <- cbind.data.frame(counts.metadata,ID=paste0(counts.metadata$Condition,".at.",str_sub(counts.metadata$Sample,str_locate(counts.metadata$Sample,"at")[,2]+2,str_length(counts.metadata$Sample))),stringsAsFactors=F)
    
    counts.metadata.1[which(counts.metadata.1$Stage == 0.75),]$Stage <- "early"
    counts.metadata.1[which(counts.metadata.1$Stage == 16),]$Stage <- "late"
    
    counts.metadata.2 <- cbind.data.frame(Sample=counts.metadata.1$ID,Condition=counts.metadata.1$Condition,Stage=counts.metadata.1$Stage,stringsAsFactors=F)
    colnames(countData) <- counts.metadata.2$Sample
    any(is.na(countData))
    
    counts.metadata.2$group <- factor(paste0(counts.metadata.2$Condition,"_",counts.metadata.2$Stage))
    
    design(dds) <- ~ group
    
    ddsFullCountTable<-DESeqDataSetFromMatrix(
        countData=countData,
        colData=counts.metadata.2,
        design= ~ group)
    
    dds<-DESeq(ddsFullCountTable,parallel=parallelFlag,betaPrior=FALSE)
    resultsNames(dds)
    
    res.late.vs.early.at.FP=results(dds, contrast = c("group","FP_late","FP_early"))
    res.late.vs.early.at.IL2=results(dds, contrast = c("group","IL2_late","IL2_early"))
    res.late.vs.early.at.PBS=results(dds, contrast = c("group","PBS_late","IL2_early"))
    
    res.FP.vs.PBS.at.early=results(dds, contrast = c("group","FP_early","PBS_early"))
    res.IL2.vs.PBS.at.early=results(dds, contrast = c("group","IL2_early","PBS_early"))
    res.FP.vs.IL2.at.early=results(dds, contrast = c("group","FP_early","IL2_early"))
    
    res.FP.vs.PBS.at.late=results(dds, contrast = c("group","FP_late","PBS_late"))
    res.IL2.vs.PBS.at.late=results(dds, contrast = c("group","IL2_late","PBS_late"))
    res.FP.vs.IL2.at.late=results(dds, contrast = c("group","FP_late","IL2_late"))
    
    res.FP_early.vs.IL2_late=results(dds, contrast = c("group","FP_early","IL2_late"))
    res.FP_late.vs.IL2_early=results(dds, contrast = c("group","FP_late","IL2_early"))
    res.FP_early.vs.PBS_late=results(dds, contrast = c("group","FP_early","PBS_late"))
    
    res.FP_late.vs.PBS_early=results(dds, contrast = c("group","FP_late","PBS_early"))
    res.IL2_early.vs.PBS_late=results(dds, contrast = c("group","IL2_early","PBS_late"))
    res.IL2_late.vs.PBS_early=results(dds, contrast = c("group","IL2_late","PBS_early"))
    
    ddsFullCountTable2<-DESeqDataSetFromMatrix(
        countData=countData,
        colData=counts.metadata.2,
        design= ~ Condition + Stage)
    
    dds2<-DESeq(ddsFullCountTable2,parallel=parallelFlag,betaPrior=FALSE)
    resultsNames(dds2)
    
    res.late.vs.early = results(dds2, contrast = c("Stage","late","early"))
    
    res.FP.vs.PBS = results(dds2, contrast = c("Condition","FP","PBS"))
    res.IL2.vs.PBS = results(dds2, contrast = c("Condition","IL2","PBS"))
    res.FP.vs.IL2 = results(dds2, contrast = c("Condition","FP","IL2"))
    
    pdf("DE_pvalue.pdf")
    par(mfrow = c(3,2))  # 3 rows and 2 columns
    
    hist(res.late.vs.early$pvalue)
    hist(res.FP.vs.PBS$pvalue)
    hist(res.IL2.vs.PBS$pvalue)
    hist(res.FP.vs.IL2$pvalue)
    
    hist(res.late.vs.early.at.FP$pvalue)
    hist(res.late.vs.early.at.IL2$pvalue)
    hist(res.late.vs.early.at.PBS$pvalue)
    
    hist(res.FP.vs.PBS.at.early$pvalue)
    hist(res.IL2.vs.PBS.at.early$pvalue)
    hist(res.FP.vs.IL2.at.early$pvalue)
    
    hist(res.FP.vs.PBS.at.late$pvalue)
    hist(res.IL2.vs.PBS.at.late$pvalue)
    hist(res.FP.vs.IL2.at.late$pvalue)
    
    hist(res.FP_early.vs.IL2_late$pvalue)
    hist(res.FP_late.vs.IL2_early$pvalue)
    hist(res.FP_early.vs.PBS_late$pvalue)
    
    hist(res.FP_late.vs.PBS_early$pvalue)
    hist(res.IL2_early.vs.PBS_late$pvalue)
    hist(res.IL2_late.vs.PBS_early$pvalue)
    
    dev.off()
    
    #Comparisions between different conditions at early stage(0.75)
    #
    de.analysis <- function(stage) {
        
        counts.metadata.early <-counts.metadata.2[counts.metadata.2$Stage == stage,1:2]
        countData.early <- countData[,which(colnames(countData) %in% counts.metadata.early$Sample)]
        
        
        ddsFullCountTable<-DESeqDataSetFromMatrix(
            countData=countData.early,
            colData=counts.metadata.early,
            design= ~Condition)
        
        dds<-DESeq(ddsFullCountTable,parallel=parallelFlag,betaPrior=FALSE)
        
        res.IL2.FP <- results(dds,contrast=c("Condition","IL2","FP"))
        res.FP.IL2 <- results(dds,contrast=c("Condition","FP","IL2"))
        
        res.IL2.PBS <- results(dds,contrast=c("Condition","IL2","PBS"))
        res.FP.PBS <- results(dds,contrast=c("Condition","FP","PBS"))
        
        res <- list(res.IL2.FP=res.IL2.FP,res.FP.IL2=res.FP.IL2,res.IL2.PBS=res.IL2.PBS,res.FP.PBS=res.FP.PBS)
        res
    }
    
    res.early <- de.analysis("early")
    res.late <- de.analysis("late")
    
    peakNames <- read.table(file.path(ppr.merged.bed.dir,"ppr.merged.bed"),header=F,colClasses = c("character"))
    alldata=data.frame(peakNames,as.data.frame(res.early$res.IL2.FP))
    # Write diffMat to disk
    write.csv(alldata,file = file.path(output.dir,"DE_IL2_FP.csv"))
    
    
    res_density=results(dds, name="Density_HD_vs_LD")
    png("maplot.density.png")
    plotMA(res.IL2.FP,ylim=c(-1,1))
    dev.off()
    
    #Treatment effect in low density samples
    res_treatment_ld= results(dds, contrast=c("Condition","IL2","PBS"))
    png("maplot.treatment.ld.png")
    plotMA(res_treatment_ld,ylim=c(-1,1))
    dev.off()
    
    #Treatment effect for high density samples
    res_treatment_hd= results(dds, list( c("Condition_treated_vs_control","DensityHD.Conditiontreated") ))
    png("maplot.treatment.hd.png")
    plotMA(res_treatment_hd,ylim=c(-1,1))
    dev.off()
    
    #THIS IS AN EXAMPLE OF A SINGLE FACTOR ANALYSIS -- SAMPLE IS THE SINGLE INDEPENDENT VARIABLE.
    
    (condition <- factor(c(rep("FP", 2), rep("PBS", 2),rep("IL2", 2))))
    
    coldata <- data.frame(row.names=colnames(countData.2), condition=colData.1$Condition)
    ddsIndividualFullCountTable<-DESeqDataSetFromMatrix(
        countData=countData.2,
        colData=coldata,
        design= ~condition)
    
    dds <- DESeq(ddsIndividualFullCountTable)
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    
    install.packages("locfit")
    rld <- rlog(dds,fitType='local')
    head(assay(rld))
    hist(assay(rld))
    library(RColorBrewer)
    (mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])
    
    sampleDists <- as.matrix(dist(t(assay(rld))))
    library(gplots)
    png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
    heatmap.2(as.matrix(sampleDists), key=F, trace="none",
              col=colorpanel(100, "black", "white"),
              ColSideColors=mycols[condition], RowSideColors=mycols[condition],
              margin=c(10, 10), main="At 16 hours")
    dev.off()
    
    rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
        require(genefilter)
        require(calibrate)
        require(RColorBrewer)
        rv = rowVars(assay(rld))
        select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        pca = prcomp(t(assay(rld)[select, ]))
        fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
        if (is.null(colors)) {
            if (nlevels(fac) >= 3) {
                colors = brewer.pal(nlevels(fac), "Paired")
            }   else {
                colors = c("black", "red")
            }
        }
        pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
        pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
        pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
        pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
        plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
        with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
        legend(legendpos, legend=levels(fac), col=colors, pch=20)
        #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
        #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
        #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
    }
    png("qc-pca.png", 1000, 1000, pointsize=20)
    rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
    dev.off()
    
    ddsIndividual<-DESeq(ddsIndividualFullCountTable,parallel=parallelFlag,betaPrior=FALSE)
    #create matrices to store the results from the PER-SAMPLE univariate comparison
    #get the individual comparisons
    conditionsToCompare <- data.frame(
        matrix(
            c('earlyG1.HD.treated', 'earlyG1.HD.controls',
              'earlyG1.LD.treated', 'earlyG1.LD.controls',
              'lateG1.HD.treated','lateG1.HD.controls',
              'lateG1.LD.treated','lateG1.LD.controls',
              'SG2M.HD.treated','SG2M.HD.controls',
              'SG2M.LD.treated','SG2M.LD.controls',
              'lateG1.HD.treated','earlyG1.HD.treated',
              'lateG1.LD.treated','earlyG1.LD.treated',
              'SG2M.HD.treated','lateG1.HD.treated',
              'SG2M.LD.treated','lateG1.LD.treated',
              'SG2M.HD.treated','earlyG1.HD.treated',
              'SG2M.LD.treated','earlyG1.LD.treated'),
            ncol = 2,
            byrow = TRUE),
        stringsAsFactors = FALSE)
    
    numCols <- 15
    numRows <- nrow(countData)
    diffMat <- matrix(, ncol = numCols, nrow = numRows)
    confidenceMatt<-matrix(,ncol=numCols, nrow=numRows)
    diffMat<-matrix(,ncol=numCols,nrow=numRows)
    diffMat[,1]=(res_density$padj<=padjCutoff)*(abs(res_density$log2FoldChange) >= foldChangeCutoff)*(res_density$log2FoldChange)
    diffMat[,2]=(res_treatment_ld$padj <=padjCutoff)*(abs(res_treatment_ld$log2FoldChange)>=foldChangeCutoff)*(res_treatment_ld$log2FoldChange)
    diffMat[,3]=(res_treatment_hd$padj <=padjCutoff)*(abs(res_treatment_hd$log2FoldChange)>=foldChangeCutoff)*(res_treatment_hd$log2FoldChange)
    
    colNameEntries=c("HD_vs_LD","LD_treated_vs_control","HD_treated_vs_control")
    for (i in 1:12){
        res <- results(ddsIndividual,
                       contrast = c("Sample",
                                    conditionsToCompare[i, 1],
                                    conditionsToCompare[i, 2]),
                       parallel = parallelFlag)
        colNameEntries=append(colNameEntries,paste(conditionsToCompare[i,1],"_vs_",conditionsToCompare[i,2],sep=""))
        #plot the MAplot
        #png(paste(conditionsToCompare[i,1],"_vs_",conditionsToCompare[i,2],".png",sep=""))
        #plotMA(res,ylim=c(-1,1))
        #dev.off()
        diffMat[, 3+i] <- (res$padj <= padjCutoff)*(abs(res$log2FoldChange) >=foldChangeCutoff)*(res$log2FoldChange)
    }
    # Replace NAs
    # Remember to take abs so that +1s and -1s don't cancel each other out
    diffMat[is.na(diffMat)] <- 0
    colnames(diffMat) <- colNameEntries
    
    #IF YOU WANT TO QC THE REPLICATES, YOU CAN PLOT RLOG-TRANSFORMED REPLICATE PAIRS AND CHECK CORRELATION
    #log-transform and plot pairs of replicates
    rld<-rlog(ddsIndividual)
    for(i in seq(1,24,2)){
        png(paste(as.character(i),"_",as.character(i+1),".png",sep=""))
        plot( assay(rld)[, i:i+1], col="#00000020", pch=20, cex=0.3, xlab=colData[i,"Sample"],ylab=colData[i+1,"Sample"])
        dev.off()
    }
}

# input.count.file <- "~/Dropbox/Alejandro_AtacSeq_out_2/count.txt"
# output.dir <- "~/Dropbox/Alejandro_AtacSeq_out_2"
# get.DE.pvalue.distribution(input.count.file,output.dir)
#
get.DE.pvalue.distribution <- function(input.count.file,output.dir) {
    
    # Metadata
    if (!dir.exists(dirname(output.dir)))
    {
        dir.create(dirname(output.dir), recursive = TRUE)
    }
    
    padjCutoff <- 0.01/20
    foldChangeCutoff <- 0.10
    register(MulticoreParam(32))
    parallelFlag <- TRUE
    
    countData=data.frame(read.table(input.count.file,header=T,sep='\t'))
    
    counts.metadata <- cbind.data.frame(Sample = colnames(countData),Condition=CutStringByNFromEnd(colnames(countData),15),Stage=trimGene3(colnames(countData)),stringsAsFactors=F)
    
    counts.metadata$Condition <- trimGene4(counts.metadata$Condition)
    
    counts.metadata[which(as.character(counts.metadata$Condition) == "IL"),]$Condition <- "IL2"
    
    counts.metadata.1 <- cbind.data.frame(counts.metadata,ID=paste0(counts.metadata$Condition,".at.",str_sub(counts.metadata$Sample,str_locate(counts.metadata$Sample,"at")[,2]+2,str_length(counts.metadata$Sample))),stringsAsFactors=F)
    
    counts.metadata.1[which(counts.metadata.1$Stage == 0.75),]$Stage <- "early"
    counts.metadata.1[which(counts.metadata.1$Stage == 16),]$Stage <- "late"
    
    counts.metadata.2 <- cbind.data.frame(Sample=counts.metadata.1$ID,Condition=counts.metadata.1$Condition,Stage=counts.metadata.1$Stage,stringsAsFactors=F)
    colnames(countData) <- counts.metadata.2$Sample
    any(is.na(countData))
    
    counts.metadata.2$group <- factor(paste0(counts.metadata.2$Condition,"_",counts.metadata.2$Stage))
    
    ddsFullCountTable<-DESeqDataSetFromMatrix(
        countData=countData,
        colData=counts.metadata.2,
        design= ~ group)
    
    dds<-DESeq(ddsFullCountTable,parallel=parallelFlag,betaPrior=FALSE)
    resultsNames(dds)
    
    res.late.vs.early.at.FP=results(dds, contrast = c("group","FP_late","FP_early"))
    res.late.vs.early.at.IL2=results(dds, contrast = c("group","IL2_late","IL2_early"))
    res.late.vs.early.at.PBS=results(dds, contrast = c("group","PBS_late","IL2_early"))
    
    res.FP.vs.PBS.at.early=results(dds, contrast = c("group","FP_early","PBS_early"))
    res.IL2.vs.PBS.at.early=results(dds, contrast = c("group","IL2_early","PBS_early"))
    res.FP.vs.IL2.at.early=results(dds, contrast = c("group","FP_early","IL2_early"))
    
    res.FP.vs.PBS.at.late=results(dds, contrast = c("group","FP_late","PBS_late"))
    res.IL2.vs.PBS.at.late=results(dds, contrast = c("group","IL2_late","PBS_late"))
    res.FP.vs.IL2.at.late=results(dds, contrast = c("group","FP_late","IL2_late"))
    
    res.FP_early.vs.IL2_late=results(dds, contrast = c("group","FP_early","IL2_late"))
    res.FP_late.vs.IL2_early=results(dds, contrast = c("group","FP_late","IL2_early"))
    res.FP_early.vs.PBS_late=results(dds, contrast = c("group","FP_early","PBS_late"))
    
    res.FP_late.vs.PBS_early=results(dds, contrast = c("group","FP_late","PBS_early"))
    res.IL2_early.vs.PBS_late=results(dds, contrast = c("group","IL2_early","PBS_late"))
    res.IL2_late.vs.PBS_early=results(dds, contrast = c("group","IL2_late","PBS_early"))
    
    ddsFullCountTable2<-DESeqDataSetFromMatrix(
        countData=countData,
        colData=counts.metadata.2,
        design= ~ Condition + Stage)
    
    dds2<-DESeq(ddsFullCountTable2,parallel=parallelFlag,betaPrior=FALSE)
    resultsNames(dds2)
    
    res.late.vs.early = results(dds2, contrast = c("Stage","late","early"))
    
    res.FP.vs.PBS = results(dds2, contrast = c("Condition","FP","PBS"))
    res.IL2.vs.PBS = results(dds2, contrast = c("Condition","IL2","PBS"))
    res.FP.vs.IL2 = results(dds2, contrast = c("Condition","FP","IL2"))
    
    pdf(file.path(output.dir,"DE_pvalue.pdf"))
    par(mfrow = c(3,2))  # 3 rows and 2 columns
    
    hist(res.late.vs.early$pvalue)
    hist(res.FP.vs.PBS$pvalue)
    hist(res.IL2.vs.PBS$pvalue)
    hist(res.FP.vs.IL2$pvalue)
    
    hist(res.late.vs.early.at.FP$pvalue)
    hist(res.late.vs.early.at.IL2$pvalue)
    hist(res.late.vs.early.at.PBS$pvalue)
    
    hist(res.FP.vs.PBS.at.early$pvalue)
    hist(res.IL2.vs.PBS.at.early$pvalue)
    hist(res.FP.vs.IL2.at.early$pvalue)
    
    hist(res.FP.vs.PBS.at.late$pvalue)
    hist(res.IL2.vs.PBS.at.late$pvalue)
    hist(res.FP.vs.IL2.at.late$pvalue)
    
    hist(res.FP_early.vs.IL2_late$pvalue)
    hist(res.FP_late.vs.IL2_early$pvalue)
    hist(res.FP_early.vs.PBS_late$pvalue)
    
    hist(res.FP_late.vs.PBS_early$pvalue)
    hist(res.IL2_early.vs.PBS_late$pvalue)
    hist(res.IL2_late.vs.PBS_early$pvalue)
    
    dev.off()
    
}

# input.atac.dir <- "/Users/axy148/.CMVolumes/AiminDropbox/Alejandro_AtacSeq_uploaded_2"
# output.dir <- "~/Dropbox/Alejandro_AtacSeq_out_3"
# ppr.merged.bed.dir <- "/Users/axy148/Dropbox/Alejandro_AtacSeq_out"
#
# ppr.merged.bed.file <- "~/Dropbox/Alejandro_AtacSeq_out/union.ppr.merged.bed"
#
# getCount2(input.atac.dir,ppr.merged.bed.file,output.dir)
#
# Example2
#
# input.atac.dir <- "/Users/axy148/.CMVolumes/AiminDropbox/Alejandro_AtacSeq_uploaded_2"
# output.dir <- "~/Dropbox/Alejandro_AtacSeq_out_4"
# ppr.merged.bed.file <- "~/Dropbox/Alejandro_AtacSeq_out/union.ppr.merged.bed"
# getCount2(input.atac.dir,ppr.merged.bed.file,output.dir)

getCount2 <- function(input.atac.dir,ppr.merged.bed.file,output.dir) {
    
    if (!dir.exists(dirname(output.dir)))
    {
        dir.create(dirname(output.dir), recursive = TRUE)
    }
    
    #file.1 <- list.files(input.atac.dir,pattern="*tn5_pooled.tagAlign.gz$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    file.1 <- list.files(input.atac.dir,pattern="*tn5.tagAlign.gz$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    file.2 <- file.1[-grep("pseudo_reps",file.1)]
    
    #file.2 <- file.1[grep("pooled_pseudo_reps",file.1)]
    
    file.3 <- cbind.data.frame(ID = basename((dirname(dirname(dirname(file.2))))),rep= basename(dirname(file.2)),fileName=basename(file.2),file.2,stringsAsFactors=FALSE)
    
    file.4 <- file.3[-which(file.3$ID %in% c("IL-2vsPBS")),]
    
    file.5 <- cbind.data.frame(ID2 = paste0(file.4$ID,"_",file.4$rep),file.4)
    
    file.id <- unique(as.character(file.5$ID2))
    
    cmd.l <- lapply(file.id, function(u,file.5){
        
        
        x <- file.5[which(as.character(file.5$ID2) == u),]
        
        cat("\n")
        cat(as.character(x$ID2),"\n")
        
        if(as.character(x$ID) == "IL-2-at-16"){
            
            cmd0 <- paste("gzip -d",x$file.2,"-c >",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),sep = " ")
            
            cmd1 <- paste("python ~/Dropbox/Aimin_project/AtacSeq/inst/deseq2_wrappers/workflow_for_raw_counts/shift_atac_tagalign.py",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),file.path(output.dir,paste0(x$ID2,"slopped")),sep = " ")
            
            cmd2 <- paste("bedtools coverage -counts -a",ppr.merged.bed.file,"-b",file.path(output.dir,paste0(x$ID2,"slopped")),"| cut -f4 >",file.path(output.dir,x$ID2),sep = " ")
            
            cmd <- paste(cmd0,"wait",cmd1,"wait",cmd2,sep = ";")
            
        }else
        {
            cmd0 <- paste("gzip -d",x$file.2,"-c >",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),sep = " ")
            cmd1 <- paste("python ~/Dropbox/Aimin_project/AtacSeq/inst/deseq2_wrappers/workflow_for_raw_counts/shift_atac_tagalign.py",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),file.path(output.dir,paste0(x$ID2,"slopped")),sep = " ")
            cmd2 <- paste("bedtools coverage -counts -a",ppr.merged.bed.file,"-b",file.path(output.dir,paste0(x$ID2,"slopped")),"| cut -f4 >",file.path(output.dir,x$ID2),sep = " ")
            cmd <- paste(cmd0,"wait",cmd1,"wait",cmd2,sep = ";")
        }
        print(cmd)
        system(cmd)
    },file.5)
}

# cmd <- "ls -lrth ~/Dropbox/Alejandro_AtacSeq_out_3/*rep? | awk '{print $9}'"
# output.dir <- "~/Dropbox/Alejandro_AtacSeq_out_3"
# makeCountTable3(cmd,output.dir)
#
makeCountTable3 <- function(cmd,output.dir) {
    
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }
    
    file.count <- system(cmd,intern = T)
    system(paste("echo",paste(basename(file.count),collapse = "\\\t"),">",file.path(output.dir,"count.txt")))
    cmd1 <- paste("paste",paste(file.count,collapse = " "), ">>",file.path(output.dir,"count.txt"))
    system(cmd1)
}

# input.count.file <- "~/Dropbox/Alejandro_AtacSeq_out_3/count.txt"
# output.dir <- "~/Dropbox/Alejandro_AtacSeq_out_3"
# res <- AtacSeq:::get.DE.pvalue.distribution.1(input.count.file,output.dir)
#
get.DE.pvalue.distribution.1 <- function(input.count.file,output.dir) {
    
    # Metadata
    if (!dir.exists(dirname(output.dir)))
    {
        dir.create(dirname(output.dir), recursive = TRUE)
    }
    
    padjCutoff <- 0.01/20
    foldChangeCutoff <- 0.10
    register(MulticoreParam(32))
    parallelFlag <- TRUE
    
    countData=data.frame(read.table(input.count.file,header=T,sep='\t'))
    
    counts.metadata <- cbind.data.frame(Sample = colnames(countData),Condition=colnames(countData),Stage=trimGene3(colnames(countData)),stringsAsFactors=F)
    
    counts.metadata$Condition <- trimGene4(counts.metadata$Condition)
    
    counts.metadata[which(as.character(counts.metadata$Condition) == "IL"),]$Condition <- "IL2"
    
    counts.metadata.1 <- cbind.data.frame(counts.metadata,ID=paste0(counts.metadata$Condition,".at.",str_sub(counts.metadata$Sample,str_locate(counts.metadata$Sample,"at")[,2]+2,str_length(counts.metadata$Sample))),stringsAsFactors=F)
    
    counts.metadata.1[which(counts.metadata.1$Stage == 0.75),]$Stage <- "early"
    counts.metadata.1[which(counts.metadata.1$Stage == 16),]$Stage <- "late"
    
    counts.metadata.2 <- cbind.data.frame(Sample=counts.metadata.1$ID,Condition=counts.metadata.1$Condition,Stage=counts.metadata.1$Stage,stringsAsFactors=F)
    colnames(countData) <- counts.metadata.2$Sample
    any(is.na(countData))
    
    counts.metadata.2$group <- factor(paste0(counts.metadata.2$Condition,"_",counts.metadata.2$Stage))
    
    ddsFullCountTable<-DESeqDataSetFromMatrix(
        countData=countData,
        colData=counts.metadata.2,
        design= ~ group)
    
    dds<-DESeq(ddsFullCountTable,parallel=parallelFlag,betaPrior=FALSE)
    resultsNames(dds)
    
    res.late.vs.early.at.FP=results(dds, contrast = c("group","FP_late","FP_early"))
    res.late.vs.early.at.IL2=results(dds, contrast = c("group","IL2_late","IL2_early"))
    res.late.vs.early.at.PBS=results(dds, contrast = c("group","PBS_late","IL2_early"))
    
    res.FP.vs.PBS.at.early=results(dds, contrast = c("group","FP_early","PBS_early"))
    res.IL2.vs.PBS.at.early=results(dds, contrast = c("group","IL2_early","PBS_early"))
    res.FP.vs.IL2.at.early=results(dds, contrast = c("group","FP_early","IL2_early"))
    
    res.FP.vs.PBS.at.late=results(dds, contrast = c("group","FP_late","PBS_late"))
    res.IL2.vs.PBS.at.late=results(dds, contrast = c("group","IL2_late","PBS_late"))
    res.FP.vs.IL2.at.late=results(dds, contrast = c("group","FP_late","IL2_late"))
    
    res.FP_early.vs.IL2_late=results(dds, contrast = c("group","FP_early","IL2_late"))
    res.FP_late.vs.IL2_early=results(dds, contrast = c("group","FP_late","IL2_early"))
    res.FP_early.vs.PBS_late=results(dds, contrast = c("group","FP_early","PBS_late"))
    
    res.FP_late.vs.PBS_early=results(dds, contrast = c("group","FP_late","PBS_early"))
    res.IL2_early.vs.PBS_late=results(dds, contrast = c("group","IL2_early","PBS_late"))
    res.IL2_late.vs.PBS_early=results(dds, contrast = c("group","IL2_late","PBS_early"))
    
    ddsFullCountTable2<-DESeqDataSetFromMatrix(
        countData=countData,
        colData=counts.metadata.2,
        design= ~ Condition + Stage)
    
    dds2<-DESeq(ddsFullCountTable2,parallel=parallelFlag,betaPrior=FALSE)
    resultsNames(dds2)
    
    res.late.vs.early = results(dds2, contrast = c("Stage","late","early"))
    
    res.FP.vs.PBS = results(dds2, contrast = c("Condition","FP","PBS"))
    res.IL2.vs.PBS = results(dds2, contrast = c("Condition","IL2","PBS"))
    res.FP.vs.IL2 = results(dds2, contrast = c("Condition","FP","IL2"))
    
    pdf(file.path(output.dir,"DE_pvalue.pdf"))
    par(mfrow = c(3,2))  # 3 rows and 2 columns
    
    hist(res.late.vs.early$pvalue)
    hist(res.FP.vs.PBS$pvalue)
    hist(res.IL2.vs.PBS$pvalue)
    hist(res.FP.vs.IL2$pvalue)
    
    hist(res.late.vs.early.at.FP$pvalue)
    hist(res.late.vs.early.at.IL2$pvalue)
    hist(res.late.vs.early.at.PBS$pvalue)
    
    hist(res.FP.vs.PBS.at.early$pvalue)
    hist(res.IL2.vs.PBS.at.early$pvalue)
    hist(res.FP.vs.IL2.at.early$pvalue)
    
    hist(res.FP.vs.PBS.at.late$pvalue)
    hist(res.IL2.vs.PBS.at.late$pvalue)
    hist(res.FP.vs.IL2.at.late$pvalue)
    
    hist(res.FP_early.vs.IL2_late$pvalue)
    hist(res.FP_late.vs.IL2_early$pvalue)
    hist(res.FP_early.vs.PBS_late$pvalue)
    
    hist(res.FP_late.vs.PBS_early$pvalue)
    hist(res.IL2_early.vs.PBS_late$pvalue)
    hist(res.IL2_late.vs.PBS_early$pvalue)
    
    dev.off()
    
    re <- list(
        dds = dds,
        dds2 = dds2,
        res.late.vs.early = res.late.vs.early,
        res.FP.vs.PBS =  res.FP.vs.PBS,
        res.IL2.vs.PBS = res.IL2.vs.PBS,
        res.FP.vs.IL2 = res.FP.vs.IL2,
        res.late.vs.early.at.FP = res.late.vs.early.at.FP,
        res.late.vs.early.at.IL2 = res.late.vs.early.at.IL2,
        res.late.vs.early.at.PBS = res.late.vs.early.at.PBS,
        res.FP.vs.PBS.at.early = res.FP.vs.PBS.at.early,
        res.IL2.vs.PBS.at.early = res.IL2.vs.PBS.at.early,
        res.FP.vs.IL2.at.early = res.FP.vs.IL2.at.early,
        res.FP.vs.PBS.at.late =  res.FP.vs.PBS.at.late,
        res.IL2.vs.PBS.at.late = res.IL2.vs.PBS.at.late,
        res.FP.vs.IL2.at.late = res.FP.vs.IL2.at.late,
        res.FP_early.vs.IL2_late = res.FP_early.vs.IL2_late,
        res.FP_late.vs.IL2_early = res.FP_late.vs.IL2_early,
        res.FP_early.vs.PBS_late = res.FP_early.vs.PBS_late,
        res.FP_late.vs.PBS_early =  res.FP_late.vs.PBS_early,
        res.IL2_early.vs.PBS_late = res.IL2_early.vs.PBS_late,
        res.IL2_late.vs.PBS_early = res.IL2_late.vs.PBS_early
    )
    
    re
    
}

#' dir.name="~/Dropbox/Alejandro_AtacSeq_out"
#' input.file.pattern="*.merged.bed"
#' out.dir.name="~/Dropbox/Alejandro_AtacSeq_out_3"
#'
#' txdb="mm10"
#' DD=5000
#'
#' AtacSeq:::anotatePeak(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000, AP=c("Promoter","Intron"))
#'
#' res.promoter <- AtacSeq:::anotatePeak(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000,AP=c("Promoter"))
#'
#' AtacSeq:::anotatePeak(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000,AP=c("Intron"))
#'
#' dir.name="~/Dropbox/Alejandro_AtacSeq_out"
#' input.file.pattern="ppr.merged.bed"
#' out.dir.name="~/Dropbox/Alejandro_AtacSeq_out_3"
#'
#'
#' AtacSeq:::anotatePeak(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000)

anotatePeak <- function(dir.name,input.file.pattern,out.dir.name,txdb=c("hg19","hg38","mm10"),DD,distanceToTSS_cutoff=NULL,assignGenomicAnnotation=TRUE,AP=c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic")) {
    
    re <- list.files(dir.name,pattern=input.file.pattern,all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    re.peaks.only.bed.2 <- re[-c(1,8)]
    
    re.peaks.only.bed.2 <- re[-c(1,8)]
    
    txdb<-match.arg(txdb)
    
    APpath <- paste(AP,collapse = "_")
    
    temp3=file.path(out.dir.name,APpath)
    
    if(!dir.exists(temp3)){dir.create(temp3,recursive = TRUE)}
    
    d=DD
    
    peaks.anno.list <- lapply(re.peaks.only.bed.2,function(u,d){
        
        peaks=readPeakFile(u,as="data.frame")
        
        print(head(peaks))
        
        if(dim(peaks)[1]>0){
            
            switch (txdb,
                    hg38 = {
                        cat("Use hg38\n")
                        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
                        peakAnno <- annotatePeak(u, tssRegion=c(-d, d),
                                                 TxDb=txdb,assignGenomicAnnotation=assignGenomicAnnotation,
                                                 genomicAnnotationPriority=AP,annoDb="org.Hs.eg.db")
                    },
                    h19 = {
                        cat("Use hg19\n")
                        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
                        peakAnno <- annotatePeak(u, tssRegion=c(-d, d),
                                                 TxDb=txdb,assignGenomicAnnotation=assignGenomicAnnotation,
                                                 genomicAnnotationPriority=AP,annoDb="org.Hs.eg.db")
                    },
                    mm10 = {
                        cat("Use mm10\n")
                        txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
                        peakAnno <- annotatePeak(u, tssRegion=c(-d, d),
                                                 TxDb=txdb,assignGenomicAnnotation=assignGenomicAnnotation,
                                                 genomicAnnotationPriority=AP,annoDb="org.Mm.eg.db")
                    }
            )
            
            dropAnnoM <- function (csAnno, distanceToTSS_cutoff)
            {
                idx <- which(abs(mcols(csAnno@anno)[["distanceToTSS"]]) <
                                 distanceToTSS_cutoff)
                csAnno@anno <- csAnno@anno[idx]
                csAnno@peakNum <- length(idx)
                if (csAnno@hasGenomicAnnotation) {
                    csAnno@annoStat <- ChIPseeker:::getGenomicAnnoStat(csAnno@anno)
                    csAnno@detailGenomicAnnotation = csAnno@detailGenomicAnnotation[idx,
                                                                                    ]
                }
                csAnno
            }
            
            if(!is.null(distanceToTSS_cutoff)){
                peakAnno <- dropAnnoM(peakAnno,distanceToTSS_cutoff = distanceToTSS_cutoff)
            }else
            {
                peakAnno <- peakAnno
            }
            
            x_name=tools::file_path_sans_ext(basename(u))
            
            cat(x_name,"\n")
            png(file.path(temp3,paste0(x_name,"_",d,".png")))
            plotAnnoPie(peakAnno)
            dev.off()
            
            peaks.anno=as.data.frame(peakAnno)
            
            print(head(peaks.anno))
            
            print(colnames(peaks.anno))
            write.table(peaks.anno,file=file.path(temp3,paste0(x_name,"_",d,".xls")),
                        row.names = FALSE,quote=FALSE,sep="\t")
            
            peaks.anno
        }
        
    },d)
    
}

# ppr.merged.bed.dir <- "~/Dropbox/Alejandro_AtacSeq_out"
# out.dir.name <- "~/Dropbox/Alejandro_AtacSeq_out_3"
# AtacSeq:::makeHeatMap(ppr.merged.bed.dir,res,out.dir.name,0.05,0.58)
#
#
makeHeatMap <- function(ppr.merged.bed.dir,res,out.dir.name,padjust_thresh,lfc_thresh) {
    
    res.1 <- res[10:15]
    
    peakNames <- read.table(file.path(ppr.merged.bed.dir,"ppr.merged.bed"),header=F,colClasses = c("character"))
    
    cmd.l <- lapply(1:length(res.1), function(u,res.1,res,out.dir.name){
        
        x <- names(res.1)[u]
        x.1 <- str_sub(x,str_locate(x,"\\.")[2]+1,str_locate(x,"vs")[1]-2)
        x.2 <- str_sub(x,str_locate(x,"vs")[2]+2,str_locate(x,"at")[1]-2)
        x.3 <- str_sub(x,str_locate(x,"at")[2]+2,str_length(x))
        
        x.name <- paste(x.1,x.2,x.3,sep="-")
        
        if(x.3 == "early"){
            y1 <- paste(x.1,"0.75",sep=".at.")
            y2 <- paste(x.2,"0.75",sep=".at.")
            
        }else
        {
            y1 <- paste(x.1,"16",sep=".at.")
            y2 <- paste(x.2,"16",sep=".at.")
        }
        
        alldata=data.frame(peakNames,as.data.frame(res.1[[u]]))
        
        index <- c(grep(y1,colnames(counts(res$dds,normalized = T))),grep(y2,colnames(counts(res$dds,normalized = T))))
        
        res.with.counts <- cbind.data.frame(peakName= seq(1,dim(alldata)[1]),alldata,counts(res$dds,normalized = T)[,index],stringsAsFactors=F)
        
        write.csv(res.with.counts,file = file.path(out.dir.name,"DE",paste0("DE_",x.name,".csv")))
        
        res.with.counts.1 <- res.with.counts[which(res.with.counts$padj<padjust_thresh),]
        
        res.with.counts.2 <- res.with.counts.1[which(abs(res.with.counts.1$log2FoldChange) > lfc_thresh),]
        
        res.with.counts.3 <- res.with.counts.2[order(res.with.counts.2$padj),]
        
        res.with.counts.4 <- as.matrix(res.with.counts.3[,-c(1:10)])
        
        n = dim(res.with.counts.4)[1]
        
        cat(x.name," ",n," ")
        
        if(n > 30){
            s.n = 30
            cat("OK","\n")
            pdf(file.path(out.dir.name,"DE",paste0("HeatMap_",x.name,".pdf")))
            heatmap.2(res.with.counts.4[1:s.n,],margins = c(15,5),dendrogram="column",Colv = T, Rowv=F,col=rev(heat.colors(16)),cexRow=0.5,cexCol=1,sepcolor="white", trace="none")
            dev.off()
        }
        else{
            cat("less than 30 to satisfy the selection thresh","\n")
        }
        
    },res.1,res,out.dir.name)
    
}

#  dir.name="~/Dropbox/Alejandro_AtacSeq_out"
#' input.file.pattern="ppr.merged.bed"
#' out.dir.name="~/Dropbox/Alejandro_AtacSeq_out_3"
#' txdb="mm10"
#' DD=5000
#'
#' anno.res <- AtacSeq:::AnntationUsingChipSeeker(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000)

AnntationUsingChipSeeker <- function(dir.name,input.file.pattern,out.dir.name,txdb=c("hg19","hg38","mm10"),DD,distanceToTSS_cutoff=5000,assignGenomicAnnotation=TRUE,AP=c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic")) {
    
    file.1 <- list.files(dir.name,pattern=input.file.pattern,all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    re.peaks.only.bed.2 <- file.1
    
    # txdb<-match.arg(txdb)
    
    txdb<-match.arg(txdb,c("hg19","hg38","mm10"),several.ok=TRUE)
    
    
    APpath <- paste(AP,collapse = "_")
    
    temp3=file.path(out.dir.name,"Annotation",APpath)
    
    if(!dir.exists(temp3)){dir.create(temp3,recursive = TRUE)}
    
    d=DD
    res <- lapply(1:length(re.peaks.only.bed.2),function(u,re.peaks.only.bed.2,d){
        
        peaks=readPeakFile(re.peaks.only.bed.2[[u]],as="data.frame")
        
        print(head(peaks))
        
        switch (txdb,
                hg38 = {
                    cat("Use hg38\n")
                    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
                    peakAnno <- annotatePeak(re.peaks.only.bed.2[[u]], tssRegion=c(-d, d),
                                             TxDb=txdb,assignGenomicAnnotation=assignGenomicAnnotation,
                                             genomicAnnotationPriority=AP,annoDb="org.Hs.eg.db")
                },
                hg19 = {
                    cat("Use hg19\n")
                    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
                    peakAnno <- annotatePeak(re.peaks.only.bed.2[[u]], tssRegion=c(-d, d),
                                             TxDb=txdb,assignGenomicAnnotation=assignGenomicAnnotation,
                                             genomicAnnotationPriority=AP,annoDb="org.Hs.eg.db")
                },
                mm10 = {
                    cat("Use mm10\n")
                    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
                    peakAnno <- annotatePeak(re.peaks.only.bed.2[[u]], tssRegion=c(-d, d),
                                             TxDb=txdb,assignGenomicAnnotation=assignGenomicAnnotation,
                                             genomicAnnotationPriority=AP,annoDb="org.Mm.eg.db")
                }
                
        )
        
        dropAnnoM <- function (csAnno, distanceToTSS_cutoff)
        {
            idx <- which(abs(mcols(csAnno@anno)[["distanceToTSS"]]) <
                             distanceToTSS_cutoff)
            csAnno@anno <- csAnno@anno[idx]
            csAnno@peakNum <- length(idx)
            if (csAnno@hasGenomicAnnotation) {
                csAnno@annoStat <- ChIPseeker:::getGenomicAnnoStat(csAnno@anno)
                csAnno@detailGenomicAnnotation = csAnno@detailGenomicAnnotation[idx,                                                                        ]
            }
            csAnno
        }
        
        peakAnno <- dropAnnoM(peakAnno,distanceToTSS_cutoff = distanceToTSS_cutoff)
        
        x_name=basename(re.peaks.only.bed.2[[u]])
        cat(x_name)
        
        png(file.path(temp3,paste0(x_name,"_",d,"_around_tss_annotation_pie.png")))
        plotAnnoPie(peakAnno)
        dev.off()
        
        peaks.anno=as.data.frame(peakAnno)
        
        print(head(peaks.anno))
        
        print(colnames(peaks.anno))
        
        peaks.anno.1 <- data.frame(peakName=rownames(peaks.anno),peaks.anno)
        
        write.table(peaks.anno.1,file=file.path(temp3,paste0(x_name,"_",d,"_around_tss_annotation_4_only_mapped_peaks.xls")),
                    row.names = FALSE,quote=FALSE,sep="\t")
        
        
        unmapped.peaks<-peaks[-which(paste0(peaks[,1],"_",peaks[,2],"_",peaks[,3]) %in% paste0(peaks.anno[,1],"_",peaks.anno[,2],"_",peaks.anno[,3])),]
        
        cat(dim(peaks)[1]," ",dim(peaks.anno)[1]," ",dim(unmapped.peaks)[1],"\n")
        
        if(dim(unmapped.peaks)[1]!=0){
            
            colnames(unmapped.peaks)=colnames(peaks.anno)[1:3]
            
            unmapped.peaks.3<-smartbind(peaks.anno,unmapped.peaks)
            
            unmapped.peaks.4<-unmapped.peaks.3[order(unmapped.peaks.3[,1],unmapped.peaks.3[,2]),]
            
            unmapped.peaks.4[which(is.na(unmapped.peaks.4$annotation)),6] <- "NotAnno"
            
            row.names(unmapped.peaks.4) <- seq(1,dim(unmapped.peaks.4)[1])
            
            unmapped.peaks.5 <- data.frame(peakName=rownames(unmapped.peaks.4),unmapped.peaks.4)
            
            write.table(unmapped.peaks.5,file=file.path(temp3,paste0(x_name,"_",d,"_around_tss_annotation_4_all_peaks.xls")),row.names = FALSE,quote=FALSE,sep="\t")
        }
        
        f.only <- file.path(temp3,paste0(x_name,"_",d,"_Pie_peaksMappedOnly.pdf"))
        f.all <- file.path(temp3,paste0(x_name,"_",d,"_Pie_peaksAll.pdf"))
        
        makePieChart(peaks.anno,f.only)
        makePieChart(unmapped.peaks.4,f.all)
        
        re <- list(peaksMappedOnly=peaks.anno,peaksAll=unmapped.peaks.4)
        re
    },re.peaks.only.bed.2,d)
    res
}

makePieChart <-function(anno.res,fileN) {
    
    pdf(fileN,width=10, height=7)
    slices <- table(anno.res$annotation)
    lbls <- names(table(anno.res$annotation))
    pct <- round(slices/sum(slices)*100)
    lbls <- paste(lbls, pct) # add percents to labels
    lbls <- paste(lbls,"%",sep="") # ad % to labels
    pie(slices,labels = lbls, col=rainbow(length(lbls)),main="Annotation")
    dev.off()
    
}

makePeakList <- function() {
    
    ppr.merged.bed <- system("ls -lrth ~/Dropbox/Alejandro_AtacSeq_out/*.merged.bed | cut -d' ' -f12",intern = T)
    
    inst.peak.list <- ppr.merged.bed[1:6]
    
    ppr.merged.bed.l <- lapply(1:length(inst.peak.list),function(u,inst.peak.list){
        
        peak.gr <- readPeakFile(inst.peak.list[u], as="GRanges")
        peak.gr
    },inst.peak.list)
    
    ppr.merged.bed.name<- unlist(lapply(1:length(inst.peak.list),function(u,inst.peak.list){
        
        peak.name <- str_sub(basename(inst.peak.list[u]),1,str_locate(basename(inst.peak.list[u]),"ppr")[1]-2)
        
        peak.name
    },inst.peak.list))
    
    names(ppr.merged.bed.l) <- ppr.merged.bed.name
    
}

Compare2Peaks <- function(ppr.merged.bed.l,i,j) {
    ol <- findOverlapsOfPeaks(ppr.merged.bed.l[[i]], ppr.merged.bed.l[[j]],connectedPeaks="merge")
    
    x <- str_sub(names(ppr.merged.bed.l)[c(i,j)],1,str_locate(names(ppr.merged.bed.l)[c(i,j)],"-at")[,1]-1)
    
    y <- unique(str_sub(names(ppr.merged.bed.l)[c(i,j)],str_locate(names(ppr.merged.bed.l)[c(i,j)],"-at")[,2]+2,str_length(names(ppr.merged.bed.l)[c(i,j)])))
    pdf(file.path(output.dir,paste0(paste(x,collapse = "_"),"_at_",y,".pdf")))
    makeVennDiagram(ol,NameOfPeaks=x,totalTest=140938,fill = c("red","blue"),main=paste("Venn diagram of",paste(x,collapse = " and "),"ATAC-seq peaks at",y,"time point",sep = " "))
    dev.off()
}

# Compare2Peaks(ppr.merged.bed.l,1,5)
# Compare2Peaks(ppr.merged.bed.l,4,5)
# Compare2Peaks(ppr.merged.bed.l,1,4)
# Compare2Peaks(ppr.merged.bed.l,2,6)
# Compare2Peaks(ppr.merged.bed.l,3,6)
# Compare2Peaks(ppr.merged.bed.l,2,3)

getOverlapOf2Peaks <- function(ppr.merged.bed.l,i,j) {
    
    ol <- findOverlapsOfPeaks(ppr.merged.bed.l[[i]], ppr.merged.bed.l[[j]],connectedPeaks="merge")
    
    peak.bed <- as.data.frame(ol$mergedPeaks)
    
    x <- str_sub(names(ppr.merged.bed.l)[c(i,j)],1,str_locate(names(ppr.merged.bed.l)[c(i,j)],"-at")[,1]-1)
    
    y <- unique(str_sub(names(ppr.merged.bed.l)[c(i,j)],str_locate(names(ppr.merged.bed.l)[c(i,j)],"-at")[,2]+2,str_length(names(ppr.merged.bed.l)[c(i,j)])))
    
    file.name = file.path(output.dir,paste0(paste(x,collapse = "_"),"_at_",y,".bed"))
    write.table(peak.bed[1:3],row.names = F,col.names = F,quote =F, file= file.name,sep= "\t")
    
}

# getOverlapOf2Peaks(ppr.merged.bed.l,1,5)
# getOverlapOf2Peaks(ppr.merged.bed.l,4,5)
# getOverlapOf2Peaks(ppr.merged.bed.l,1,4)
# getOverlapOf2Peaks(ppr.merged.bed.l,2,6)
# getOverlapOf2Peaks(ppr.merged.bed.l,3,6)
# getOverlapOf2Peaks(ppr.merged.bed.l,2,3)


getUniqueOf2Peaks <- function(ppr.merged.bed.l,i,j) {
    
    # i=1
    #  j=5
    
    ol <- findOverlapsOfPeaks(ppr.merged.bed.l[[i]], ppr.merged.bed.l[[j]],connectedPeaks="merge")
    
    peak.bed <- as.data.frame(ol$uniquePeaks)
    
    x <- str_sub(names(ppr.merged.bed.l)[c(i,j)],1,str_locate(names(ppr.merged.bed.l)[c(i,j)],"-at")[,1]-1)
    
    y <- unique(str_sub(names(ppr.merged.bed.l)[c(i,j)],str_locate(names(ppr.merged.bed.l)[c(i,j)],"-at")[,2]+2,str_length(names(ppr.merged.bed.l)[c(i,j)])))
    
    file.name = file.path(output.dir,paste0(paste(x,collapse = "_"),"_un_",y,".bed"))
    write.table(peak.bed[1:3],row.names = F,col.names = F,quote =F, file= file.name,sep= "\t")
    
}

# getUniqueOf2Peaks(ppr.merged.bed.l,1,5)
# getUniqueOf2Peaks(ppr.merged.bed.l,4,5)
# getUniqueOf2Peaks(ppr.merged.bed.l,1,4)
# getUniqueOf2Peaks(ppr.merged.bed.l,2,6)
# getUniqueOf2Peaks(ppr.merged.bed.l,3,6)
# getUniqueOf2Peaks(ppr.merged.bed.l,2,3)


getAllOf2Peaks <- function(ppr.merged.bed.l,i,j,output.dir) {
    
    # i=1
    #  j=5
    
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }
    
    ol <- findOverlapsOfPeaks(ppr.merged.bed.l[[i]], ppr.merged.bed.l[[j]],connectedPeaks="merge")
    
    peak.uniq <- as.data.frame(ol$uniquePeaks)
    
    peak.uniq.1 <- cbind.data.frame(peak.uniq,peakNames=rownames(peak.uniq))
    
    peak.merge <- as.data.frame(ol$mergedPeaks)
    
    peak.all <- rbind.data.frame(peak.merge,peak.uniq.1)
    
    peak.bed <- peak.all[order(peak.all[,1],peak.all[,2]),]
    
    x <- str_sub(names(ppr.merged.bed.l)[c(i,j)],1,str_locate(names(ppr.merged.bed.l)[c(i,j)],"-at")[,1]-1)
    
    y <- unique(str_sub(names(ppr.merged.bed.l)[c(i,j)],str_locate(names(ppr.merged.bed.l)[c(i,j)],"-at")[,2]+2,str_length(names(ppr.merged.bed.l)[c(i,j)])))
    
    file.name = file.path(output.dir,paste0(paste(x,collapse = "_"),"_al_",y,".bed"))
    write.table(peak.bed[1:3],row.names = F,col.names = F,quote =F, file= file.name,sep= "\t")
    
}

# getAllOf2Peaks(ppr.merged.bed.l,1,5,"~/Dropbox/Alejandro_AtacSeq_out_All_Peaks")
# getAllOf2Peaks(ppr.merged.bed.l,4,5,"~/Dropbox/Alejandro_AtacSeq_out_All_Peaks")
# getAllOf2Peaks(ppr.merged.bed.l,1,4,"~/Dropbox/Alejandro_AtacSeq_out_All_Peaks")
# getAllOf2Peaks(ppr.merged.bed.l,2,6,"~/Dropbox/Alejandro_AtacSeq_out_All_Peaks")
# getAllOf2Peaks(ppr.merged.bed.l,3,6,"~/Dropbox/Alejandro_AtacSeq_out_All_Peaks")
# getAllOf2Peaks(ppr.merged.bed.l,2,3,"~/Dropbox/Alejandro_AtacSeq_out_All_Peaks")


# input.atac.dir <- "/Users/axy148/.CMVolumes/AiminDropbox/Alejandro_AtacSeq_uploaded_2"
# output.dir <- "~/Dropbox/Alejandro_AtacSeq_out_5"
# ppr.merged.bed.file <- "/Users/axy148/Dropbox/Alejandro_AtacSeq_out/FP_PBS_at_0.75.bed"
# getCount3(input.atac.dir,ppr.merged.bed.file,output.dir)

getCount3 <- function(input.atac.dir,ppr.merged.bed.file,output.dir) {
    
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }
    
    #file.1 <- list.files(input.atac.dir,pattern="*tn5_pooled.tagAlign.gz$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    file.1 <- list.files(input.atac.dir,pattern="*tn5.tagAlign.gz$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    file.2 <- file.1[-grep("pseudo_reps",file.1)]
    
    #file.2 <- file.1[grep("pooled_pseudo_reps",file.1)]
    
    file.3 <- cbind.data.frame(ID = basename((dirname(dirname(dirname(file.2))))),rep= basename(dirname(file.2)),fileName=basename(file.2),file.2,stringsAsFactors=FALSE)
    
    file.4 <- file.3[-which(file.3$ID %in% c("IL-2vsPBS")),]
    
    file.5 <- cbind.data.frame(ID2 = paste0(file.4$ID,"_",file.4$rep),file.4,stringsAsFactors=FALSE)
    
    file.5[grep("IL-2-at-16_rep1",file.5$ID2),]$ID2 <-  "IL2-at-16_rep1"
    file.5[grep("IL-2-at-16_rep2",file.5$ID2),]$ID2 <-  "IL2-at-16_rep2"
    file.5[grep("IL-2-at-16_rep3",file.5$ID2),]$ID2 <-  "IL2-at-16_rep3"
    
    file.5[grep("IL-2-at-16",file.5$ID),]$ID <-  "IL2-at-16"
    
    
    file.id <- unique(as.character(file.5$ID2))
    
    #ppr.merged.bed.file="/Users/axy148/Dropbox/Alejandro_AtacSeq_out_All_Peaks/FP_IL2_al_16.bed"
    
    XName <- basename(ppr.merged.bed.file)
    
    XName.1 <- str_sub(XName,1,str_locate_all(XName,"_")[[1]][1,1]-1)
    
    XName.2 <- str_sub(XName,str_locate_all(XName,"_")[[1]][1,2]+1,str_locate_all(XName,"_")[[1]][2,1]-1)
    
    XName.3 <- str_sub(XName,str_locate_all(XName,"_")[[1]][3,2]+1,str_locate_all(XName,".bed")[[1]][1,1]-1)
    
    p1 <- paste(XName.1,"at",XName.3,sep="-")
    p2 <- paste(XName.2,"at",XName.3,sep="-")
    
    file.id <- file.id[c(grep(p1,file.id),grep(p2,file.id))]
    #  output.dir <- "/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/outFPIL216-"
    
    cmd.l <- lapply(file.id, function(u,file.5){
        
        
        x <- file.5[which(as.character(file.5$ID2) == u),]
        
        cat("\n")
        cat(as.character(x$ID2),"\n")
        
        if(as.character(x$ID) == "IL-2-at-16"){
            
            cmd0 <- paste("gzip -d",x$file.2,"-c >",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),sep = " ")
            
            cmd1 <- paste("python ~/Dropbox/Aimin_project/AtacSeq/inst/deseq2_wrappers/workflow_for_raw_counts/shift_atac_tagalign.py",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),file.path(output.dir,paste0(x$ID2,"slopped")),sep = " ")
            
            cmd2 <- paste("bedtools coverage -counts -a",ppr.merged.bed.file,"-b",file.path(output.dir,paste0(x$ID2,"slopped")),"| cut -f4 >",file.path(output.dir,x$ID2),sep = " ")
            
            cmd <- paste(cmd0,"wait",cmd1,"wait",cmd2,sep = ";")
            
        }else
        {
            cmd0 <- paste("gzip -d",x$file.2,"-c >",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),sep = " ")
            cmd1 <- paste("python ~/Dropbox/Aimin_project/AtacSeq/inst/deseq2_wrappers/workflow_for_raw_counts/shift_atac_tagalign.py",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),file.path(output.dir,paste0(x$ID2,"slopped")),sep = " ")
            cmd2 <- paste("bedtools coverage -counts -a",ppr.merged.bed.file,"-b",file.path(output.dir,paste0(x$ID2,"slopped")),"| cut -f4 >",file.path(output.dir,x$ID2),sep = " ")
            cmd <- paste(cmd0,"wait",cmd1,"wait",cmd2,sep = ";")
        }
        print(cmd)
        system(cmd)
    },file.5)
}


# input.count.file <- "~/Dropbox/Alejandro_AtacSeq_out_5/count.txt"
# output.dir <- "~/Dropbox/Alejandro_AtacSeq_out_5"
# res <- AtacSeq:::get.DE.pvalue.4.pair.comparision(input.count.file,output.dir)
#
get.DE.pvalue.4.pair.comparision <- function(input.count.file,output.dir) {
    
    # Metadata
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }
    
    register(MulticoreParam(32))
    parallelFlag <- TRUE
    
    countData=data.frame(read.table(input.count.file,header=T,sep='\t'))
    
    counts.metadata <- cbind.data.frame(Sample = colnames(countData),Condition=colnames(countData),Stage=trimGene3(colnames(countData)),stringsAsFactors=F)
    
    counts.metadata$Condition <- trimGene4(counts.metadata$Condition)
    
    counts.metadata.1 <- cbind.data.frame(counts.metadata,ID=paste0(counts.metadata$Condition,".at.",str_sub(counts.metadata$Sample,str_locate(counts.metadata$Sample,"at")[,2]+2,str_length(counts.metadata$Sample))),stringsAsFactors=F)
    
    if(unique(counts.metadata.1$Stage) == 0.75){
        counts.metadata.1[which(counts.metadata.1$Stage == 0.75),]$Stage <- "early"
    }else
    {
        counts.metadata.1[which(counts.metadata.1$Stage == 16),]$Stage <- "late"
    }
    
    counts.metadata.2 <- cbind.data.frame(Sample=counts.metadata.1$ID,Condition=counts.metadata.1$Condition,Stage=counts.metadata.1$Stage,stringsAsFactors=F)
    colnames(countData) <- counts.metadata.2$Sample
    any(is.na(countData))
    
    counts.metadata.2$group <- factor(paste0(counts.metadata.2$Condition,"_",counts.metadata.2$Stage))
    
    ddsFullCountTable<-DESeqDataSetFromMatrix(
        countData=countData,
        colData=counts.metadata.2,
        design= ~ group)
    
    dds<-DESeq(ddsFullCountTable,parallel=parallelFlag,betaPrior=FALSE)
    resultsNames(dds)
    
    ct <- unique(counts.metadata.1$Condition)
    st <- unique(counts.metadata.1$Stage)
    
    ct1 <- paste0(ct[1],"_",st)
    ct2 <- paste0(ct[2],"_",st)
    
    res.DE=results(dds, contrast = c("group",ct1,ct2))
    
    pdf(file.path(output.dir,"Pair_DE_pvalue.pdf"))
    hist(res.DE$pvalue)
    dev.off()
    
    re <- list(dds = dds,res.DE = res.DE)
    
    re
    
}

# peakFile <- "~/Dropbox/Alejandro_AtacSeq_out"

# peakFile <-  "/Users/axy148/Dropbox/Alejandro_AtacSeq_out_All_Peaks"
# out.dir.name <- "~/Dropbox/Alejandro_AtacSeq_out_3"
# AtacSeq:::makeHeatMap2(peakFile,res,out.dir.name,0.05,0.58)

makeHeatMap2 <- function(peakFile,res,out.dir.name,padjust_thresh,lfc_thresh) {
    
    if (!dir.exists(out.dir.name))
    {
        dir.create(out.dir.name, recursive = TRUE)
    }
    
    res.1 <- res[2]
    
    peakNames <- read.table(peakFile,header=F,colClasses = c("character"))
    
    #u <- 1
    
    cmd.l <- lapply(1:length(res.1), function(u,res.1,res,out.dir.name){
        
        x <- names(res.1)[u]
        x.1 <- str_sub(x,str_locate(x,"\\.")[2]+1,str_locate(x,"vs")[1]-2)
        x.2 <- str_sub(x,str_locate(x,"vs")[2]+2,str_locate(x,"at")[1]-2)
        x.3 <- str_sub(x,str_locate(x,"at")[2]+2,str_length(x))
        
        x.name <- paste(x.1,x.2,x.3,sep="-")
        
        if(x.3 == "early"){
            y1 <- paste(x.1,"0.75",sep=".at.")
            y2 <- paste(x.2,"0.75",sep=".at.")
            
        }else
        {
            y1 <- paste(x.1,"16",sep=".at.")
            y2 <- paste(x.2,"16",sep=".at.")
        }
        
        alldata=data.frame(peakNames,as.data.frame(res.1[[u]]))
        
        index <- c(grep(y1,colnames(counts(res$dds,normalized = T))),grep(y2,colnames(counts(res$dds,normalized = T))))
        
        res.with.counts <- cbind.data.frame(peakName= seq(1,dim(alldata)[1]),alldata,counts(res$dds,normalized = T)[,index],stringsAsFactors=F)
        
        res.with.counts <- res.with.counts[order(res.with.counts$pvalue),]
        
        write.csv(res.with.counts,file = file.path(out.dir.name,paste0("DE_",x.name,".csv")))
        
        res.with.counts.1 <- res.with.counts[which(res.with.counts$padj<padjust_thresh),]
        
        res.with.counts.2 <- res.with.counts.1[which(abs(res.with.counts.1$log2FoldChange) > lfc_thresh),]
        
        res.with.counts.3 <- res.with.counts.2[order(res.with.counts.2$padj),]
        
        res.with.counts.4 <- as.matrix(res.with.counts.3[,-c(1:10)])
        
        n = dim(res.with.counts.4)[1]
        
        cat(x.name," ",n," ")
        
        if(n > 30){
            s.n = 30
            cat("OK","\n")
            pdf(file.path(out.dir.name,paste0("HeatMap_",x.name,".pdf")))
            heatmap.2(res.with.counts.4[1:s.n,],margins = c(15,5),dendrogram="column",Colv = T, Rowv=F,col=rev(heat.colors(16)),cexRow=0.5,cexCol=1,sepcolor="white", trace="none")
            dev.off()
        }
        else{
            cat("less than 30 to satisfy the selection thresh","\n")
        }
        
    },res.1,res,out.dir.name)
    
}


# input.atac.dir <- "/Users/axy148/.CMVolumes/AiminDropbox/Alejandro_AtacSeq_uploaded_2"
# output.dir <- "~/Dropbox/Alejandro_AtacSeq_out_All_Peaks"
# ppr.merged.bed.dir <- "/Users/axy148/Dropbox/Alejandro_AtacSeq_out_All_Peaks"
#

AnalysisDEall <- function(input.atac.dir,ppr.merged.bed.dir,outputDir) {
    
    ppr.merged.bed.file <- list.files( path = ppr.merged.bed.dir,pattern=".bed$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    #peakFile <- ppr.merged.bed.file[-c(1,3)]
    
    peakFile <- ppr.merged.bed.file[c(4,5,6)]
    
    cmdL <- lapply(peakFile, function(u){
        
        #u <- "/Users/axy148/Dropbox/Alejandro_AtacSeq_out_All_Peaks/FP_IL2_al_0.75.bed"
        XName <- basename(u)
        XName.1 <- str_sub(XName,1,str_locate_all(XName,"_")[[1]][1,1]-1)
        XName.2 <- str_sub(XName,str_locate_all(XName,"_")[[1]][1,2]+1,str_locate_all(XName,"_")[[1]][2,1]-1)
        XName.3 <- str_sub(XName,str_locate_all(XName,"_")[[1]][3,2]+1,str_locate_all(XName,".bed")[[1]][1,1]-1)
        
        output.dir <- file.path(outputDir,paste0("out",XName.1,XName.2,XName.3,sep="-"))
        
        if (!dir.exists(output.dir))
        {
            dir.create(output.dir, recursive = TRUE)
        }
        
        getCount3(input.atac.dir,u,output.dir)
        
        cmd <- paste("ls -lrth",file.path(output.dir,"*rep?"),"| awk '{print $9}'",sep=" ")
        makeCountTable3(cmd,output.dir)
        
        input.count.file <- file.path(output.dir,"count.txt")
        res <- get.DE.pvalue.4.pair.comparision(input.count.file,output.dir)
        
        makeHeatMap3(u,res,output.dir,0.05,0.58)
        
    })
    
}


# peakFile <- "~/Dropbox/Alejandro_AtacSeq_out_4/FP_PBS_at_0.75.bed"

# out.dir.name <- "~/Dropbox/Alejandro_AtacSeq_common_peaks"

# AtacSeq:::makeHeatMap3(peakFile,res,out.dir.name,0.05,0.58)

makeHeatMap3 <- function(peakFile,res,out.dir.name,padjust_thresh,lfc_thresh) {
    
    if (!dir.exists(out.dir.name))
    {
        dir.create(out.dir.name, recursive = TRUE)
    }
    
    res.1 <- res[2]
    
    peakNames <- read.table(peakFile,header=F,colClasses = c("character"))
    
    names(res.1) <- basename(peakFile)
    
    u <- 1
    
    cmd.l <- lapply(1:length(res.1), function(u,res.1,res,out.dir.name){
        
        x <- names(res.1)[u]
        
        x.1 <- str_sub(x,1,str_locate_all(x,"_")[[1]][1,1]-1)
        
        x.2 <- str_sub(x,str_locate_all(x,"_")[[1]][1,2]+1,str_locate_all(x,"_")[[1]][2,1]-1)
        
        x.3 <- str_sub(x,str_locate_all(x,"_")[[1]][3,1]+1,str_locate_all(x,".bed")[[1]][1,1]-1)
        
        x.name <- paste(x.1,x.2,x.3,sep="-")
        
        if(x.3 == "0.75"){
            y1 <- paste(x.1,"0.75",sep=".at.")
            y2 <- paste(x.2,"0.75",sep=".at.")
            
        }else
        {
            y1 <- paste(x.1,"16",sep=".at.")
            y2 <- paste(x.2,"16",sep=".at.")
        }
        
        alldata=data.frame(peakNames,as.data.frame(res.1[[u]]))
        
        index <- c(grep(y1,colnames(counts(res$dds,normalized = T))),grep(y2,colnames(counts(res$dds,normalized = T))))
        
        res.with.counts <- cbind.data.frame(peakName= seq(1,dim(alldata)[1]),alldata,counts(res$dds,normalized = T)[,index],stringsAsFactors=F)
        
        res.with.counts <- res.with.counts[order(res.with.counts$pvalue),]
        
        write.csv(res.with.counts,file = file.path(out.dir.name,paste0("DE_",x.name,".csv")))
        
        res.with.counts.1 <- res.with.counts[which(res.with.counts$padj<padjust_thresh),]
        
        res.with.counts.2 <- res.with.counts.1[which(abs(res.with.counts.1$log2FoldChange) > lfc_thresh),]
        
        res.with.counts.3 <- res.with.counts.2[order(res.with.counts.2$padj),]
        
        res.with.counts.4 <- as.matrix(res.with.counts.3[,-c(1:10)])
        
        n = dim(res.with.counts.4)[1]
        
        cat(x.name," ",n," ")
        
        if(n > 30){
            s.n = 30
            cat("OK","\n")
            pdf(file.path(out.dir.name,paste0("HeatMap_",x.name,".pdf")))
            heatmap.2(res.with.counts.4[1:s.n,],margins = c(15,5),dendrogram="column",Colv = T, Rowv=F,col=rev(heat.colors(16)),cexRow=0.5,cexCol=1,sepcolor="white", trace="none")
            dev.off()
        }
        else{
            cat("less than 30 to satisfy the selection thresh","\n")
        }
        
    },res.1,res,out.dir.name)
    
}

# input.atac.dir <- "/Users/axy148/.CMVolumes/AiminDropbox/Alejandro_AtacSeq_uploaded_2"
# output.dir <- "~/Dropbox/Alejandro_AtacSeq_out_5"
# ppr.merged.bed.file <- "/Users/axy148/Dropbox/Alejandro_AtacSeq_out/FP_PBS_at_0.75.bed"
# getCount3.2(input.atac.dir,ppr.merged.bed.file,output.dir)

getCount3.2 <- function(input.atac.dir,ppr.merged.bed.file,output.dir) {
    
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }
    
    #file.1 <- list.files(input.atac.dir,pattern="*tn5_pooled.tagAlign.gz$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    file.1 <- list.files(input.atac.dir,pattern="*tn5.tagAlign.gz$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    file.2 <- file.1[-grep("pseudo_reps",file.1)]
    
    #file.2 <- file.1[grep("pooled_pseudo_reps",file.1)]
    
    file.3 <- cbind.data.frame(ID = basename((dirname(dirname(dirname(file.2))))),rep= basename(dirname(file.2)),fileName=basename(file.2),file.2,stringsAsFactors=FALSE)
    
    file.4 <- file.3[-which(file.3$ID %in% c("IL-2vsPBS")),]
    
    file.5 <- cbind.data.frame(ID2 = paste0(file.4$ID,"_",file.4$rep),file.4,stringsAsFactors=FALSE)
    
    file.5[grep("IL-2-at-16_rep1",file.5$ID2),]$ID2 <-  "IL2-at-16_rep1"
    file.5[grep("IL-2-at-16_rep2",file.5$ID2),]$ID2 <-  "IL2-at-16_rep2"
    file.5[grep("IL-2-at-16_rep3",file.5$ID2),]$ID2 <-  "IL2-at-16_rep3"
    
    file.5[grep("IL-2-at-16",file.5$ID),]$ID <-  "IL2-at-16"
    
    
    file.id <- unique(as.character(file.5$ID2))
    
    ppr.merged.bed.file="/Users/axy148/Dropbox/Alejandro_AtacSeq_out_All_Peaks/FP_IL2_al_16.bed"
    
    XName <- basename(ppr.merged.bed.file)
    
    XName.1 <- str_sub(XName,1,str_locate_all(XName,"_")[[1]][1,1]-1)
    
    XName.2 <- str_sub(XName,str_locate_all(XName,"_")[[1]][1,2]+1,str_locate_all(XName,"_")[[1]][2,1]-1)
    
    XName.3 <- str_sub(XName,str_locate_all(XName,"_")[[1]][3,2]+1,str_locate_all(XName,".bed")[[1]][1,1]-1)
    
    p1 <- paste(XName.1,"at",XName.3,sep="-")
    p2 <- paste(XName.2,"at",XName.3,sep="-")
    
    file.id <- file.id[c(grep(p1,file.id),grep(p2,file.id))]
    output.dir <- "/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/outFPIL216-"
    
    cmd.l <- lapply(file.id, function(u,file.5){
        
        
        x <- file.5[which(as.character(file.5$ID2) == u),]
        
        cat("\n")
        cat(as.character(x$ID2),"\n")
        
        if(as.character(x$ID) == "IL-2-at-16"){
            
            cmd0 <- paste("gzip -d",x$file.2,"-c >",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),sep = " ")
            
            cmd1 <- paste("python ~/Dropbox/Aimin_project/AtacSeq/inst/deseq2_wrappers/workflow_for_raw_counts/shift_atac_tagalign.py",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),file.path(output.dir,paste0(x$ID2,"slopped")),sep = " ")
            
            cmd2 <- paste("bedtools coverage -counts -a",ppr.merged.bed.file,"-b",file.path(output.dir,paste0(x$ID2,"slopped")),"| cut -f4 >",file.path(output.dir,x$ID2),sep = " ")
            
            cmd <- paste(cmd0,"wait",cmd1,"wait",cmd2,sep = ";")
            
        }else
        {
            cmd0 <- paste("gzip -d",x$file.2,"-c >",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),sep = " ")
            cmd1 <- paste("python ~/Dropbox/Aimin_project/AtacSeq/inst/deseq2_wrappers/workflow_for_raw_counts/shift_atac_tagalign.py",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),file.path(output.dir,paste0(x$ID2,"slopped")),sep = " ")
            cmd2 <- paste("bedtools coverage -counts -a",ppr.merged.bed.file,"-b",file.path(output.dir,paste0(x$ID2,"slopped")),"| cut -f4 >",file.path(output.dir,x$ID2),sep = " ")
            cmd <- paste(cmd0,"wait",cmd1,"wait",cmd2,sep = ";")
        }
        print(cmd)
        system(cmd)
    },file.5)
}

# anno.dir = "/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/Annotation/Promoter_5UTR_3UTR_Exon_Intron_Downstream_Intergenic"
# de.dir = "/Volumes/Bioinformatics/Aimin_project/ATAC-Seq"
#
# res.anno.de <- mergeDeWithAnno(anno.dir,de.dir,"/Volumes/Bioinformatics/Aimin_project/ATAC-Seq")

mergeDeWithAnno <- function(anno.dir,de.dir,out.dir.name){
    
    file.1 <- list.files(anno.dir,pattern="*all_peaks.xls$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    file.2 <- list.files(de.dir,pattern="DE.*csv",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    temp1 <- cbind.data.frame(anno.file=file.1,anno.file.name=basename(file.1),ID=str_replace(str_sub(basename(file.1),1,str_locate(basename(file.1),"bed_5000")[,1]-2),"_al",""),stringsAsFactors=FALSE)
    
    temp2 <- cbind.data.frame(de.file=file.2,de.file.name=basename(file.2),ID=str_replace_all(str_sub(basename(file.2),4,str_locate(basename(file.2),"csv")[,1]-2),"-","_"),stringsAsFactors=FALSE)
    
    temp3 <- merge.data.frame(temp1,temp2,sort=FALSE,by="ID")
    
    res <- apply(temp3,1,function(u){
        
        x <- as.data.frame(t(u))
        
        anno <- read.table(as.character(x$anno.file),header=T,sep="\t",quote = "",row.names = NULL, stringsAsFactors = FALSE)
        
        de <- read.csv(as.character(x$de.file))
        
        #colnames(de)[2] = "peakName"
        anno.de <-  merge.data.frame(anno,de,by="peakName",sort=FALSE)
        anno.de.1 <- anno.de[,-19]
        
        sub("X.","",colnames(anno.de.1)[19:31])
        colnames(anno.de.1)[19:21] <- c("chr","start","end")
        
        print(head(anno.de.1))
        
        x_name <-as.character(x$ID)
        
        write.table(anno.de.1,file=file.path(out.dir.name,paste0(x_name,"_annotation_DE.xls")),row.names = FALSE,quote=FALSE,sep="\t")
        
        anno.de.1
        
    })
    
    names(res) <- temp3$ID
    res
}

filterPeaks <- function(peak.set,padjust_thresh,lfc_thresh) {
    
    peak.set.1 <- peak.set[which(peak.set$padj<padjust_thresh),]
    
    peak.set.2 <- peak.set.1[which(abs(peak.set.1$log2FoldChange) > lfc_thresh),]
    
    peak.set.3 <- peak.set.2[order(peak.set.2$padj),]
    
    peak.set.3
}

# res.violin <- AtacSeq:::makeViolinplot(res.anno.de,"/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/Violinplot")
#
makeViolinplot <- function(res.anno.de,output.dir){
    
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }
    
    name.trt <- as.character(names(res.anno.de))
    
    cmd.l <- lapply(1:length(res.anno.de), function(u,res.anno.de,name.trt){
        
        temp <- filterPeaks(res.anno.de[[u]],0.05,0.58)
        
        cat(name.trt[u],"\t",dim(temp)[1],"\n")
        
        if(dim(temp)[1]>0){
            
            df <- melt(temp[,28:33])
            colnames(df) <- c("trt","openness")
            
            print(unique(df$trt))
            
            x <- df$trt
            
            #x.2 <- str_sub(x,str_locate_all(x,"_")[[1]][1,2]+1,str_locate_all(x,"_")[[1]][2,1]-1)
            
            x.1 <- str_sub(x,1,str_locate_all(x,"\\.")[[1]][1,2]-1)
            
            x.2 <- str_sub(x,unlist(lapply(str_locate_all(x,"_"),"[",c(2))),str_length(x))
            
            x.3 <- str_sub(x,unlist(lapply(str_locate_all(x,"\\."),"[",c(2)))+1,unlist(lapply(str_locate_all(x,"_"),"[",c(1)))-1)
            
            #print(x.1)
            #print(x.2)
            #print(x.3)
            
            x.4 <- paste0(x.1,x.2)
            
            df.1 <- data.frame(gr=x.4,time=x.3,openness=df$openness,cell=x.1)
            
            png(file.path(output.dir,paste0(names(res.anno.de)[u],".png")))
            
            #plot <- ggplot2.violinplot(data=df.1,xName='gr', yName='openness',
            #                           mainTitle=unique(df.1$time),xtitle="trt",ytitle="openness")
            library(dplyr)
            df.1Summary <- df.1 %>%
                group_by(gr) %>%
                summarize(openness_mean = mean(openness),
                          openness_se = sqrt(var(openness)/length(openness)))
            
            plotViolins <- ggplot(df.1, aes(x = factor(gr), y = openness, fill = factor(gr))) +
                geom_violin(aes(y = openness, fill = factor(cell))) +
                geom_point(aes(y = openness_mean), color = "black", size = 2, data = df.1Summary) +
                geom_errorbar(aes(y = openness_mean, ymin = openness_mean-openness_se, ymax = openness_mean+openness_se),
                              color = "black", width = 0.2, data = df.1Summary) +
                ylim(0, 2000) +
                theme(legend.position = "none") +
                ggtitle(paste0("At_",unique(df.1$time),"_Time_Point")) + theme(plot.title=element_text(hjust=0.5))+xlab("Trt") + ylab("Chromatin Accessibility")
            print(plotViolins)
            
            dev.off()
        }
        
    },res.anno.de,name.trt)
    
}

sortPeaks <- function(peak.set,top.peak.num) {
    
    peak.set.1 <- peak.set[order(abs(peak.set$log2FoldChange),decreasing = TRUE),]
    
    peak.set.2 <- peak.set.1[order(peak.set.1$padj),]
    
    peak.set.3 <- peak.set.2[1:top.peak.num,]
    
    peak.set.3
}

# res.violin2 <- AtacSeq:::makeViolinplot2(res.anno.de,882,"/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/Violinplot2")

makeViolinplot2 <- function(res.anno.de,top.peak.num,output.dir){
    
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }
    
    name.trt <- as.character(names(res.anno.de))
    
    cmd.l <- lapply(1:length(res.anno.de), function(u,res.anno.de,name.trt,top.peak.num,output.dir){
        
        temp <- sortPeaks(res.anno.de[[u]],top.peak.num)
        
        cat(name.trt[u],"\t",dim(temp)[1],"\n")
        
        if(dim(temp)[1]>0){
            
            df <- melt(temp[,28:33])
            colnames(df) <- c("trt","openness")
            
            print(unique(df$trt))
            
            x <- df$trt
            
            #x.2 <- str_sub(x,str_locate_all(x,"_")[[1]][1,2]+1,str_locate_all(x,"_")[[1]][2,1]-1)
            
            x.1 <- str_sub(x,1,str_locate_all(x,"\\.")[[1]][1,2]-1)
            
            x.2 <- str_sub(x,unlist(lapply(str_locate_all(x,"_"),"[",c(2))),str_length(x))
            
            x.3 <- str_sub(x,unlist(lapply(str_locate_all(x,"\\."),"[",c(2)))+1,unlist(lapply(str_locate_all(x,"_"),"[",c(1)))-1)
            
            #print(x.1)
            #print(x.2)
            #print(x.3)
            
            x.4 <- paste0(x.1,x.2)
            
            df.1 <- data.frame(gr=x.4,time=x.3,openness=df$openness,cell=x.1)
            
            png(file.path(output.dir,paste0(names(res.anno.de)[u],".png")))
            
            #plot <- ggplot2.violinplot(data=df.1,xName='gr', yName='openness',
            #                           mainTitle=unique(df.1$time),xtitle="trt",ytitle="openness")
            library(dplyr)
            df.1Summary <- df.1 %>%
                group_by(gr) %>%
                summarize(openness_mean = mean(openness),
                          openness_se = sqrt(var(openness)/length(openness)))
            
            plotViolins <- ggplot(df.1, aes(x = factor(gr), y = openness, fill = factor(gr))) +
                geom_violin(aes(y = openness, fill = factor(cell))) +
                geom_point(aes(y = openness_mean), color = "black", size = 2, data = df.1Summary) +
                geom_errorbar(aes(y = openness_mean, ymin = openness_mean-openness_se, ymax = openness_mean+openness_se),
                              color = "black", width = 0.2, data = df.1Summary) +
                ylim(0, 2000) +
                theme(legend.position = "none") +
                ggtitle(paste0("At_",unique(df.1$time),"_Time_Point")) + theme(plot.title=element_text(hjust=0.5))+xlab("Trt") + ylab("Chromatin Accessibility")
            print(plotViolins)
            
            dev.off()
        }
        
    },res.anno.de,name.trt,top.peak.num,output.dir)
    
}


# input.rna.seq.dir <- "/Volumes/Bioinformatics-1/Aimin_project/ATAC-Seq/RNA-Seq"
# input.file.pattern <- "*xlsx$"
# rna.data <- processRNASeq(input.rna.seq.dir,input.file.pattern)

processRNASeq  <- function(input.rna.seq.dir,input.file.pattern) {
    
    #file.1 <- list.files(input.atac.dir,pattern="*tn5_pooled.tagAlign.gz$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    file.1 <- list.files(input.rna.seq.dir,pattern=input.file.pattern,all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    file.2 <- gsub(" ", "-",basename(file.1))
    
    #file.2 <- file.1[grep("pooled_pseudo_reps",file.1)]
    ID <- str_sub(file.2,unlist(lapply(str_locate_all(file.2,"-"),"[",c(1)))+1,unlist(lapply(str_locate_all(file.2,"xlsx"),"[",c(1)))-2)
    
    file.3 <- cbind.data.frame(ID = ID,fileName=file.1,stringsAsFactors=FALSE)
    
    print(file.3)
    library(gdata)
    
    ID <- file.3$ID
    fileName <- file.3$fileName
    
    cmd.l <- lapply(1:length(fileName), function(u,fileName,ID){
        
        rna.data <- read.xls(fileName[u])
        
        rna.data
        
    },fileName,ID)
    
    names(cmd.l) <- ID
    cmd.l
}

# FP.PBS.at.16.rna.atac <- mergeRnaSeqAtacSeq(rna.data,res.anno.de,4,4)
# IL2.PBS.at.16.rna.atac <- mergeRnaSeqAtacSeq(rna.data,res.anno.de,2,6)
#
# IL2.PBS.at.0.75.rna.atac <- mergeRnaSeqAtacSeq(rna.data,res.anno.de,1,5)
#
# FP.PBS.at.0.75.rna.atac <- mergeRnaSeqAtacSeq(rna.data,res.anno.de,3,3)
#
#
#
mergeRnaSeqAtacSeq <- function(rna.data, res.anno.de,rna.index,atac.index) {
    x <- rna.data[[rna.index]]$ensembl_gene_id
    FP.PBS.at.16 <- cbind.data.frame(ENSEMBL = str_sub(x,1,str_locate_all(x,"\\.")[[1]][1,2]-1),rna.data[[rna.index]])
    FP.PBS.at.16.rna.atac <- merge(res.anno.de[[atac.index]],FP.PBS.at.16,by="ENSEMBL")
}


#
# FP at 16
# makeContour(FP.PBS.at.16.rna.atac,index=c(46,49,28,30))
#
# PBS at 16
# makeContour(FP.PBS.at.16.rna.atac,index=c(50,53,31,33))
#
# IL2 at 16
# makeContour(IL2.PBS.at.16.rna.atac,index=c(46,49,28,30))
#
# PBS at 16
# makeContour(IL2.PBS.at.16.rna.atac,index=c(50,53,31,33))
#
# IL2 at 0.75
# makeContour(IL2.PBS.at.0.75.rna.atac,index=c(46,49,28,30))
#
# PBS at 0.75
# makeContour(IL2.PBS.at.0.75.rna.atac,index=c(50,53,31,33))
#
# FP at 0.75
# makeContour(FP.PBS.at.0.75.rna.atac,index=c(46,49,28,30))
#
#
makeContour <- function(FP.PBS.at.16.rna.atac,index,rpkm=FALSE) {
    
    i1 <- index[1]
    i2 <- index[2]
    i3 <- index[3]
    i4 <- index[4]
    
    countmatrix.rna <- as.matrix(FP.PBS.at.16.rna.atac[,i1:i2]) + 1
    cpm <- apply(countmatrix.rna,2, function(x) (x/sum(x))*1000000)
    log.cpm.rna <- log10(cpm)
    
    # if using rpkm
    if (rpkm == TRUE){
        expression <- data.frame(FP.PBS.at.16.rna.atac[,i1:i2]+1,geneLength=FP.PBS.at.16.rna.atac$V5-FP.PBS.at.16.rna.atac$V4)
        geneLength <- expression$geneLength
        expression.rpkm <- apply(expression[,1:8],2,function(u,geneLength){
            x <- (u/geneLength)*(10^9/sum(u))
        },geneLength)
        log.cpm.rna <-log10(expression.rpkm)
    }
    
    countmatrix.atac <- as.matrix(FP.PBS.at.16.rna.atac[,i3:i4]) + 1
    log.cpm.atac <- log10(countmatrix.atac)
    
    
    Y <- cbind.data.frame(ENSEMBL=FP.PBS.at.16.rna.atac$ENSEMBL,log.cpm.rna,log.cpm.atac)
    
    FP.rna.atac <- cbind.data.frame(FP.rna=apply(Y[,2:5],1,mean), Fp.atac = apply(Y[,6:8],1,mean))
    
    plot(FP.rna.atac$FP.rna~FP.rna.atac$Fp.atac)
    cor.test(FP.rna.atac$FP.rna,FP.rna.atac$Fp.atac)
    print(cor(FP.rna.atac))
    
    dataset <- FP.rna.atac
    colnames(dataset)=c("x","y")
    dataset2 <- with(dataset, dataset[order(x,y), ])
    
    ggplot(dataset2, aes(y, x)) + ylim(min(dataset2$x),max(dataset2$x)) + xlim(min(dataset2$y),max(dataset2$y))+ stat_density2d(aes(fill = ..level..),geom="polygon") + ggtitle("RNA-seq vs chromatin openness for FP At 16 hours") + theme(plot.title=element_text(hjust=0.5)) + xlab("Openness [log10 (ATAC cut sites+1)]") + ylab("Expression [log10 (CPM)]")
    
}

#ggplot(dataset2, aes(y, x)) + stat_density2d(aes(fill = ..level..),geom="polygon")
# ggplot(dataset2, aes(y, x)) +
#   stat_density2d(aes(alpha=..level.., fill=..level..), size=2, bins=10, geom="polygon") +
#   scale_fill_gradient(low = "yellow", high = "red") +
#   scale_alpha(range = c(0.00, 0.5), guide = FALSE) +
#   geom_density2d(colour="black") +
#   geom_point(data = dataset2) +
#   guides(alpha=FALSE)


compareAtacPvalueWithRnaSeqPvalue <- function() {
    IL2.PBS.at.0.75.rna.atac.2 <- IL2.PBS.at.0.75.rna.atac[order(IL2.PBS.at.0.75.rna.atac$padj.x),]
    plot(IL2.PBS.at.0.75.rna.atac.2$padj.x,IL2.PBS.at.0.75.rna.atac.2$padj.y)
    
    cor(IL2.PBS.at.0.75.rna.atac.2$padj.x,IL2.PBS.at.0.75.rna.atac.2$padj.y,use="pairwise.complete.obs",method="spearman")
}


filterPeaks2 <- function(peak.set,padjust_thresh,lfc_thresh) {
    
    peak.set.1 <- peak.set[which(peak.set$padj.x<padjust_thresh),]
    
    peak.set.2 <- peak.set.1[which(abs(peak.set.1$log2FoldChange.x) > lfc_thresh),]
    
    peak.set.3 <- peak.set.2[order(peak.set.2$padj.x),]
    
    peak.set.3
}



# IL2.PBS.at.0.75.rna.atac[IL2.PBS.at.0.75.rna.atac$peakName==17338,]

# makeContour4Gene(IL2.PBS.at.0.75.rna.atac,atac.index=c(28:33),rna.index=c(46:53),peak.index=17338,rpkm=FALSE,output.dir)

# output.dir <- "/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/RNASeq-AtacSeq/geneBased"

# IL2.PBS.at.0.75.rna.atac.countour.single.gene <-  makeContour4Gene(IL2.PBS.at.0.75.rna.atac,atac.index=c(28:33),rna.index=c(46:53),peak.index=17338,rpkm=FALSE,output.dir)

makeContour4Gene <- function(IL2.PBS.at.0.75.rna.atac,atac.index,rna.index,peak.index,rpkm=FALSE,output.dir) {
    
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }
    
    countmatrix.rna <- as.matrix(IL2.PBS.at.0.75.rna.atac[,rna.index]) + 1
    cpm <- apply(countmatrix.rna,2, function(x) (x/sum(x))*1000000)
    log.cpm.rna <- log10(cpm)
    
    # if using rpkm
    if (rpkm == TRUE){
        expression <- data.frame(IL2.PBS.at.0.75.rna.atac[,rna.index]+1,geneLength=IL2.PBS.at.0.75.rna.atac$V5-IL2.PBS.at.0.75.rna.atac$V4)
        geneLength <- expression$geneLength
        expression.rpkm <- apply(expression[,1:8],2,function(u,geneLength){
            x <- (u/geneLength)*(10^9/sum(u))
        },geneLength)
        log.cpm.rna <-log10(expression.rpkm)
    }
    
    countmatrix.atac <- as.matrix(IL2.PBS.at.0.75.rna.atac[,atac.index]) + 1
    log.cpm.atac <- log10(countmatrix.atac)
    
    
    Y <- cbind.data.frame(ENSEMBL=IL2.PBS.at.0.75.rna.atac$ENSEMBL,GeneName =IL2.PBS.at.0.75.rna.atac$gene_id,peakName=IL2.PBS.at.0.75.rna.atac$peakName,log.cpm.rna,log.cpm.atac)
    
    print(Y)
    Y
    YY <- Y[Y$peakName==peak.index,]
    YY
    
    FP.rna.atac <- data.frame(measurement=as.numeric(YY[,4:17]),type=c(rep("rna",8),rep("atac",6)),pt=c(rep("IL2",4),rep("PBS",4),rep("IL2",3),rep("PBS",3)),time=rep("0.75",14),gene=rep("Phc3",14))
    
    #FP.rna.atac
    df.1 <- data.frame(gr=paste0(FP.rna.atac$pt,"_",FP.rna.atac$type),time=FP.rna.atac$time,openness=FP.rna.atac$measurement,cell=FP.rna.atac$type)
    
    pdf(file.path(output.dir,paste0(unique(as.character(FP.rna.atac$gene)),".pdf")))
    
    #plot <- ggplot2.violinplot(data=df.1,xName='gr', yName='openness',
    #                           mainTitle=unique(df.1$time),xtitle="trt",ytitle="openness")
    library(dplyr)
    df.1Summary <- df.1 %>%
        group_by(gr) %>%
        summarize(openness_mean = mean(openness),
                  openness_se = sqrt(var(openness)/length(openness)))
    
    plotViolins <- ggplot(df.1, aes(x = factor(gr), y = openness, fill = factor(gr))) +
        geom_violin(aes(y = openness, fill = factor(cell))) +
        geom_point(aes(y = openness_mean), color = "black", size = 2, data = df.1Summary) +
        geom_errorbar(aes(y = openness_mean, ymin = openness_mean-openness_se, ymax = openness_mean+openness_se),
                      color = "black", width = 0.2, data = df.1Summary) +
        ylim(1, 3) +
        theme(legend.position = "none") +
        ggtitle(paste0(unique(as.character(FP.rna.atac$gene)),"_At_",unique(df.1$time),"_Time_Point")) + theme(plot.title=element_text(hjust=0.5))+xlab("Trt") + ylab("Gene exprepression difference regulated by different chromatin accessibility")
    print(plotViolins)
    
    dev.off()
    
}


# IL2.PBS.at.0.75.rna.atac.2.filtered <- filterPeaks2(IL2.PBS.at.0.75.rna.atac.2,0.05,0.58)
# IL2.PBS.at.0.75.rna.atac.2.filtered.2 <- IL2.PBS.at.0.75.rna.atac.2.filtered[which(IL2.PBS.at.0.75.rna.atac.2.filtered$log2FoldChange.x>0&IL2.PBS.at.0.75.rna.atac.2.filtered$log2FoldChange.y>0),]
# IL2.PBS.at.0.75.rna.atac.2.filtered.3 <- IL2.PBS.at.0.75.rna.atac.2.filtered.2[which(IL2.PBS.at.0.75.rna.atac.2.filtered.2$padj.y<0.05),]
# output.dir <- "~/ATAC-Seq/RNASeq-AtacSeq/geneBased"
# IL2.PBS.at.0.75.rna.atac.countour.mutiple.gene <-  makeContour4Gene2(IL2.PBS.at.0.75.rna.atac.2.filtered.3,atac.index=c(28:33),rna.index=c(46:53),peak.index=IL2.PBS.at.0.75.rna.atac.2.filtered.3$peakName,rpkm=FALSE,output.dir)

makeContour4Gene2 <- function(IL2.PBS.at.0.75.rna.atac,atac.index,rna.index,peak.index,rpkm=FALSE,output.dir,s1,s2,tp) {
    
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }
    
    countmatrix.rna <- as.matrix(IL2.PBS.at.0.75.rna.atac[,rna.index]) + 1
    cpm <- apply(countmatrix.rna,2, function(x) (x/sum(x))*1000000)
    log.cpm.rna <- log10(cpm)
    
    # if using rpkm
    if (rpkm == TRUE){
        expression <- data.frame(IL2.PBS.at.0.75.rna.atac[,rna.index]+1,geneLength=IL2.PBS.at.0.75.rna.atac$V5-IL2.PBS.at.0.75.rna.atac$V4)
        geneLength <- expression$geneLength
        expression.rpkm <- apply(expression[,1:8],2,function(u,geneLength){
            x <- (u/geneLength)*(10^9/sum(u))
        },geneLength)
        log.cpm.rna <-log10(expression.rpkm)
    }
    
    countmatrix.atac <- as.matrix(IL2.PBS.at.0.75.rna.atac[,atac.index]) + 1
    log.cpm.atac <- log10(countmatrix.atac)
    
    
    Y <- cbind.data.frame(ENSEMBL=IL2.PBS.at.0.75.rna.atac$ENSEMBL,GeneName =IL2.PBS.at.0.75.rna.atac$gene_id,peakName=IL2.PBS.at.0.75.rna.atac$peakName,log.cpm.rna,log.cpm.atac)
    
    print(Y)
    Y
    YY <- Y[which(Y$peakName %in% peak.index),]
    
    
    gene.id <- YY$GeneName
    
    FP.rna.atac.L <- lapply(1:length(gene.id),function(u,YY,output.dir){
        
        #print(YY[which(YY$GeneName == gene.id[u]),])
        
        n= dim(YY[which(YY$GeneName == gene.id[u]),])[1]
        
        print(c(YY[which(YY$GeneName == gene.id[u]),4:17]))
        
        FP.rna.atac <- data.frame(measurement=as.numeric(unlist(c(YY[which(YY$GeneName == gene.id[u]),4:17]))),type=c(rep("rna",8*n),rep("atac",6*n)),pt=c(rep(s1,4*n),rep(s2,4*n),rep(s1,3*n),rep(s2,3*n)),time=rep(tp,14*n),gene=rep(gene.id[u],14*n))
        
        #FP.rna.atac
        df.1 <- data.frame(gr=paste0(FP.rna.atac$pt,"_",FP.rna.atac$type),time=FP.rna.atac$time,openness=FP.rna.atac$measurement,cell=FP.rna.atac$type)
        
        #pdf(file.path(output.dir,paste0(gene.id[u],".pdf")))
        
        #plot <- ggplot2.violinplot(data=df.1,xName='gr', yName='openness',
        #                           mainTitle=unique(df.1$time),xtitle="trt",ytitle="openness")
        library(dplyr)
        df.1Summary <- df.1 %>%
            group_by(gr) %>%
            summarize(openness_mean = mean(openness),
                      openness_se = sqrt(var(openness)/length(openness)))
        
        plotViolins <- ggplot(df.1, aes(x = factor(gr), y = openness, fill = factor(gr))) +
            geom_violin(aes(y = openness, fill = factor(cell))) +
            geom_point(aes(y = openness_mean), color = "black", size = 2, data = df.1Summary) +
            geom_errorbar(aes(y = openness_mean, ymin = openness_mean-openness_se, ymax = openness_mean+openness_se),
                          color = "black", width = 0.2, data = df.1Summary) +
            ylim(min(df.1$openness), max(df.1$openness)) +
            theme(legend.position = "none") +
            ggtitle(paste0(gene.id[u],"_At_",unique(df.1$time),"_Time_Point")) + theme(plot.title=element_text(hjust=0.5))+xlab("T")+ylab("M")
        #print(plotViolins)
        
        #dev.off()
        
        plotViolins
        
    },YY,output.dir)
    
    names(FP.rna.atac.L)=gene.id
    #FP.rna.atac <- do.call(rbind.data.frame,FP.rna.atac.L)
    
    #gene.id <- unique(as.character(FP.rna.atac$gene))
    
    #lapply()
    
    # df.1 <- data.frame(gr=paste0(FP.rna.atac$pt,"_",FP.rna.atac$type),time=FP.rna.atac$time,openness=FP.rna.atac$measurement,cell=FP.rna.atac$type)
    #
    # pdf(file.path(output.dir,paste0(length(unique(as.character(FP.rna.atac$gene))),".pdf")))
    #
    # #plot <- ggplot2.violinplot(data=df.1,xName='gr', yName='openness',
    # #                           mainTitle=unique(df.1$time),xtitle="trt",ytitle="openness")
    # library(dplyr)
    # df.1Summary <- df.1 %>%
    #   group_by(gr) %>%
    #   summarize(openness_mean = mean(openness),
    #             openness_se = sqrt(var(openness)/length(openness)))
    #
    # plotViolins <- ggplot(df.1, aes(x = factor(gr), y = openness, fill = factor(gr))) +
    #   geom_violin(aes(y = openness, fill = factor(cell))) +
    #   geom_point(aes(y = openness_mean), color = "black", size = 2, data = df.1Summary) +
    #   geom_errorbar(aes(y = openness_mean, ymin = openness_mean-openness_se, ymax = openness_mean+openness_se),
    #                 color = "black", width = 0.2, data = df.1Summary) +
    #   ylim(min(y), max(y)) +
    #   theme(legend.position = "none") +
    #   ggtitle(paste0(length(unique(as.character(FP.rna.atac$gene))),"_At_",unique(df.1$time),"_Time_Point")) + theme(plot.title=element_text(hjust=0.5))+xlab("Trt") + ylab("Gene exprepression difference regulated by different chromatin accessibility")
    # print(plotViolins)
    #
    # dev.off()
    # #df.1
    FP.rna.atac.L
}

makeAtacRnaPlot2 <- function(IL2.PBS.at.0.75.rna.atac.countour.mutiple.gene, output.dir) {
    multi.page <- ggarrange(plotlist=IL2.PBS.at.0.75.rna.atac.countour.mutiple.gene,nrow = 2, ncol = 2)
    ggexport(multi.page, filename = file.path(output.dir,"53_genes.pdf"))
    
    multi.4.page <- ggarrange(plotlist=IL2.PBS.at.0.75.rna.atac.countour.mutiple.gene[c(3,5,9,12)],nrow = 2, ncol = 2)
    ggexport(multi.4.page, filename = file.path(output.dir,"4_genes.pdf"))
}


# IL2.PBS.at.0.75.rna.atac.countour.mutiple.gene <-  makeContour4SelectedGenes(IL2.PBS.at.0.75.rna.atac.2.filtered.3,atac.index=c(28:33),rna.index=c(46:53),peak.index=IL2.PBS.at.0.75.rna.atac.2.filtered.3$peakName,rpkm=FALSE,output.dir)

makeContour4SelectedGenes <- function(IL2.PBS.at.0.75.rna.atac,atac.index,rna.index,peak.index,rpkm=FALSE,output.dir) {
    
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }
    
    countmatrix.rna <- as.matrix(IL2.PBS.at.0.75.rna.atac[,rna.index]) + 1
    cpm <- apply(countmatrix.rna,2, function(x) (x/sum(x))*1000000)
    log.cpm.rna <- log10(cpm)
    
    # if using rpkm
    if (rpkm == TRUE){
        expression <- data.frame(IL2.PBS.at.0.75.rna.atac[,rna.index]+1,geneLength=IL2.PBS.at.0.75.rna.atac$V5-IL2.PBS.at.0.75.rna.atac$V4)
        geneLength <- expression$geneLength
        expression.rpkm <- apply(expression[,1:8],2,function(u,geneLength){
            x <- (u/geneLength)*(10^9/sum(u))
        },geneLength)
        log.cpm.rna <-log10(expression.rpkm)
    }
    
    countmatrix.atac <- as.matrix(IL2.PBS.at.0.75.rna.atac[,atac.index]) + 1
    log.cpm.atac <- log10(countmatrix.atac)
    
    
    Y <- cbind.data.frame(ENSEMBL=IL2.PBS.at.0.75.rna.atac$ENSEMBL,GeneName =IL2.PBS.at.0.75.rna.atac$gene_id,peakName=IL2.PBS.at.0.75.rna.atac$peakName,log.cpm.rna,log.cpm.atac)
    
    print(Y)
    Y
    YY <- Y[which(Y$peakName %in% peak.index),]
    
    
    #gene.id <- YY$GeneName
    
    #FP.rna.atac.L <- lapply(1:length(gene.id),function(u,YY,output.dir){
    
    #print(YY[which(YY$GeneName == gene.id[u]),])
    
    n= dim(YY)[1]
    
    #print(c(YY[which(YY$GeneName == gene.id[u]),4:17]))
    
    FP.rna.atac <- data.frame(measurement=as.numeric(unlist(c(YY[,4:17]))),type=c(rep("rna",8*n),rep("atac",6*n)),pt=c(rep("IL2",4*n),rep("PBS",4*n),rep("IL2",3*n),rep("PBS",3*n)),time=rep("0.75",14*n),gene=rep(gene.id[u],14*n))
    
    #FP.rna.atac
    df.1 <- data.frame(gr=paste0(FP.rna.atac$pt,"_",FP.rna.atac$type),time=FP.rna.atac$time,openness=FP.rna.atac$measurement,cell=FP.rna.atac$type)
    
    #pdf(file.path(output.dir,paste0(gene.id[u],".pdf")))
    
    #plot <- ggplot2.violinplot(data=df.1,xName='gr', yName='openness',
    #                           mainTitle=unique(df.1$time),xtitle="trt",ytitle="openness")
    library(dplyr)
    df.1Summary <- df.1 %>%
        group_by(gr) %>%
        summarize(openness_mean = mean(openness),
                  openness_se = sqrt(var(openness)/length(openness)))
    
    plotViolins <- ggplot(df.1, aes(x = factor(gr), y = openness, fill = factor(gr))) +
        geom_violin(aes(y = openness, fill = factor(cell))) +
        geom_point(aes(y = openness_mean), color = "black", size = 2, data = df.1Summary) +
        geom_errorbar(aes(y = openness_mean, ymin = openness_mean-openness_se, ymax = openness_mean+openness_se),
                      color = "black", width = 0.2, data = df.1Summary) +
        ylim(min(df.1$openness), max(df.1$openness)) +
        theme(legend.position = "none") +
        ggtitle(paste0(gene.id[u],"_At_",unique(df.1$time),"_Time_Point")) + theme(plot.title=element_text(hjust=0.5))+xlab("T")+ylab("M")
    #print(plotViolins)
    
    #dev.off()
    
    plotViolins
    
    #  },YY, output.dir)
    
    #  names(FP.rna.atac.L)=gene.id
    #FP.rna.atac <- do.call(rbind.data.frame,FP.rna.atac.L)
    
    #gene.id <- unique(as.character(FP.rna.atac$gene))
    
    #lapply()
    
    # df.1 <- data.frame(gr=paste0(FP.rna.atac$pt,"_",FP.rna.atac$type),time=FP.rna.atac$time,openness=FP.rna.atac$measurement,cell=FP.rna.atac$type)
    #
    # pdf(file.path(output.dir,paste0(length(unique(as.character(FP.rna.atac$gene))),".pdf")))
    #
    # #plot <- ggplot2.violinplot(data=df.1,xName='gr', yName='openness',
    # #                           mainTitle=unique(df.1$time),xtitle="trt",ytitle="openness")
    # library(dplyr)
    # df.1Summary <- df.1 %>%
    #   group_by(gr) %>%
    #   summarize(openness_mean = mean(openness),
    #             openness_se = sqrt(var(openness)/length(openness)))
    #
    # plotViolins <- ggplot(df.1, aes(x = factor(gr), y = openness, fill = factor(gr))) +
    #   geom_violin(aes(y = openness, fill = factor(cell))) +
    #   geom_point(aes(y = openness_mean), color = "black", size = 2, data = df.1Summary) +
    #   geom_errorbar(aes(y = openness_mean, ymin = openness_mean-openness_se, ymax = openness_mean+openness_se),
    #                 color = "black", width = 0.2, data = df.1Summary) +
    #   ylim(min(y), max(y)) +
    #   theme(legend.position = "none") +
    #   ggtitle(paste0(length(unique(as.character(FP.rna.atac$gene))),"_At_",unique(df.1$time),"_Time_Point")) + theme(plot.title=element_text(hjust=0.5))+xlab("Trt") + ylab("Gene exprepression difference regulated by different chromatin accessibility")
    # print(plotViolins)
    #
    # dev.off()
    # #df.1
    #  FP.rna.atac.L
}

getOneRnaSeq <- function(res.anno.de) {
    fp.il2.16 <- read.csv("/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/RNA-Seq/FP_16vsIL2_16.csv")
    x <- fp.il2.16$ensembl_gene_id
    FP.il2.at.16 <- cbind.data.frame(ENSEMBL = str_sub(x,1,str_locate_all(x,"\\.")[[1]][1,2]-1),fp.il2.16)
    FP.il2.at.16.rna.atac <- merge(res.anno.de[[2]],FP.il2.at.16,by="ENSEMBL")
    FP.il2.at.16.rna.atac
}

mergeRnaWithAtac <- function(fp.il2.16,res.anno.de) {
    x <- fp.il2.16$ensembl_gene_id
    FP.il2.at.16 <- cbind.data.frame(ENSEMBL = str_sub(x,1,str_locate_all(x,"\\.")[[1]][1,2]-1),fp.il2.16)
    FP.il2.at.16.rna.atac <- merge(res.anno.de[[2]],FP.il2.at.16,by="ENSEMBL")
    FP.il2.at.16.rna.atac
}

mergeRnaWithAtac2 <- function(fp.il2.16,res.anno.de) {
    x <- fp.il2.16$ensembl_gene_id
    FP.il2.at.16 <- cbind.data.frame(ENSEMBL = str_sub(x,1,str_locate_all(x,"\\.")[[1]][1,2]-1),fp.il2.16)
    FP.il2.at.16.rna.atac <- merge(res.anno.de,FP.il2.at.16,by="ENSEMBL")
    FP.il2.at.16.rna.atac
}

matchGeneAtc <- function(IL2.PBS.at.0.75.rna.atac) {
    IL2.PBS.at.0.75.rna.atac.2 <- IL2.PBS.at.0.75.rna.atac[order(IL2.PBS.at.0.75.rna.atac$padj.x),]
    IL2.PBS.at.0.75.rna.atac.2.filtered <- filterPeaks2(IL2.PBS.at.0.75.rna.atac.2,0.05,0.58)
    
    if(dim(IL2.PBS.at.0.75.rna.atac.2.filtered)[1]>0){
        IL2.PBS.at.0.75.rna.atac.2.filtered.2 <- IL2.PBS.at.0.75.rna.atac.2.filtered[which(IL2.PBS.at.0.75.rna.atac.2.filtered$log2FoldChange.x>0&IL2.PBS.at.0.75.rna.atac.2.filtered$log2FoldChange.y>0),]
        
        if(dim(IL2.PBS.at.0.75.rna.atac.2.filtered.2)[1]>0){
            IL2.PBS.at.0.75.rna.atac.2.filtered.3 <- IL2.PBS.at.0.75.rna.atac.2.filtered.2[which(IL2.PBS.at.0.75.rna.atac.2.filtered.2$padj.y<0.05),]
            if(dim(IL2.PBS.at.0.75.rna.atac.2.filtered.3)[1]<=0){
                IL2.PBS.at.0.75.rna.atac.2.filtered.3 <- NULL
            }
        }else{
            IL2.PBS.at.0.75.rna.atac.2.filtered.3 <- NULL
        }
    }else
    {
        cat("After applying filtering, nothing is left\nn")
        IL2.PBS.at.0.75.rna.atac.2.filtered.3 <- NULL
    }
    
    IL2.PBS.at.0.75.rna.atac.2.filtered.3
}

makeAtacRnaPlot <- function(FP.il2.at.16.rna.atac) {
    FP.il2.at.16.rna.atac.3 <- matchGeneAtc(FP.il2.at.16.rna.atac)
    
    output.dir <- "~/ATAC-Seq/RNASeq-AtacSeq/geneBased/FP.il2.at.16"
    
    FP.il2.at.16.rna.atac.3.countour.mutiple.gene <-  makeContour4Gene2(FP.il2.at.16.rna.atac.3,atac.index=c(28:33),rna.index=c(46:53),peak.index=FP.il2.at.16.rna.atac.3$peakName,rpkm=FALSE,output.dir,"FP","IL2","16h")
    
    multi.page <- ggarrange(plotlist=FP.il2.at.16.rna.atac.3.countour.mutiple.gene[c(2,4,6,15)],nrow = 2, ncol = 2)
    ggexport(multi.page, filename = file.path(output.dir,"FP.il2.at.16_genes.pdf"))
}

# output.dir <- "~/ATAC-Seq/RNASeq-AtacSeq/geneBased/FP.il2.at.16"
# file.name <- "FP.il2.at.16_all_genes.pdf"

# makeAtacRnaPlot3(FP.il2.at.16.rna.atac,"FP","IL2","16h",output.dir,file.name)

# file.name <- "FP.il2.at.16_all_genes_rpkm.pdf"
# makeAtacRnaPlot3(FP.il2.at.16.rna.atac,"FP","IL2","16h",output.dir,file.name,rpkm=TRUE)


# FP.il2.at.16.rna.atac.with.fc <- data.frame(FP.il2.at.16.rna.atac,FoldChange.x = 2^FP.il2.at.16.rna.atac$log2FoldChange.x)

# FP.il2.at.16.rna.atac.with.fc.same.direction.1<- FP.il2.at.16.rna.atac.with.fc[which(FP.il2.at.16.rna.atac.with.fc$FoldChange.x>1&FP.il2.at.16.rna.atac.with.fc$FoldChange>1),]

# FP.il2.at.16.rna.atac.with.fc.same.direction.2<- FP.il2.at.16.rna.atac.with.fc[which(FP.il2.at.16.rna.atac.with.fc$FoldChange.x<=1&FP.il2.at.16.rna.atac.with.fc$FoldChange<=1),]

# FP.il2.at.16.rna.atac.with.fc.same.direction <- rbind.data.frame(FP.il2.at.16.rna.atac.with.fc.same.direction.1,FP.il2.at.16.rna.atac.with.fc.same.direction.2)

# plot(log(FP.il2.at.16.rna.atac.with.fc.same.direction$FoldChange.x),log(FP.il2.at.16.rna.atac.with.fc.same.direction$FoldChange),xlab = "logFC(chromatin accessibilty)",ylab = "logFC(Gene Expression)")

# cor(log(FP.il2.at.16.rna.atac.with.fc.same.direction$FoldChange.x),log(FP.il2.at.16.rna.atac.with.fc.same.direction$FoldChange))


# FP.il2.at.16.rna.atac.with.fc.same.direction.sig <- FP.il2.at.16.rna.atac.with.fc.same.direction[which(FP.il2.at.16.rna.atac.with.fc.same.direction$pvalue.x<0.01&FP.il2.at.16.rna.atac.with.fc.same.direction$pvalue.y<0.01),]

# plot(log(FP.il2.at.16.rna.atac.with.fc.same.direction.sig$FoldChange.x),log(FP.il2.at.16.rna.atac.with.fc.same.direction.sig$FoldChange),xlab = "logFC(chromatin accessibilty)",ylab = "logFC(Gene Expression)")

# cor(log(FP.il2.at.16.rna.atac.with.fc.same.direction.sig$FoldChange.x),log(FP.il2.at.16.rna.atac.with.fc.same.direction.sig$FoldChange))

# file.name <- "FP.il2.at.16_all_genes_rpkm_same_FC_direction.pdf"
# makeAtacRnaPlot3(FP.il2.at.16.rna.atac.with.fc.same.direction,"FP","IL2","16h",output.dir,file.name,rpkm=TRUE)
# FP.il2.at.16.rna.atac.with.fc.same.direction.2 <- matchGeneAtc(FP.il2.at.16.rna.atac.with.fc.same.direction)

# interested.gene = c("Il2ra","Il2rb","Foxp3","Ap1","Il10","Ctla4","Lag3","Gata3","Tbet","Prdm1","Helios","Granzyme B","Mir223")
# gene.id <- unique(as.character(fp.il2.16$gene_id))

# gene.id[which(gene.id %in% interested.gene)]


# which(FP.il2.at.16.rna.atac.with.fc.same.direction$SYMBOL%in%interested.gene)
# which(FP.il2.at.16.rna.atac$SYMBOL%in%interested.gene)

makeAtacRnaPlot3 <- function(FP.il2.at.16.rna.atac,s1,s2,time.point,output.dir,file.name,select.gene=NULL,rpkm=FALSE) {
    FP.il2.at.16.rna.atac.3 <- matchGeneAtc(FP.il2.at.16.rna.atac)
    
    FP.il2.at.16.rna.atac.3.countour.mutiple.gene <-  makeContour4Gene2(FP.il2.at.16.rna.atac.3,atac.index=c(28:33),rna.index=c(46:53),peak.index=FP.il2.at.16.rna.atac.3$peakName,rpkm=rpkm,output.dir,s1,s2,time.point)
    
    if(!is.null(select.gene))
    {
        multi.page <- ggarrange(plotlist=FP.il2.at.16.rna.atac.3.countour.mutiple.gene[select.gene],nrow = 2, ncol = 2)
        ggexport(multi.page, filename = file.path(output.dir,file.name))
    }else
    {
        multi.page <- ggarrange(plotlist=FP.il2.at.16.rna.atac.3.countour.mutiple.gene,nrow = 2, ncol = 2)
        ggexport(multi.page, filename = file.path(output.dir,file.name))
    }
}

footPrint <- function() {
    
    cmd0 = "pythonw /Users/axy148/anaconda/bin/example_footprint_scores.py"
    
    cmd1 = "python /Users/axy148/anaconda/bin/wellington_footprints.py"
    
    cmd2 = "python /Users/axy148/anaconda/bin/wellington_footprints.py regions /Users/axy148/AiminDropbox/Alejandro_AtacSeq_out_All_Peaks/FP_IL2_al_16.bed reads /Users/axy148/AiminDropBox/Alejandro_AtacSeq_uploaded_2/FP-at-16/align/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.bam outputdir ~/wellington_footprints_out_FP-at-16"
    
    cmd3 = "samtools view /Users/axy148/AiminDropBox/Alejandro_AtacSeq_uploaded_2/FP-at-16/align/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.bam -b /Volumes/Bioinformatics/Aimin_project/ATAC-Seq/Bam_New/2017-10-26-03_S6_R1.PE2SE.nodup_c.bam"
    
    cmd4 = "python /Users/axy148/anaconda/bin/example_footprint_scores.py"
    
    cmd5 = "python /Users/axy148/anaconda/bin/dnase_cut_counter.py /Users/axy148/AiminDropbox/Alejandro_AtacSeq_out_All_Peaks/FP_IL2_al_16.bed /Users/axy148/AiminDropBox/Alejandro_AtacSeq_uploaded_2/FP-at-16/align/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.bam -A TRUE ~/Fp_count.txt"
    
    cmd6 = "python /Users/axy148/anaconda/bin/dnase_average_profile.py /Users/axy148/AiminDropbox/Alejandro_AtacSeq_out_All_Peaks/FP_IL2_al_16.bed /Users/axy148/AiminDropBox/Alejandro_AtacSeq_uploaded_2/FP-at-16/align/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.bam ~/Fp_profile.pdf"
    
    cmd7 = "python /Users/axy148/anaconda/bin/dnase_wig_tracks.py  /Users/axy148/AiminDropbox/Alejandro_AtacSeq_out_All_Peaks/FP_IL2_al_16.bed /Users/axy148/AiminDropBox/Alejandro_AtacSeq_uploaded_2/FP-at-16/align/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.bam FP_FS.bw FP_RS.bw"
    
    cmd8 ="python /Users/axy148/anaconda/bin/dnase_wig_tracks.py  FP_IL2_al_16_test.bed /Users/axy148/AiminDropBox/Alejandro_AtacSeq_uploaded_2/FP-at-16/align/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.bam FP_FS.bw FP_RS.bw"
    
    cmd9 = "python /Users/axy148/anaconda/bin/wellington_footprints.py FP_IL2_al_16_test.bed /Users/axy148/AiminDropBox/Alejandro_AtacSeq_uploaded_2/FP-at-16/align/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.bam ~/wellington_footprints_out_FP-at-16"
    
    cmd10 = "python /Users/axy148/anaconda/bin/wellington_footprints.py FP_IL2_al_16_test2.bed /Users/axy148/AiminDropBox/Alejandro_AtacSeq_uploaded_2/FP-at-16/align/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.bam ~/wellington_footprints_out_FP-at-16-2"
    
    cmd11 = "python /Users/axy148/anaconda/bin/dnase_cut_counter.py FP_IL2_al_16_test2.bed /Users/axy148/AiminDropBox/Alejandro_AtacSeq_uploaded_2/FP-at-16/align/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.bam -A TRUE ~/Fp_count_pp.txt"
    
    
}


# To get normalized RNA-Seq data using Deseq2
#
# fp.il2.16 <- read.csv("/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/RNA-Seq/FP_16vsIL2_16.csv")

# write.table(fp.il2.16[,13:20],file = "~/Dropbox/Alejandro_AtacSeq_RnaSeq/count.fp.il2.txt",quote = FALSE,row.names = FALSE,sep ="\t")

# input.count.file <- "~/Dropbox/Alejandro_AtacSeq_RnaSeq/count.fp.il2.txt"

# res <- DeAnalysis4RnaSeq(input.count.file,output.dir)

DeAnalysis4RnaSeq <- function(input.count.file,output.dir) {
    
    # Metadata
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }
    
    
    countData=data.frame(read.table(input.count.file,header=T,sep='\t'))
    
    counts.metadata <- cbind.data.frame(Sample = colnames(countData),Condition=c(rep("FP",4),rep("IL2",4)),Stage=rep("16",8),stringsAsFactors=F)
    
    counts.metadata.1 <- counts.metadata
    
    counts.metadata.2 <- cbind.data.frame(Sample=counts.metadata.1$Sample,Condition=counts.metadata.1$Condition,Stage=counts.metadata.1$Stage,stringsAsFactors=F)
    colnames(countData) <- counts.metadata.2$Sample
    any(is.na(countData))
    
    counts.metadata.2$group <- factor(paste0(counts.metadata.2$Condition,"_",counts.metadata.2$Stage))
    
    ddsFullCountTable<-DESeqDataSetFromMatrix(
        countData=countData,
        colData=counts.metadata.2,
        design= ~ group)
    
    dds <- estimateSizeFactors(ddsFullCountTable)
    
    dds.normalized <- counts(dds, normalized=TRUE)
    
    dds<-DESeq(dds)
    resultsNames(dds)
    
    ct <- unique(counts.metadata.1$Condition)
    st <- unique(counts.metadata.1$Stage)
    
    ct1 <- paste0(ct[1],"_",st)
    ct2 <- paste0(ct[2],"_",st)
    
    res.DE=results(dds, contrast = c("group",ct1,ct2))
    
    res <- data.frame(res.DE,dds.normalized)
    
    res
    
}

#fp.il2.16.with.normalized <- data.frame(fp.il2.16[,1:12],res[,-c(1:6)])

#fp.il2.16.with.normalized.2 <- mergeRnaWithAtac(fp.il2.16.with.normalized,res.anno.de)

# getLogNormaAtacRnaSeq(FP.il2.at.16.rna.atac.3,atac.index=c(28:33),rna.index=c(46:53),normalize.option="DESeq2")

getLogNormaAtacRnaSeq <- function(IL2.PBS.at.0.75.rna.atac,atac.index,rna.index,normalize.option=c("cpm","rpkm","DESeq2")) {
    
    normalize.option <- match.arg(normalize.option)
    
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }
    
    switch (job.option,
            cpm = {
                countmatrix.rna <- as.matrix(IL2.PBS.at.0.75.rna.atac[,rna.index]) + 1
                cpm <- apply(countmatrix.rna,2, function(x) (x/sum(x))*1000000)
                log.cpm.rna <- log10(cpm)
            },
            rpkm = {
                expression <- data.frame(IL2.PBS.at.0.75.rna.atac[,rna.index]+1,geneLength=IL2.PBS.at.0.75.rna.atac$V5-IL2.PBS.at.0.75.rna.atac$V4)
                geneLength <- expression$geneLength
                expression.rpkm <- apply(expression[,1:8],2,function(u,geneLength){
                    x <- (u/geneLength)*(10^9/sum(u))
                },geneLength)
                log.cpm.rna <-log10(expression.rpkm)
            },
            DESeq2 = {
                
                countData <- IL2.PBS.at.0.75.rna.atac[,rna.index] + 1
                counts.metadata <- cbind.data.frame(Sample = colnames(countData),Condition=c(rep("FP",4),rep("IL2",4)),Stage=rep("16",8),stringsAsFactors=F)
                
                counts.metadata.1 <- counts.metadata
                
                counts.metadata.2 <- cbind.data.frame(Sample=counts.metadata.1$Sample,Condition=counts.metadata.1$Condition,Stage=counts.metadata.1$Stage,stringsAsFactors=F)
                colnames(countData) <- counts.metadata.2$Sample
                any(is.na(countData))
                
                counts.metadata.2$group <- factor(paste0(counts.metadata.2$Condition,"_",counts.metadata.2$Stage))
                
                ddsFullCountTable<-DESeqDataSetFromMatrix(
                    countData=countData,
                    colData=counts.metadata.2,
                    design= ~ group)
                
                dds <- estimateSizeFactors(ddsFullCountTable)
                
                dds.normalized <- counts(dds, normalized=TRUE)
                log.cpm.rna <-log10(dds.normalized)
            }
    )
    
    countmatrix.atac <- as.matrix(IL2.PBS.at.0.75.rna.atac[,atac.index]) + 1
    log.cpm.atac <- log10(countmatrix.atac)
    
    Y <- cbind.data.frame(ENSEMBL=IL2.PBS.at.0.75.rna.atac$ENSEMBL,GeneName =IL2.PBS.at.0.75.rna.atac$gene_id,peakName=IL2.PBS.at.0.75.rna.atac$peakName,log.cpm.rna,log.cpm.atac)
    
    print(Y)
    Y
}

# output.dir <- "~/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks"
# output.file <- "peakAll.bed"

# peak.all <- putAllPeaksTogether(ppr.merged.bed,output.dir,output.file)

putAllPeaksTogether <- function(ppr.merged.bed,output.dir,output.file){
    
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }
    
    peakL <- lapply(ppr.merged.bed[1:6], function(u){
        
        sampleName <- basename(u)
        
        sampleName1 <- gsub("fp-at","FP-at",sampleName)
        
        peak <- read.table(u,header=F)
        
        peak2 <- data.frame(peak,sampleName=rep(sampleName1,dim(peak)[1]))
        
        peak2
    })
    
    peakAll <- do.call(rbind.data.frame,peakL)
    
    colnames(peakAll)[1:3] = c("seqnames","start","end")
    chrOrder<-c(paste("chr",1:19,sep=""),"chrX","chrY")
    
    peakAll$seqnames <- factor(peakAll$seqnames, levels=chrOrder)
    
    peakAll2 <- peakAll[order(peakAll$seqnames,peakAll$start),]
    
    peakID <- c(paste("peak",seq(1:dim(peakAll2)[1]),sep=""))
    
    peakAll3 <- cbind.data.frame(peakAll2,peakID=peakID)
    
    write.table(peakAll3[,1:3],file = file.path(output.dir,output.file),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = F)
    
    peakAll3
    
}

# ppr.merged.bed.file <- "/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll.bed"

getCountUseAllPeaks <- function(input.atac.dir,ppr.merged.bed.file,output.dir) {
    
    if (!dir.exists(dirname(output.dir)))
    {
        dir.create(dirname(output.dir), recursive = TRUE)
    }
    
    file.1 <- list.files(input.atac.dir,pattern="*tn5.tagAlign.gz$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    file.2 <- file.1[-grep("pseudo_reps",file.1)]
    
    #file.2 <- file.1[grep("pooled_pseudo_reps",file.1)]
    
    file.3 <- cbind.data.frame(ID = basename((dirname(dirname(dirname(file.2))))),rep= basename(dirname(file.2)),fileName=basename(file.2),file.2,stringsAsFactors=FALSE)
    
    file.4 <- file.3[-which(file.3$ID %in% c("IL-2vsPBS")),]
    
    file.5 <- cbind.data.frame(ID2 = paste0(file.4$ID,"_",file.4$rep),file.4,stringsAsFactors=FALSE)
    
    file.5[grep("IL-2-at-16_rep1",file.5$ID2),]$ID2 <-  "IL2-at-16_rep1"
    file.5[grep("IL-2-at-16_rep2",file.5$ID2),]$ID2 <-  "IL2-at-16_rep2"
    file.5[grep("IL-2-at-16_rep3",file.5$ID2),]$ID2 <-  "IL2-at-16_rep3"
    
    file.5[grep("IL-2-at-16",file.5$ID),]$ID <-  "IL2-at-16"
    
    file.id <- unique(as.character(file.5$ID2))
    
    cmd.l <- lapply(file.id, function(u,file.5,ppr.merged.bed.file,output.dir){
        
        
        x <- file.5[which(as.character(file.5$ID2) == u),]
        
        cat("\n")
        cat(as.character(x$ID2),"\n")
        
        if(as.character(x$ID) == "IL-2-at-16"){
            
            cmd0 <- paste("gzip -d",x$file.2,"-c >",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),sep = " ")
            
            cmd1 <- paste("python ~/Dropbox/Aimin_project/AtacSeq/inst/deseq2_wrappers/workflow_for_raw_counts/shift_atac_tagalign.py",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),file.path(output.dir,paste0(x$ID2,"slopped")),sep = " ")
            
            cmd2 <- paste("bedtools coverage -counts -a",ppr.merged.bed.file,"-b",file.path(output.dir,paste0(x$ID2,"slopped")),"| cut -f4 >",file.path(output.dir,x$ID2),sep = " ")
            
            cmd <- paste(cmd0,"wait",cmd1,"wait",cmd2,sep = ";")
            
        }else
        {
            cmd0 <- paste("gzip -d",x$file.2,"-c >",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),sep = " ")
            cmd1 <- paste("python ~/Dropbox/Aimin_project/AtacSeq/inst/deseq2_wrappers/workflow_for_raw_counts/shift_atac_tagalign.py",file.path(output.dir,CutStringByNFromEnd(x$fileName,3)),file.path(output.dir,paste0(x$ID2,"slopped")),sep = " ")
            cmd2 <- paste("bedtools coverage -counts -a",ppr.merged.bed.file,"-b",file.path(output.dir,paste0(x$ID2,"slopped")),"| cut -f4 >",file.path(output.dir,x$ID2),sep = " ")
            cmd <- paste(cmd0,"wait",cmd1,"wait",cmd2,sep = ";")
        }
        print(cmd)
        #system(cmd)
    },file.5,ppr.merged.bed.file,output.dir)
}

# getCountUseAllPeaks(input.atac.dir,ppr.merged.bed.file,output.dir)

# cmd <- "ls -lrth /Volumes/Bioinformatics/Aimin_project/ATAC-Seq/CountUseAllPeaks/*rep? | awk '{print $9}'"
# output.dir <- "/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/CountUseAllPeaks"
# makeCountTable3(cmd,output.dir)

# input.count.file <- "/Volumes/Bioinformatics-1/Aimin_project/ATAC-Seq/CountUseAllPeaks/count.txt"
# output.dir <- "/Volumes/Bioinformatics-1/Aimin_project/ATAC-Seq/CountUseAllPeaks"
#
# vst
# res.qc <- qcAllCount(input.count.file,output.dir)
# res.qc.1 <- vst(res.qc)
# plotPCA(res.qc.1,intgroup=c("Sample","Condition","Stage","group"))
# sampleDists <- dist(t(assay(res.qc.1)))
# hc.vst <- hclust(sampleDists)
# plot(hc.vst)
#
# normalized by size factor
# dds <- estimateSizeFactors(res.qc)
# res.qc.dis <- dist(t(as.matrix(counts(dds, normalized=TRUE))))
# hc <- hclust(res.qc.dis)
# plot(hc)
#
qcAllCount <- function(input.count.file,output.dir) {
    
    # Metadata
    if (!dir.exists(dirname(output.dir)))
    {
        dir.create(dirname(output.dir), recursive = TRUE)
    }
    
    countData=data.frame(read.table(input.count.file,header=T,sep='\t'))
    
    counts.metadata <- cbind.data.frame(Sample = colnames(countData),Condition=colnames(countData),Stage=trimGene3(colnames(countData)),stringsAsFactors=F)
    
    counts.metadata$Condition <- trimGene4(counts.metadata$Condition)
    
    counts.metadata.1 <- cbind.data.frame(counts.metadata,ID=paste0(counts.metadata$Condition,".at.",str_sub(counts.metadata$Sample,str_locate(counts.metadata$Sample,"at")[,2]+2,str_length(counts.metadata$Sample))),stringsAsFactors=F)
    
    counts.metadata.1[which(counts.metadata.1$Stage == 0.75),]$Stage <- "early"
    counts.metadata.1[which(counts.metadata.1$Stage == 16),]$Stage <- "late"
    
    counts.metadata.2 <- cbind.data.frame(Sample=counts.metadata.1$ID,Condition=counts.metadata.1$Condition,Stage=counts.metadata.1$Stage,stringsAsFactors=F)
    colnames(countData) <- counts.metadata.2$Sample
    any(is.na(countData))
    
    counts.metadata.2$group <- factor(paste0(counts.metadata.2$Condition,"_",counts.metadata.2$Stage))
    
    counts.metadata.2$Condition <- factor(counts.metadata.2$Condition)
    counts.metadata.2$Stage <- factor(counts.metadata.2$Stage)
    
    ddsFullCountTable<-DESeqDataSetFromMatrix(
        countData=countData,
        colData=counts.metadata.2,
        design= ~ group)
    
    res.qc.1 <- vst(ddsFullCountTable)
    
    sampleDists <- dist(t(assay(res.qc.1)))
    hc.vst <- hclust(sampleDists)
    pdf(file.path(output.dir,"QC_vst_clustering_samples.pdf"))
    plot(hc.vst)
    dev.off()
    
    dds <- estimateSizeFactors(ddsFullCountTable)
    res.qc.dis <- dist(t(as.matrix(counts(dds, normalized=TRUE))))
    hc <- hclust(res.qc.dis)
    
    pdf(file.path(output.dir,"QC_normalized__by_size_factor_clustering_samples.pdf"))
    plot(hc)
    dev.off()
    
    pdf(file.path(output.dir,"PCA_by_Sample.pdf"))
    plotPCA(res.qc.1,intgroup=c("Sample"))
    dev.off()
    
    pdf(file.path(output.dir,"PCA_by_Condition.pdf"))
    plotPCA(res.qc.1,intgroup=c("Condition"))
    dev.off()
    
    pdf(file.path(output.dir,"PCA_by_Stage.pdf"))
    plotPCA(res.qc.1,intgroup=c("Stage"))
    dev.off()
    
    pdf(file.path(output.dir,"PCA_by_group.pdf"))
    plotPCA(res.qc.1,intgroup=c("group"))
    dev.off()
    
    ddsFullCountTable
}

# input.count.file <- "/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/CountUseAllPeaks/count.txt"
# output.dir <- "/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/CountUseAllPeaks"
# res <- doDeAnalysis4AllCount(input.count.file,output.dir)
#
doDeAnalysis4AllCount <- function(input.count.file,output.dir) {
    
    # Metadata
    if (!dir.exists(dirname(output.dir)))
    {
        dir.create(dirname(output.dir), recursive = TRUE)
    }
    
    countData=data.frame(read.table(input.count.file,header=T,sep='\t'))
    
    counts.metadata <- cbind.data.frame(Sample = colnames(countData),Condition=colnames(countData),Stage=trimGene3(colnames(countData)),stringsAsFactors=F)
    
    counts.metadata$Condition <- trimGene4(counts.metadata$Condition)
    
    counts.metadata.1 <- cbind.data.frame(counts.metadata,ID=paste0(counts.metadata$Condition,".at.",str_sub(counts.metadata$Sample,str_locate(counts.metadata$Sample,"at")[,2]+2,str_length(counts.metadata$Sample))),stringsAsFactors=F)
    
    counts.metadata.1[which(counts.metadata.1$Stage == 0.75),]$Stage <- "early"
    counts.metadata.1[which(counts.metadata.1$Stage == 16),]$Stage <- "late"
    
    counts.metadata.2 <- cbind.data.frame(Sample=counts.metadata.1$ID,Condition=counts.metadata.1$Condition,Stage=counts.metadata.1$Stage,stringsAsFactors=F)
    colnames(countData) <- counts.metadata.2$Sample
    any(is.na(countData))
    
    counts.metadata.2$group <- factor(paste0(counts.metadata.2$Condition,"_",counts.metadata.2$Stage))
    
    counts.metadata.2$Condition <- factor(counts.metadata.2$Condition)
    counts.metadata.2$Stage <- factor(counts.metadata.2$Stage)
    
    ddsFullCountTable<-DESeqDataSetFromMatrix(
        countData=countData,
        colData=counts.metadata.2,
        design= ~ group)
    
    dds <- estimateSizeFactors(dds)
    counts(dds, normalized=TRUE)
    
    dds<-DESeq(ddsFullCountTable)
    resultsNames(dds)
    
    raw.counts <- counts(dds)
    normalized.counts <-  counts(dds, normalized=TRUE)
    
    
    # FP late vs FP early
    res.late.vs.early.at.FP=results(dds, contrast = c("group","FP_late","FP_early"))
    
    mergDeWithCount <- function(res.late.vs.early.at.FP,raw.counts,normalized.counts) {
        c1 <- str_sub(res.late.vs.early.at.FP@elementMetadata[2,2],str_locate(res.late.vs.early.at.FP@elementMetadata[2,2],"group")[2]+2,str_locate(res.late.vs.early.at.FP@elementMetadata[2,2],"vs")[1]-2)
        
        c2 <- str_sub(res.late.vs.early.at.FP@elementMetadata[2,2],str_locate(res.late.vs.early.at.FP@elementMetadata[2,2],"vs")[1]+3,str_length(res.late.vs.early.at.FP@elementMetadata[2,2]))
        
        res.late.vs.early.at.FP.with.count <- cbind.data.frame(res.late.vs.early.at.FP,raw.counts[,which(colnames(raw.counts) %in% counts.metadata.2[which(counts.metadata.2$group == c1),]$Sample)],raw.counts[,which(colnames(raw.counts) %in% counts.metadata.2[which(counts.metadata.2$group == c2),]$Sample)],normalized.counts[,which(colnames(normalized.counts) %in% counts.metadata.2[which(counts.metadata.2$group == c1),]$Sample)],normalized.counts[,which(colnames(normalized.counts) %in% counts.metadata.2[which(counts.metadata.2$group == c2),]$Sample)])
        
        res.late.vs.early.at.FP.with.count
    }
    
    res.late.vs.early.at.FP.with.count <- mergDeWithCount(res.late.vs.early.at.FP,raw.counts,normalized.counts)
    
    # IL2 late vs IL2 early
    res.late.vs.early.at.IL2=results(dds, contrast = c("group","IL2_late","IL2_early"))
    res.late.vs.early.at.IL2.with.count <- mergDeWithCount(res.late.vs.early.at.IL2,raw.counts,normalized.counts)
    
    # PBS late vs PBS early
    res.late.vs.early.at.PBS=results(dds, contrast = c("group","PBS_late","PBS_early"))
    res.late.vs.early.at.PBS.with.count <- mergDeWithCount(res.late.vs.early.at.PBS,raw.counts,normalized.counts)
    
    # early
    res.FP.vs.PBS.at.early=results(dds, contrast = c("group","FP_early","PBS_early"))
    res.FP.vs.PBS.at.early.with.count <- mergDeWithCount(res.FP.vs.PBS.at.early,raw.counts,normalized.counts)
    
    res.IL2.vs.PBS.at.early=results(dds, contrast = c("group","IL2_early","PBS_early"))
    res.IL2.vs.PBS.at.early.with.count <- mergDeWithCount(res.IL2.vs.PBS.at.early,raw.counts,normalized.counts)
    
    res.FP.vs.IL2.at.early=results(dds, contrast = c("group","FP_early","IL2_early"))
    res.FP.vs.IL2.at.early.with.count <- mergDeWithCount(res.FP.vs.IL2.at.early,raw.counts,normalized.counts)
    
    # late
    res.FP.vs.PBS.at.late=results(dds, contrast = c("group","FP_late","PBS_late"))
    res.FP.vs.PBS.at.late.with.count <- mergDeWithCount(res.FP.vs.PBS.at.late,raw.counts,normalized.counts)
    
    res.IL2.vs.PBS.at.late=results(dds, contrast = c("group","IL2_late","PBS_late"))
    res.IL2.vs.PBS.at.late.with.count <- mergDeWithCount(res.IL2.vs.PBS.at.late,raw.counts,normalized.counts)
    
    res.FP.vs.IL2.at.late=results(dds, contrast = c("group","FP_late","IL2_late"))
    res.FP.vs.IL2.at.late.with.count <- mergDeWithCount(res.FP.vs.IL2.at.late,raw.counts,normalized.counts)
    
    
    # FP early vs IL2 late
    res.FP_early.vs.IL2_late=results(dds, contrast = c("group","FP_early","IL2_late"))
    res.FP_early.vs.IL2_late.with.count <- mergDeWithCount(res.FP_early.vs.IL2_late,raw.counts,normalized.counts)
    
    # FP late vs IL2 early
    res.FP_late.vs.IL2_early=results(dds, contrast = c("group","FP_late","IL2_early"))
    res.FP_late.vs.IL2_early.with.count <- mergDeWithCount(res.FP_late.vs.IL2_early,raw.counts,normalized.counts)
    
    # FP early vs PBS late
    res.FP_early.vs.PBS_late=results(dds, contrast = c("group","FP_early","PBS_late"))
    res.FP_early.vs.PBS_late.with.count <- mergDeWithCount(res.FP_early.vs.PBS_late,raw.counts,normalized.counts)
    
    # FP late vs PBS early
    res.FP_late.vs.PBS_early=results(dds, contrast = c("group","FP_late","PBS_early"))
    res.FP_late.vs.PBS_early.with.count <- mergDeWithCount(res.FP_late.vs.PBS_early,raw.counts,normalized.counts)
    
    # IL2 early vs PBS late
    res.IL2_early.vs.PBS_late=results(dds, contrast = c("group","IL2_early","PBS_late"))
    res.IL2_early.vs.PBS_late.with.count <- mergDeWithCount(res.IL2_early.vs.PBS_late,raw.counts,normalized.counts)
    
    # IL2 late vs PBS early
    res.IL2_late.vs.PBS_early=results(dds, contrast = c("group","IL2_late","PBS_early"))
    res.IL2_late.vs.PBS_early.with.count <- mergDeWithCount(res.IL2_late.vs.PBS_early,raw.counts,normalized.counts)
    
    ddsFullCountTable2<-DESeqDataSetFromMatrix(
        countData=countData,
        colData=counts.metadata.2,
        design= ~ Condition + Stage)
    
    dds2<-DESeq(ddsFullCountTable2)
    resultsNames(dds2)
    
    res.late.vs.early = results(dds2, contrast = c("Stage","late","early"))
    
    mergDeWithCount.1 <- function(res.late.vs.early,raw.counts,normalized.counts) {
        c1 <- str_sub(res.late.vs.early@elementMetadata[2,2],str_locate(res.late.vs.early@elementMetadata[2,2],"Stage")[2]+2,str_locate(res.late.vs.early@elementMetadata[2,2],"vs")[1]-2)
        
        c2 <- str_sub(res.late.vs.early@elementMetadata[2,2],str_locate(res.late.vs.early@elementMetadata[2,2],"vs")[1]+3,str_length(res.late.vs.early@elementMetadata[2,2]))
        
        res.late.vs.early.with.count <- cbind.data.frame(res.late.vs.early,raw.counts[,which(colnames(raw.counts) %in% counts.metadata.2[which(counts.metadata.2$Stage == c1),]$Sample)],raw.counts[,which(colnames(raw.counts) %in% counts.metadata.2[which(counts.metadata.2$Stage == c2),]$Sample)],normalized.counts[,which(colnames(normalized.counts) %in% counts.metadata.2[which(counts.metadata.2$Stage == c1),]$Sample)],normalized.counts[,which(colnames(normalized.counts) %in% counts.metadata.2[which(counts.metadata.2$Stage == c2),]$Sample)])
        
        res.late.vs.early.with.count
    }
    
    res.late.vs.early.with.count <- mergDeWithCount.1(res.late.vs.early,raw.counts,normalized.counts)
    
    res.FP.vs.PBS = results(dds2, contrast = c("Condition","FP","PBS"))
    
    
    mergDeWithCount.2 <- function(res.FP.vs.PBS,raw.counts,normalized.counts) {
        c1 <- str_sub(res.FP.vs.PBS@elementMetadata[2,2],str_locate(res.FP.vs.PBS@elementMetadata[2,2],"Condition")[2]+2,str_locate(res.FP.vs.PBS@elementMetadata[2,2],"vs")[1]-2)
        
        c2 <- str_sub(res.FP.vs.PBS@elementMetadata[2,2],str_locate(res.FP.vs.PBS@elementMetadata[2,2],"vs")[1]+3,str_length(res.FP.vs.PBS@elementMetadata[2,2]))
        
        res.FP.vs.PBS.with.count <- cbind.data.frame(res.FP.vs.PBS,raw.counts[,which(colnames(raw.counts) %in% counts.metadata.2[which(counts.metadata.2$Condition == c1),]$Sample)],raw.counts[,which(colnames(raw.counts) %in% counts.metadata.2[which(counts.metadata.2$Condition == c2),]$Sample)],normalized.counts[,which(colnames(normalized.counts) %in% counts.metadata.2[which(counts.metadata.2$Condition == c1),]$Sample)],normalized.counts[,which(colnames(normalized.counts) %in% counts.metadata.2[which(counts.metadata.2$Condition == c2),]$Sample)])
        
        res.FP.vs.PBS.with.count
    }
    
    res.FP.vs.PBS.with.count <- mergDeWithCount.2(res.FP.vs.PBS,raw.counts,normalized.counts)
    
    res.IL2.vs.PBS = results(dds2, contrast = c("Condition","IL2","PBS"))
    res.IL2.vs.PBS.with.count <- mergDeWithCount.2(res.IL2.vs.PBS,raw.counts,normalized.counts)
    
    res.FP.vs.IL2 = results(dds2, contrast = c("Condition","FP","IL2"))
    res.FP.vs.IL2.with.count <- mergDeWithCount.2(res.FP.vs.IL2,raw.counts,normalized.counts)
    
    pdf(file.path(output.dir,"DE_pvalue.pdf"))
    par(mfrow = c(3,2))  # 3 rows and 2 columns
    
    hist(res.late.vs.early$pvalue)
    hist(res.FP.vs.PBS$pvalue)
    hist(res.IL2.vs.PBS$pvalue)
    hist(res.FP.vs.IL2$pvalue)
    
    hist(res.late.vs.early.at.FP$pvalue)
    hist(res.late.vs.early.at.IL2$pvalue)
    hist(res.late.vs.early.at.PBS$pvalue)
    
    hist(res.FP.vs.PBS.at.early$pvalue)
    hist(res.IL2.vs.PBS.at.early$pvalue)
    hist(res.FP.vs.IL2.at.early$pvalue)
    
    hist(res.FP.vs.PBS.at.late$pvalue)
    hist(res.IL2.vs.PBS.at.late$pvalue)
    hist(res.FP.vs.IL2.at.late$pvalue)
    
    hist(res.FP_early.vs.IL2_late$pvalue)
    hist(res.FP_late.vs.IL2_early$pvalue)
    hist(res.FP_early.vs.PBS_late$pvalue)
    
    hist(res.FP_late.vs.PBS_early$pvalue)
    hist(res.IL2_early.vs.PBS_late$pvalue)
    hist(res.IL2_late.vs.PBS_early$pvalue)
    
    dev.off()
    
    re <- list(
        raw.counts = raw.counts,
        normalized.counts = normalized.counts,
        dds = dds,
        dds2 = dds2,
        res.late.vs.early = res.late.vs.early.with.count,
        res.FP.vs.PBS =  res.FP.vs.PBS.with.count,
        res.IL2.vs.PBS = res.IL2.vs.PBS.with.count,
        res.FP.vs.IL2 = res.FP.vs.IL2.with.count,
        res.late.vs.early.at.FP = res.late.vs.early.at.FP.with.count,
        res.late.vs.early.at.IL2 = res.late.vs.early.at.IL2.with.count,
        res.late.vs.early.at.PBS = res.late.vs.early.at.PBS.with.count,
        res.FP.vs.PBS.at.early = res.FP.vs.PBS.at.early.with.count,
        res.IL2.vs.PBS.at.early = res.IL2.vs.PBS.at.early.with.count,
        res.FP.vs.IL2.at.early = res.FP.vs.IL2.at.early.with.count,
        res.FP.vs.PBS.at.late =  res.FP.vs.PBS.at.late.with.count,
        res.IL2.vs.PBS.at.late = res.IL2.vs.PBS.at.late.with.count,
        res.FP.vs.IL2.at.late = res.FP.vs.IL2.at.late.with.count,
        res.FP_early.vs.IL2_late = res.FP_early.vs.IL2_late.with.count,
        res.FP_late.vs.IL2_early = res.FP_late.vs.IL2_early.with.count,
        res.FP_early.vs.PBS_late = res.FP_early.vs.PBS_late.with.count,
        res.FP_late.vs.PBS_early =  res.FP_late.vs.PBS_early.with.count,
        res.IL2_early.vs.PBS_late = res.IL2_early.vs.PBS_late.with.count,
        res.IL2_late.vs.PBS_early = res.IL2_late.vs.PBS_early.with.count
    )
    
    re
    
}

#  dir.name="/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/CountUseAllPeaks"
#' input.file.pattern="peakAll.bed"
#' out.dir.name="/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/CountUseAllPeaks"
#' txdb="mm10"
#' DD=5000
#'
#' anno.res <- AtacSeq:::AnntationUsingChipSeeker(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000)

# annoted.all.peak <- orderAndLabelPeak(anno.res[[1]]$peaksAll)

orderAndLabelPeak <- function(peakAll){
    
    chrOrder<-c(paste("chr",1:19,sep=""),"chrX","chrY")
    
    peakAll$seqnames <- factor(peakAll$seqnames, levels=chrOrder)
    
    peakAll2 <- peakAll[order(peakAll$seqnames,peakAll$start),]
    
    peakID <- c(paste("peak",seq(1:dim(peakAll2)[1]),sep=""))
    
    peakAll3 <- cbind.data.frame(peakAll2,peakID=peakID)
    
    peakAll3
    
}


# output.dir <- "/Volumes/Bioinformatics/Aimin_project/ATAC-Seq/CountUseAllPeaks/DeWithAnnotation"
# res.de.anno <- mergDeWithCount.3(annoted.all.peak,res,output.dir)
# save.image(file="inst/extdata/AtacSeq_3_11_2018.Rdata")

mergDeWithCount.3 <- function(annoted.all.peak,res,output.dir){
    
    if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
    
    res.1 <- res[-c(1:4)]
    
    file.name <- names(res.1)
    
    res.de.anno <- lapply(1:length(res.1), function(u,res.1,file.name,annoted.all.peak,output.dir){
        
        temp <- res.1[[u]]
        peakID <- c(paste("peak",seq(1:dim(temp)[1]),sep=""))
        temp <- cbind.data.frame(temp,peakID=peakID)
        
        temp <-  merge.data.frame(annoted.all.peak,temp,by="peakID",sort=FALSE)
        
        x_name <-as.character(file.name[u])
        
        sub("late",16,x_name)
        sub("early",075,x_name)
        
        write.table(temp,file=file.path(output.dir,paste0(x_name,"_annotation_DE.xls")),row.names = FALSE,quote=FALSE,sep="\t")
        
        temp
        
    },res.1,file.name,annoted.all.peak,output.dir)
    
    names(res.de.anno) <- file.name
    
    res.de.anno
    
}

# rna.data.normalized <- getNormalizedRnaSeq(rna.data,c(13:20),normalize.option="DESeq2")
#
getNormalizedRnaSeq <- function(rna.data,rna.index,normalize.option=c("cpm","rpkm","DESeq2")) {
    
    normalize.option <- match.arg(normalize.option)
    
    names(rna.data) <- gsub("IL-2", "IL2",names(rna.data))
    
    comparision.name <- names(rna.data)
    
    rna.list <- lapply(1:length(rna.data),function(u,rna.index,rna.data,comparision.name,normalize.option){
        
        s1 <- str_sub(comparision.name[u],1,str_locate(comparision.name[u],"vs")[,1]-2)
        s2 <- str_sub(comparision.name[u],str_locate(comparision.name[u],"vs")[,2]+2,str_locate(comparision.name[u],"_")[,1]-1)
        tp <- str_sub(comparision.name[u],str_locate(comparision.name[u],"_")[,2]+1,str_locate(comparision.name[u],"hrs")[,1]-2)
        
        switch (normalize.option,
                cpm = {
                    countmatrix.rna <- as.matrix(rna.data[[u]][,rna.index])
                    cpm <- apply(countmatrix.rna,2, function(x) (x/sum(x))*1000000)
                    log.cpm.rna <- log10(cpm+1)
                },
                rpkm = {
                    expression <- data.frame(rna.data[[u]][,rna.index],geneLength=rna.data[[u]]$V5-rna.data[[u]]$V4)
                    geneLength <- expression$geneLength
                    expression.rpkm <- apply(expression[,1:8],2,function(x,geneLength){
                        y <- (x/geneLength)*(10^9/sum(x))
                    },geneLength)
                    log.cpm.rna <-log10(expression.rpkm+1)
                },
                DESeq2 = {
                    countData <- rna.data[[u]][,rna.index]
                    counts.metadata <- cbind.data.frame(Sample = colnames(countData),Condition=c(rep(s1,4),rep(s2,4)),Stage=rep(tp,8),stringsAsFactors=F)
                    counts.metadata.1 <- counts.metadata
                    counts.metadata.2 <- cbind.data.frame(Sample=counts.metadata.1$Sample,Condition=counts.metadata.1$Condition,Stage=counts.metadata.1$Stage,stringsAsFactors=F)
                    colnames(countData) <- counts.metadata.2$Sample
                    any(is.na(countData))
                    
                    counts.metadata.2$group <- factor(paste0(counts.metadata.2$Condition,"_",counts.metadata.2$Stage))
                    
                    ddsFullCountTable<-DESeqDataSetFromMatrix(
                        countData=countData,
                        colData=counts.metadata.2,
                        design= ~ group)
                    
                    dds <- estimateSizeFactors(ddsFullCountTable)
                    dds.normalized <- counts(dds, normalized=TRUE)
                    log.cpm.rna <-log10(dds.normalized+1)
                }
        )
        
        Y <- cbind.data.frame(rna.data[[u]],log.cpm.rna)
        Y
    },rna.index,rna.data,comparision.name,normalize.option)
    
    names(rna.list) <-names(rna.data)
    
    rna.list
    
}


# FP.IL2.16.rna.atac.merge.using.all.peaks <- mergeRnaWithAtac2(fp.il2.16.normalized,res.de.anno$res.FP.vs.IL2.at.late)
#
# write.table(FP.IL2.16.rna.atac.merge.using.all.peaks,file=file.path(output.dir,paste0(x_name,"_annotation_DE.xls")),row.names = FALSE,quote=FALSE,sep="\t")
#
#
# FP.IL2.16.rna.atac.merge.using.all.peaks.3 <- matchGeneAtc(FP.IL2.16.rna.atac.merge.using.all.peaks)
# makeContour4GeneUseAllPeaks(FP.IL2.16.rna.atac.merge.using.all.peaks.3,atac.index=c(31:36),rna.index=c(57:64),peak.index=FP.IL2.16.rna.atac.merge.using.all.peaks.3$peakID,"FP","IL2","16h")

makeContour4GeneUseAllPeaks <- function(IL2.PBS.at.0.75.rna.atac,atac.index,rna.index,peak.index,s1,s2,tp) {
    
    countmatrix.rna <- as.matrix(IL2.PBS.at.0.75.rna.atac[,rna.index])
    log.cpm.rna <- countmatrix.rna
    
    countmatrix.atac <- as.matrix(IL2.PBS.at.0.75.rna.atac[,atac.index]) + 1
    log.cpm.atac <- log10(countmatrix.atac)
    
    
    Y <- cbind.data.frame(ENSEMBL=IL2.PBS.at.0.75.rna.atac$ENSEMBL,GeneName =IL2.PBS.at.0.75.rna.atac$gene_id,peakName=IL2.PBS.at.0.75.rna.atac$peakID,log.cpm.rna,log.cpm.atac)
    
    YY <- Y[which(Y$peakName %in% peak.index),]
    
    gene.id <- YY$GeneName
    
    FP.rna.atac.L <- lapply(1:length(gene.id),function(u,YY,output.dir){
        
        
        n= dim(YY[which(YY$GeneName == gene.id[u]),])[1]
        
        print(c(YY[which(YY$GeneName == gene.id[u]),4:17]))
        
        FP.rna.atac <- data.frame(measurement=as.numeric(unlist(c(YY[which(YY$GeneName == gene.id[u]),4:17]))),type=c(rep("rna",8*n),rep("atac",6*n)),pt=c(rep(s1,4*n),rep(s2,4*n),rep(s1,3*n),rep(s2,3*n)),time=rep(tp,14*n),gene=rep(gene.id[u],14*n))
        
        df.1 <- data.frame(gr=paste0(FP.rna.atac$pt,"_",FP.rna.atac$type),time=FP.rna.atac$time,openness=FP.rna.atac$measurement,cell=FP.rna.atac$type)
        
        library(dplyr)
        library(ggplot2)
        df.1Summary <- df.1 %>%
            group_by(gr) %>%
            summarize(openness_mean = mean(openness),
                      openness_se = sqrt(var(openness)/length(openness)))
        
        plotViolins <- ggplot(df.1, aes(x = factor(gr), y = openness, fill = factor(gr))) +
            geom_violin(aes(y = openness, fill = factor(cell))) +
            geom_point(aes(y = openness_mean), color = "black", size = 2, data = df.1Summary) +
            geom_errorbar(aes(y = openness_mean, ymin = openness_mean-openness_se, ymax = openness_mean+openness_se),
                          color = "black", width = 0.2, data = df.1Summary) +
            ylim(min(df.1$openness), max(df.1$openness)) +
            theme(legend.position = "none") +
            ggtitle(paste0(gene.id[u],"_At_",unique(df.1$time))) + theme(plot.title=element_text(hjust=0.5))+xlab("T")+ylab("M")
        
        plotViolins
        
    },YY)
    
    names(FP.rna.atac.L)=gene.id
    
    FP.rna.atac.L
}

# output.dir <- "~/ATAC-Seq/RNASeq-AtacSeq/geneBased/FP.il2.at.16"
# file.name <- "FP.il2.at.16_all_genes_using_all_peak.pdf"
# makeAtacRnaPlot4(FP.IL2.16.rna.atac.merge.using.all.peaks,"FP","IL2","16h",output.dir,file.name)

makeAtacRnaPlot4 <- function(FP.il2.at.16.rna.atac,s1,s2,time.point,output.dir,file.name,select.gene=NULL) {
    FP.il2.at.16.rna.atac.3 <- matchGeneAtc(FP.il2.at.16.rna.atac)
    
    if(!is.null(FP.il2.at.16.rna.atac.3)){
        
        FP.il2.at.16.rna.atac.3.countour.mutiple.gene <-  makeContour4GeneUseAllPeaks(FP.il2.at.16.rna.atac.3,atac.index=c(31:36),rna.index=c(57:64),peak.index=FP.il2.at.16.rna.atac.3$peakID,s1,s2,time.point)
        
        library(ggpubr)
        
        if(!is.null(select.gene))
        {
            multi.page <- ggarrange(plotlist=FP.il2.at.16.rna.atac.3.countour.mutiple.gene[select.gene],nrow = 2, ncol = 2)
            ggexport(multi.page, filename = file.path(output.dir,file.name))
        }else
        {
            multi.page <- ggarrange(plotlist=FP.il2.at.16.rna.atac.3.countour.mutiple.gene,nrow = 2, ncol = 2)
            ggexport(multi.page, filename = file.path(output.dir,file.name))
        }
    }
    
}

# rna.atac.merge.list <- geneRnaAtacMergeList(res.de.anno, rna.data.normalized)
#
geneRnaAtacMergeList <- function(res.de.anno, rna.data.normalized) {
    
    changeDeDirection <- function(temp){
        temp2 <- data.frame(temp[,1:7],log2FoldChange=(-temp$log2FoldChange),FoldChange=1/temp$FoldChange,temp[,c(10:12,17:20,13:16,25:28,21:24)])
        temp2
    }
    
    rna.data.normalized$`16-vs-1.5_IL2`<- changeDeDirection(rna.data.normalized$`1.5-vs-16_IL2`)
    rna.data.normalized$`16-vs-1.5_PBS`<- changeDeDirection(rna.data.normalized$`1.5-vs-16_PBS`)
    rna.data.normalized$`16-vs-1.5_FP` <- changeDeDirection(rna.data.normalized$`1.5-vs-16_FP`)
    
    
    res.de.anno.1 <- res.de.anno
    
    rna.data.normalized.1 <- rna.data.normalized
    
    names(rna.data.normalized.1) <- gsub("-hrs","",names(rna.data.normalized.1))
    
    names(res.de.anno.1) <- gsub("res.","",names(res.de.anno.1))
    names(res.de.anno.1) <- gsub(".vs.","-vs-",names(res.de.anno.1))
    names(res.de.anno.1) <- gsub("early","1.5",names(res.de.anno.1))
    names(res.de.anno.1) <- gsub("late","16",names(res.de.anno.1))
    names(res.de.anno.1) <- gsub(".at.","_",names(res.de.anno.1))
    rna.atac.match <- data.frame(rna.id = match(names(rna.data.normalized.1),names(rna.data.normalized.1)), atac.id = match(names(rna.data.normalized.1),names(res.de.anno.1)))
    
    rna.atac.match <- rna.atac.match[-which(is.na(rna.atac.match$atac.id)),]
    
    rna.atac.merge.list <- apply(rna.atac.match, 1, function(u,rna.data.normalized,res.de.anno){
        
        if(!is.na(names(rna.data.normalized)[u[1]])&&!is.na(names(res.de.anno)[u[2]]))
        {
            cat(names(rna.data.normalized)[u[1]],"\t",names(res.de.anno)[u[2]],"\n")
            
            temp <- mergeRnaWithAtac2(rna.data.normalized[[u[1]]],res.de.anno[[u[2]]])
            temp
        }
        
    },rna.data.normalized,res.de.anno)
    
    names(rna.atac.merge.list) <- names(rna.data.normalized)[rna.atac.match[,1]]
    rna.atac.merge.list
}

# output.dir <- "~/ATAC-Seq/RNASeqAndAtacSeq/Ccplot"
# makeAtacRnaPlot5(rna.atac.merge.list,output.dir)
#
makeAtacRnaPlot5 <- function(rna.atac.merge.list,...){
    
    if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
    print(output.dir)
    
    fName <- names(rna.atac.merge.list)
    
    cc.L <- lapply(1:length(rna.atac.merge.list),function(u,rna.atac.merge.list,fName){
        
        s1 <- str_sub(fName[u],1,str_locate(fName[u],"vs")[,1]-2)
        s2 <- str_sub(fName[u],str_locate(fName[u],"vs")[,2]+2,str_locate(fName[u],"_")[,1]-1)
        
        if (!is.na(str_locate(fName[u],"hrs")[,1])){
            tp <- str_sub(fName[u],str_locate(fName[u],"_")[,2]+1,str_locate(fName[u],"hrs")[,1]-2)
        }else
        {
            tp <- str_sub(fName[u],str_locate(fName[u],"_")[,2]+1,str_length(fName[u]))
        }
        
        file.name <- paste0(fName[u],"_accessibility_gene.pdf")
        makeAtacRnaPlot4(rna.atac.merge.list[[u]],s1,s2,tp,output.dir,file.name)
        
    },rna.atac.merge.list,fName)
    
}

# input.file.dir <- "~/Dropbox/Alejandro_AtacSeq_uploaded_2"
# input.pattern <- "*nodup.bam$"
# output.file.dir <- "~/Dropbox/Atac-Seq"

# res.tmp <- generateCountFragment(input.file.dir, input.pattern, output.file.dir)

generateCountFragment <- function(input.file.dir, input.pattern, output.file.dir)
{
    
    
    if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
    
    file.1 <- list.files(input.file.dir, pattern = input.pattern,
                         all.files = TRUE, full.names = TRUE, recursive = TRUE,
                         include.dirs = TRUE)
    
    file.2 <- data.frame(ID= basename(dirname(dirname(dirname(file.1)))),rep=basename(dirname(file.1)),file.1)
    
    
    file.3 <- file.2[-which(file.2$ID %in% "IL-2vsPBS"),]
    
    file.3$ID <- gsub("IL-2","IL2",file.3$ID)
    
    
    if (!dir.exists(output.file.dir)) {dir.create(output.file.dir, recursive = TRUE)}
    
    cmd.L <- apply(file.3,1,function(u,file.3,output.file.dir) {
        
        #cat(u[1],"\t",u[2],"\t",u[3],"\n")
        x_name = paste0(u[1],"_",u[2])
        
        cmd = paste("samtools view",u[3],"| awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' >",file.path(output.file.dir,paste0(x_name,"_fragment_length_count.txt")),sep= " ")
        
        #cat(cmd,"\n")
        system(cmd)
        cmd
    },file.3,output.file.dir)
    cmd.L
}

#cmd ="samtools view /home/aiminyan/Dropbox/Alejandro_AtacSeq_uploaded_2/FP-at-0.75/align/rep1/2017-10-27-03_S4_R1.PE2SE.nodup.bam  | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > /home/aiminyan/Dropbox/Atac-Seq/fragment_length_count.txt"

#system(cmd)

# input.file.dir <- "~/Dropbox/Atac-Seq"
# input.pattern <- "*fragment_length_count.txt$"
# output.file.dir <- "~/Dropbox/Atac-Seq/plot"
# file.name <- "FragmentLenDistribution.pdf"

# res.tmp <- plotCountFragment(input.file.dir, input.pattern, output.file.dir,file.name)

plotCountFragment <- function(input.file.dir, input.pattern, output.file.dir,file.name){
    
    if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
    
    file.1 <- list.files(input.file.dir, pattern = input.pattern,
                         all.files = TRUE, full.names = TRUE, recursive = TRUE,
                         include.dirs = TRUE)
    
    file.2 <- data.frame(ID= str_sub(basename(file.1),1,str_locate(basename(file.1),"fragment")[,1]-2),file.1)
    
    cmd.L <- apply(file.2,1,function(u,file.2,output.file.dir) {
        
        x_name = u[1]
        
        count.fl <- read.table(u[2])
        
        min.value <- min(as.numeric(count.fl$V2))
        max.value <- max(as.numeric(count.fl$V2))
        #pdf(file.path(output.file.dir,paste0(x_name,"_fragment_length_distribution.pdf")))
        #fig <- plot(count.fl$V1~count.fl$V2,xaxt="n",ylab = "count",xlab="fragment length",main=x_name,type="h")
        
        fig <- ggplot(count.fl, aes(x=V2, y=V1)) + geom_point() +
            ggtitle(paste0(x_name)) + theme(plot.title=element_text(hjust=0.5))+xlab("Fragment")+ylab("Count") +  scale_x_continuous(breaks=c(0,100,500))
        #axis(1, at = seq(0,max(count.fl$V2), by = 100), las=2)
        #dev.off()
        
    },file.2,output.file.dir)
    
    
    multi.page <- ggarrange(plotlist=cmd.L,nrow = 2, ncol = 2)
    ggexport(multi.page, filename = file.path(output.file.dir,file.name))
    
    #plot(count.fl$V1~count.fl$V2,xaxt="n")
    #axis(1, at = seq(0,max(count.fl$V2), by = 100), las=2)
    
}

#save.image(file="inst/extdata/AtacSeq_3_13_2018_2.Rdata")
# samtools view ATAC_f2q30_sorted.bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > fragment_length_count.txt


generatePhastCons60w.UCSC.mm10 <- function(Mmusculus, i) {
    
    system("wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons.bw -P ~/Dropbox/Atac-Seq/")
    
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(GenomeInfoDb)
    library(rtracklayer)
    
    seqlen <- seqlengths(Mmusculus)
    gr <- GRanges(names(seqlen), IRanges(1, seqlen))
    for(i in 1:length(gr)){
        g <- gr[i]
        d <- import("mm10.60way.phastCons.bw", selection=BigWigSelection(ranges=g), as="RleList")
        saveRDS(d[[as.character(seqnames(g))]], paste0("phastCons60way.UCSC.mm10.", as.character(seqnames(g)), ".rds"))
    }
    
    refgenomeGD <- GenomeDescription(organism=organism(Mmusculus),
                                     common_name=commonName(Mmusculus),
                                     provider=provider(Mmusculus),
                                     provider_version=providerVersion(Mmusculus),
                                     release_date=releaseDate(Mmusculus),
                                     release_name=releaseName(Mmusculus),
                                     seqinfo=Seqinfo(seqnames=seqnames(Mmusculus),
                                                     seqlengths=seqlengths(Mmusculus),
                                                     isCircular=isCircular(Mmusculus),
                                                     genome=releaseName(Mmusculus)))
    
    saveRDS(refgenomeGD, "refgenomeGD.rds")
    unlink("mm10.60way.phastCons.bw")
    
}



qcWithATACseqQC <- function(GenomicScores, ATACseqQC, BSgenome.Mmusculus.UCSC.mm10, GenomeInfoDb, rtracklayer, ChIPpeakAnno, Mmusculus, i, TxDb.Mmusculus.UCSC.mm10.knownGene, MotifDb) {
    
    
    
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Mmusculus.UCSC.mm10")
    biocLite("scater")
    biocLite("scde")
    install.packages("flexmix")
    biocLite("DelayedMatrixStats")
    biocLite("GenomicScores")
    biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
    
    library(GenomicScores)
    
    gsco <- getGScores("phastCons60way.UCSC.mm10")
    gsco
    
    library(ATACseqQC)
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(GenomeInfoDb)
    library(rtracklayer)
    library(ChIPpeakAnno)
    
    
    tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
    ## files will be output into outPath
    outPath <- "splited"
    dir.create(outPath)
    
    seqlev <- "chr1" ## subsample data for quick run
    which <- as(seqinfo(Mmusculus)[seqlev], "GRanges")
    
    bamfile <- "/home/aiminyan/Dropbox/Alejandro_AtacSeq_uploaded_2/FP-at-16/align/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.bam"
    bamfile.labels <- "FP-at-16-rep1"
    fragSize <- fragSizeDist(bamfile, bamfile.labels)
    
    gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE)
    gal1 <- shiftGAlignmentsList(gal)
    shiftedBamfile <- file.path(outPath, "shifted.bam")
    export(gal1, shiftedBamfile)
    
    #library(phastCons100way.UCSC.hg19)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
    ## run program for chromosome 1 only
    txs <- txs[seqnames(txs) %in% "chr1"]
    genome <- Mmusculus
    ## split the reads into NucleosomeFree, mononucleosome,
    ## dinucleosome and trinucleosome.
    objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome,
                                  conservation=gsco)
    
    null <- writeListOfGAlignments(objs, outPath)
    dir(outPath)
    
    
    bamfiles <- file.path(outPath,
                          c("NucleosomeFree.bam",
                            "mononucleosome.bam",
                            "dinucleosome.bam",
                            "trinucleosome.bam"))
    
    cumulativePercentage(bamfiles[1:2], as(seqinfo(Mmusculus)["chr1"], "GRanges"))
    TSS <- promoters(txs, upstream=0, downstream=1)
    TSS <- unique(TSS)
    
    (librarySize <- estLibSize(bamfiles))
    
    NTILE <- 101
    dws <- ups <- 1010
    sigs <- enrichedFragments(gal=objs[c("NucleosomeFree",
                                         "mononucleosome",
                                         "dinucleosome",
                                         "trinucleosome")],
                              TSS=TSS,
                              librarySize=librarySize,
                              seqlev=seqlev,
                              TSS.filter=0.5,
                              n.tile = NTILE,
                              upstream = ups,
                              downstream = dws)
    
    sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
    featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                          zeroAt=.5, n.tile=NTILE)
    
    out <- featureAlignedDistribution(sigs,
                                      reCenterPeaks(TSS, width=ups+dws),
                                      zeroAt=.5, n.tile=NTILE, type="l",
                                      ylab="Averaged coverage")
    
    
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    out <- apply(out, 2, range01)
    matplot(out, type="l", xaxt="n",
            xlab="Position (bp)",
            ylab="Fraction of signal")
    axis(1, at=seq(0, 100, by=10)+1,
         labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
    abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
    
    library(MotifDb)
    CTCF <- query(MotifDb, c("CTCF"))
    CTCF <- as.list(CTCF)
    print(CTCF[[1]], digits=2)
    
    # For chr1
    sigs <- factorFootprints(shiftedBamfile, pfm=CTCF[[1]],
                             genome=genome,
                             min.score="90%", seqlev=seqlev,
                             upstream=100, downstream=100)
    
    
    # For all chr
    sigs <- factorFootprints(shiftedBamfile, pfm=CTCF[[1]],
                             genome=genome,
                             min.score="90%",seqlev = paste0("chr", c(1:19, "X","Y")),
                             upstream=100, downstream=100)
    
    write.table(data.frame(sigs$bindingSites),"~/Dropbox/Atac-Seq/FP-at-16-rep1_CTCT_all_chrs.txt",sep="\t",row.names = F,quote = F,col.names = T)
    
    
    STAT5 <- query(MotifDb, c("STAT5"))
    STAT5 <- as.list(STAT5)
    print(STAT5[[1]], digits=2)
    
    
    # For all chr
    sigs <- factorFootprints(shiftedBamfile, pfm=STAT5[[1]],
                             genome=genome,
                             min.score="90%",seqlev = paste0("chr", c(1:19, "X","Y")),
                             upstream=100, downstream=100)
    
    write.table(data.frame(sigs$bindingSites),"~/Dropbox/Atac-Seq/mm10/FP-at16-rep1_STAT5_all_chrs.txt",sep="\t",row.names = F,quote = F,col.names = T)
    
    
    sigs <- factorFootprints(shiftedBamfile, pfm=STAT5[[1]],
                             genome=genome,
                             min.score="90%",seqlev = paste0("chr", c(1:19, "X","Y")),
                             upstream=100, downstream=100)
    
    write.table(data.frame(sigs$bindingSites),"~/Dropbox/Atac-Seq/mm10/FP-at16-rep1_STAT5_all_chrs.txt",sep="\t",row.names = F,quote = F,col.names = T)
    
    
    FOXP3 <- query(MotifDb, c("FOXP3"))
    FOXP3 <- as.list(FOXP3)
    print(FOXP3[[1]], digits=2)
    
    sigs <- factorFootprints(shiftedBamfile, pfm=FOXP3[[2]],
                             genome=genome,
                             min.score="90%",seqlev = paste0("chr", c(1)),
                             upstream=100, downstream=100)
    
    write.table(data.frame(sigs$bindingSites),"~/Dropbox/Atac-Seq/mm10/FP-at16-rep1_STAT5_all_chrs.txt",sep="\t",row.names = F,quote = F,col.names = T)
    
    
    
    
    featureAlignedHeatmap(sigs$signal,
                          feature.gr=reCenterPeaks(sigs$bindingSites,
                                                   width=200+width(sigs$bindingSites[1])),
                          annoMcols="score",
                          sortBy="score",
                          n.tile=ncol(sigs$signal[[1]]))
    sigs$spearman.correlation
    
    sessionInfo()
}

useCentipede4hg19 <- function() {
    syetem("wget http://cisbp.ccbr.utoronto.ca/data/1.02/DataFiles/PWMs/Files/M6496_1.02.txt -P ~/Dropbox/Atac-Seq/Test/")
    
    ststem("matrix2meme < <(tail -n+2 ~/Dropbox/Atac-Seq/Test/M6496_1.02.txt | cut -f2-) > ~/Dropbox/Atac-Seq/Test/M6496_1.02.meme")
    
    system("cat ~/Dropbox/Atac-Seq/Test/M6496_1.02.meme")
    
    
    syste("wget 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz' -P ~/Dropbox/Atac-Seq/Test/")
    
    system("cd ~/Dropbox/Atac-Seq/mm10;tar -xzvf chromFaMasked.tar.gz")
    
    
    system("cd ~/Dropbox/Atac-Seq/mm10;cat chr*.fa.masked > mm10.fa")
    
    system("wget --no-check-certificate  https://www.encodeproject.org/files/ENCFF001UUQ/@@download/ENCFF001UUQ.bed.gz")
    
    # 4.4 GB
    system("wget --no-check-certificate https://www.encodeproject.org/files/ENCFF000SHS/@@download/ENCFF000SHS.bam")
    
    
    system("cd ~/Dropbox/Atac-Seq/Test;dnase=ENCFF001UUQ.bed.gz;dnase_gt8=ENCFF001UUQ_gt8.narrowPeak.bed;zcat $dnase | awk '{if ($8 > 8) print}' | gzip > $dnase_gt8;genome=hg19.fa;dnase_fasta=ENCFF001UUQ_gt8.fa;bedtools getfasta -fi $genome -bed $dnase_gt8 -fo $dnase_fasta;meme=M6496_1.02.meme;sites=M6496_1.02.fimo.txt.gz;fimo --text --parse-genomic-coord $meme $dnase_fasta | gzip > $sites;zcat $sites | head")
    
    library(Rsamtools)
    library(CENTIPEDE)
    library(CENTIPEDE.tutorial)
    
    
    cen <- centipede_data(
        bam_file = "ENCFF000SHS.bam",
        fimo_file = "M6496_1.02.fimo.txt",
        pvalue = 1e-4,
        flank_size = 100
    )
    head(cen$regions)
    
    rowSums(cen$mat)[1:10]
    
    library(CENTIPEDE)
    
    
    plot(cen$mat[2096,], xlab = "Position", ylab = "Read Start Sites", type = "h",
         col = rep(c("blue", "red"), each = 213))
    abline(v = c(100, 113, 313, 326) + 0.5, lty = 2)
    abline(v = 213 + 0.5)
    
    
    fit <- fitCentipede(
        Xlist = list(DNase = cen$mat),
        Y = as.matrix(data.frame(
            Intercept = rep(1, nrow(cen$mat))
        ))
    )
    
    sum(fit$PostPr == 1)
    imageCutSitesCombined(cen$mat[fit$PostPr == 1,])
    
    plotProfile(fit$LambdaParList[[1]], Mlen = 13)
    
    cons <- read_bedGraph('ENCFF001UUQ_gt8_phastCons.bed.gz')
    flank_size <- 100L
    sites <- GRanges(
        seqnames = Rle(cen$regions$sequence.name),
        ranges = IRanges(
            start = cen$regions$start,
            end = cen$regions$stop
        ),
        strand = Rle(cen$regions$strand)
    )
    sites <- resize(sites, width(sites) - flank_size, fix = "end")
    sites <- resize(sites, width(sites) - flank_size, fix = "start")
    
    
    xs <- findOverlaps(sites, cons)
    site_cons <- sapply(1:length(sites), function(i) {
        # Conservation scores for each positions in a PWM match.
        ys <- cons[subjectHits(xs[queryHits(xs) == i])]
        vals <- rep(ys$score, width(ys))
        idx <- seq(
            from = start(sites[i]) - min(start(ys)) + 1,
            length.out = width(sites[i])
        )
        vals <- vals[idx]
        mean(vals)
    })
    
    hist(site_cons)
    
    plot(site_cons, log10(rowSums(cen$mat) + 1),
         ylab = "Log10 Read Starts",
         xlab = "phastCons across 100 vertebrates")
    
    
    
    fit2 <- fitCentipede(
        Xlist = list(DNase = cen$mat),
        Y = as.matrix(data.frame(
            Intercept = rep(1, nrow(cen$mat)),
            Conservation = site_cons
        ))
    )
    plotProfile(fit2$LambdaParList[[1]], Mlen = 13)
    all.equal(fit2$PostPr == 1, fit$PostPr == 1)
    
    
    range(fit2$PostPr - fit$PostPr)
    
    plot(sort(fit2$PostPr - fit$PostPr), ylab = "delta PostPr")
    
    
    idx <- site_cons > 0.999
    fit3 <- fitCentipede(
        Xlist = list(DNase = cen$mat[idx, ]),
        Y = as.matrix(data.frame(
            Intercept = rep(1, nrow(cen$mat[idx, ])),
            Conservation = site_cons[idx]
        ))
    )
    
    plotProfile(fit3$LambdaParList[[1]], Mlen = 13)
}

useCentipede4mm10 <- function() {
    system("wget http://cisbp.ccbr.utoronto.ca/data/1.02/DataFiles/PWMs/Files/M6496_1.02.txt -P ~/Dropbox/Atac-Seq/Test/;matrix2meme < <(tail -n+2 ~/Dropbox/Atac-Seq/Test/M6496_1.02.txt | cut -f2-) > ~/Dropbox/Atac-Seq/Test/M6496_1.02.meme;cat ~/Dropbox/Atac-Seq/Test/M6496_1.02.meme")
    
    
    system("wget 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFaMasked.tar.gz' -P ~/Dropbox/Atac-Seq/mm10/")
    
    system("cd ~/Dropbox/Atac-Seq/mm10;tar -xzvf chromFaMasked.tar.gz;gunzip -c chr*.fa.masked > hg19.fa;wget --no-check-certificate https://www.encodeproject.org/files/ENCFF001UUQ/@@download/ENCFF001UUQ.bed.gz")
    
    # 4.4 GB
    system("wget --no-check-certificate https://www.encodeproject.org/files/ENCFF000SHS/@@download/ENCFF000SHS.bam")
    
    matrices.mouse.stat5 <- query(query(MotifDb, 'Mmusculus'),'STAT5')
    matrices.mouse.stat5
    matrix.output.file ="~/Dropbox/Atac-Seq/mm10/stat5.meme"
    meme.text = export (matrices.mouse.stat5, matrix.output.file,'meme')
    #query(MotifDb, 'STAT5'')
    
    
    system("cd ~/Dropbox/Atac-Seq/mm10;dnase=~/Dropbox/Alejandro_AtacSeq_uploaded_2/FP-at-16/peak/macs2/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.tn5.pf.narrowPeak.gz;dnase_gt8=FP-at-16-rep1_gt8.narrowPeak.bed;zcat $dnase | awk '{if ($8 > 8) print}' | gzip > $dnase_gt8;genome=mm10.fa;dnase_fasta=FP-at-16-rep1_gt8.fa;bedtools getfasta -fi $genome -bed $dnase_gt8 -fo $dnase_fasta;meme=~/Dropbox/Atac-Seq/mm10/stat5.meme;sites=stat5.fimo.txt;fimo --text --parse-genomic-coord $meme $dnase_fasta > $sites;head $sites")
    
    
    matrices.mouse.ctcf <- query(query(MotifDb, 'Mmusculus'),'CTCF')
    matrices.mouse.ctcf
    matrix.output.file ="~/Dropbox/Atac-Seq/mm10/ctcf.meme"
    meme.text = export (matrices.mouse.ctcf, matrix.output.file,'meme')
    #query(MotifDb, 'STAT5'')
    
    system("cd ~/Dropbox/Atac-Seq/mm10;dnase=~/Dropbox/Alejandro_AtacSeq_uploaded_2/FP-at-16/peak/macs2/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.tn5.pf.narrowPeak.gz;dnase_gt8=FP-at-16-rep1_gt8.narrowPeak.bed;zcat $dnase | awk '{if ($8 > 8) print}' | gzip > $dnase_gt8;genome=mm10.fa;dnase_fasta=FP-at-16-rep1_gt8.fa;bedtools getfasta -fi $genome -bed $dnase_gt8 -fo $dnase_fasta;meme=~/Dropbox/Atac-Seq/mm10/ctcf.meme;sites=ctcf.fimo.txt;fimo --text --parse-genomic-coord $meme $dnase_fasta > $sites;head $sites")
    
    
    library(Rsamtools)
    library(CENTIPEDE)
    library(CENTIPEDE.tutorial)
    
    bamfile <- "/home/aiminyan/Dropbox/Alejandro_AtacSeq_uploaded_2/FP-at-16/align/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.bam"
    bamfile.labels <- "FP-at-16-rep1"
    #fragSize <- fragSizeDist(bamfile, bamfile.labels)
    
    gal <- readBamFile(bamfile, tag=tags,asMates=TRUE)
    gal1 <- shiftGAlignmentsList(gal)
    outPath <- "~/Dropbox/Atac-Seq/mm10"
    shiftedBamfile <- file.path(outPath, "FP-at-16-rep1_shifted.bam")
    export(gal1, shiftedBamfile)
    
    
    
    
    cen <- centipede_data(
        bam_file = "~/Dropbox/Atac-Seq/mm10/FP-at-16-rep1_shifted.bam",
        fimo_file = "~/Dropbox/Atac-Seq/mm10/ctcf.fimo.txt",
        pvalue = 1e-4,
        flank_size = 100
    )
    head(cen$regions)
    
    rowSums(cen$mat)[1:10]
    
    library(CENTIPEDE)
    
    
    plot(cen$mat[2096,], xlab = "Position", ylab = "Read Start Sites", type = "h",
         col = rep(c("blue", "red"), each = 213))
    abline(v = c(100, 113, 313, 326) + 0.5, lty = 2)
    abline(v = 213 + 0.5)
    
    
    fit <- fitCentipede(
        Xlist = list(DNase = cen$mat),
        Y = as.matrix(data.frame(
            Intercept = rep(1, nrow(cen$mat))
        ))
    )
    
    sum(fit$PostPr == 1 )
    imageCutSitesCombined(cen$mat[fit$PostPr ==1,])
    
    plotProfile(fit$LambdaParList[[1]], Mlen = 20)
    
    cons <- read_bedGraph('~/Dropbox/Atac-Seq/mm10/FP-at-16-rep1_gt8_phastCons.bed.gz')
    flank_size <- 50L
    sites <- GRanges(
        seqnames = Rle(cen$regions$sequence.name),
        ranges = IRanges(
            start = cen$regions$start,
            end = cen$regions$stop
        ),
        strand = Rle(cen$regions$strand)
    )
    sites <- resize(sites, width(sites) - flank_size, fix = "end")
    sites <- resize(sites, width(sites) - flank_size, fix = "start")
    
    
    xs <- findOverlaps(sites, cons)
    site_cons <- sapply(1:length(sites), function(i) {
        # Conservation scores for each positions in a PWM match.
        ys <- cons[subjectHits(xs[queryHits(xs) == i])]
        vals <- rep(ys$score, width(ys))
        idx <- seq(
            from = start(sites[i]) - min(start(ys)) + 1,
            length.out = width(sites[i])
        )
        vals <- vals[idx]
        mean(vals)
    })
    
    hist(site_cons)
    
    plot(site_cons, log10(rowSums(cen$mat) + 1),
         ylab = "Log10 Read Starts",
         xlab = "phastCons across 100 vertebrates")
    
    
    
    fit2 <- fitCentipede(
        Xlist = list(DNase = cen$mat),
        Y = as.matrix(data.frame(
            Intercept = rep(1, nrow(cen$mat)),
            Conservation = site_cons
        ))
    )
    plotProfile(fit2$LambdaParList[[1]], Mlen = 13)
    all.equal(fit2$PostPr == 1, fit$PostPr == 1)
    
    
    range(fit2$PostPr - fit$PostPr)
    
    plot(sort(fit2$PostPr - fit$PostPr), ylab = "delta PostPr")
    
    
    idx <- site_cons > 0.999
    fit3 <- fitCentipede(
        Xlist = list(DNase = cen$mat[idx, ]),
        Y = as.matrix(data.frame(
            Intercept = rep(1, nrow(cen$mat[idx, ])),
            Conservation = site_cons[idx]
        ))
    )
    
    plotProfile(fit3$LambdaParList[[1]], Mlen = 13)
}


# output.dir <- "~/Dropbox/Atac-Seq/homer"
# peaks.selected <- generateHOMERpeakFiles(rna.atac.merge.list,output.dir)
#

generateHOMERpeakFiles <- function(rna.atac.merge.list,output.dir){
    
    if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
    
    
    peaks.selected <- lapply(1:length(rna.atac.merge.list), function(u,rna.atac.merge.list){
        
        FP.vs.IL2.16.hrs.4.homer <- filterPeaks2(rna.atac.merge.list[[u]],0.05,0.58)
        
        x_name <- names(rna.atac.merge.list)[u]
        
        cat(x_name,"\t",dim(FP.vs.IL2.16.hrs.4.homer)[1],"\n")
        
        #x_name <- names(rna.atac.merge.list)[u]
        
        if(dim(FP.vs.IL2.16.hrs.4.homer)[1] > 0){
            
            write.table(FP.vs.IL2.16.hrs.4.homer[,c(2,3,4,5,42)],file=file.path(output.dir,paste0(x_name,".txt")),row.names = FALSE,col.names=FALSE,quote=FALSE,sep="\t")
            
        }
        
        FP.vs.IL2.16.hrs.4.homer
        
    },rna.atac.merge.list)
    
    names(peaks.selected) <- names(rna.atac.merge.list)
    peaks.selected
}

# input.peak.dir <- "~/Dropbox/Atac-Seq/homer"
# input.pattern <- "*.txt$"
# output.dir <- "~/Dropbox/Atac-Seq/homer"
#
# useHomer(input.peak.dir,input.pattern,output.dir,"mm10")

useHomer <- function(input.peak.dir,input.pattern,output.dir,genome){
    
    if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
    
    file.1 <- list.files(input.peak.dir, pattern = input.pattern,
                         all.files = TRUE, full.names = TRUE, recursive = FALSE,
                         include.dirs = TRUE)
    
    null <- lapply(file.1,function(u){
        
        x_name <- tools::file_path_sans_ext(basename(u))
        
        cmd = paste("findMotifsGenome.pl",u,genome,file.path(output.dir,paste0("Analysis_",x_name)),"-size 200",sep=" ")
        cat(cmd,"\n")
        system(cmd)
    })
}

useGoldmine <- function(){
    #
    source("http://bioconductor.org/biocLite.R")
    biocLite(c("GenomicRanges","IRanges","devtools"))
    library(devtools)
    install_github("jeffbhasin/goldmine")
    library(goldmine)
    
    csvpath <- system.file("extdata", "dmrs.csv", package = "goldmine")
    query <- read.csv(csvpath)
    
    head(query)
    cachedir <- "gbcache"
    genome <- "hg19"
    gm <- goldmine(query=query,genome=genome,cachedir=cachedir)
    
    summary(gm)
    
    nrow(gm$context)
    colnames(gm$context)
    
    nrow(gm$genes)
    colnames(gm$genes)
    
    genes <- getGenes("refseq",genome=genome,cachedir=cachedir)
    genes <- genes[str_detect(genes$isoform.id,"NM"),]
    gm <- goldmine(query=query,genes=genes,genome=genome,cachedir=cachedir)
    nrow(gm$genes)
    
    features <- getFeatures(tables=c("wgEncodeRegDnaseClusteredV3",
                                     "wgEncodeRegTfbsClusteredV3"),
                            genome=genome, cachedir=cachedir)
    
    features <- c(features,getCpgFeatures(genome=genome,cachedir=cachedir))
    summary(features)
    
    cd8_states <- fread("http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E047_15_coreMarks_segments.bed")
    
    cd4_states <- fread("http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E039_15_coreMarks_segments.bed")
    
    enh <- rbind(cd8_states,cd4_states)
    enh <- enh[V4=="E7",]
    setnames(enh,c("chr","start","end","state"))
    enh.gr <- reduce(makeGRanges(enh))
    features$enhancers <- enh.gr
    
    gm <- goldmine(query=query,genes=getGenes("ensembl",genome=genome,cachedir)
                   ,features=features,genome=genome,cachedir=cachedir)
    
    colnames(gm$context)
    
    summary(gm$features)
    
    colnames(gm$features$wgEncodeRegTfbsClusteredV3)
    gmWrite(gm, path="gm_csv")
    
    gencon <- gm$context[,list(count=length(chr)),by=c("pattern","call")]
    gencon$call <- factor(gencon$call,levels=c("promoter","exon","intron","3' end","intergenic"))
    gencon <- gencon[,list(call=call,count=count,total=sum(count),
                           fraction=count/sum(count)),by="pattern"]
    gencon
    
    ggplot(gencon,aes(x=call,y=fraction,fill=pattern)) + geom_bar(stat="identity",
                                                                  position="dodge") + ggnice() + labs(title="Genomic Context of DMRs")
    
    
    featcon <- gm$context[,list(CPGisland=sum(cpgIsland_per>0)/length(chr),
                                CPGshore=sum(cpgShore_per>0)/length(chr),
                                CPGshelf=sum(cpgShelf_per>0)/length(chr),
                                TFBS=sum(wgEncodeRegTfbsClusteredV3_per>0)/length(chr),
                                DNaseI=sum(wgEncodeRegDnaseClusteredV3_per>0)/length(chr),
                                Enhancers=sum(enhancers_per>0)/length(chr)),
                          by=c("pattern")]
    featcon <- melt(featcon,id.vars=c("pattern"))
    setnames(featcon,c("variable","value"),c("call","percent"))
    featcon
    ggplot(featcon,aes(x=call,y=percent,fill=pattern)) + geom_bar(stat="identity",
                                                                  position="dodge") + ggnice() + labs(title="Feature Context of DMRs")
    
    etf <- getFeatures("wgEncodeRegTfbsClusteredV3",
                       genome=genome,cachedir=cachedir)[[1]]
    length(etf)
    length(unique(etf$name))
    
    con <- goldmine(etf,genome=genome,cachedir=cachedir,contextonly=TRUE)
    
    
    agg <- con[,list(n=length(chr)),by=c("name","call")]
    agg <- agg[,list(call=call,n=n,frac=n/sum(n)),by="name"]
    head(agg)
    
    levs <- agg[call=="promoter",][order(frac,decreasing=T),]$name
    agg$name <- factor(agg$name,levels=levs)
    agg <- agg[order(agg$name),]
    agg$call <- factor(agg$call,
                       levels=c("promoter","exon","intron","3' end","intergenic"))
    agg <- agg[order(agg$call),]
    ggplot(agg,aes(x=name,y=frac,fill=call)) + geom_bar(stat="identity",
                                                        width=1.3) + ggnice() + coord_flip() + theme(
                                                            axis.ticks.y=element_blank(),
                                                            axis.text.y=element_blank()) + scale_fill_manual(values=
                                                                                                                 c("promoter"="#e41a1c","exon"="#4daf4a","intron"="#377eb8",
                                                                                                                   "3' end"="#ff7f00","intergenic"="#984ea3")) + theme(legend.position="bottom")
    
    tab <- getUCSCTable(table="wgEncodeAwgTfbsUwHct116CtcfUniPk",
                        genome=genome, cachedir=cachedir)
}


useGoldmine4mm10 <- function(){
    #
    # source("http://bioconductor.org/biocLite.R")
    # biocLite(c("GenomicRanges","IRanges","devtools"))
    # library(devtools)
    # install_github("jeffbhasin/goldmine")
    # library(goldmine)
    
    csvpath <- system.file("extdata", "dmrs.csv", package = "goldmine")
    query <- read.csv(csvpath)
    
    query <- peak.all
    
    colnames(query)[1] = "chr"
    
    head(query)
    cachedir <- "gbcache"
    genome <- "mm10"
    gm <- goldmine(query=query,genome=genome,cachedir=cachedir)
    
    summary(gm)
    
    nrow(gm$context)
    colnames(gm$context)
    
    nrow(gm$genes)
    colnames(gm$genes)
    
    genes <- getGenes("refseq",genome=genome,cachedir=cachedir)
    genes <- genes[str_detect(genes$isoform.id,"NM"),]
    gm <- goldmine(query=query,genes=genes,genome=genome,cachedir=cachedir)
    nrow(gm$genes)
    
    features <- getFeatures(tables=c("wgEncodeRegDnaseClusteredV3",
                                     "wgEncodeRegTfbsClusteredV3"),
                            genome=genome, cachedir=cachedir)
    
    features <- c(features,getCpgFeatures(genome=genome,cachedir=cachedir))
    summary(features)
    
    cd8_states <- fread("http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E047_15_coreMarks_segments.bed")
    
    cd4_states <- fread("http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E039_15_coreMarks_segments.bed")
    
    enh <- rbind(cd8_states,cd4_states)
    enh <- enh[V4=="E7",]
    setnames(enh,c("chr","start","end","state"))
    enh.gr <- reduce(makeGRanges(enh))
    features$enhancers <- enh.gr
    
    gm <- goldmine(query=query,genes=getGenes("ensembl",genome=genome,cachedir)
                   ,features=features,genome=genome,cachedir=cachedir)
    
    colnames(gm$context)
    
    summary(gm$features)
    
    colnames(gm$features$wgEncodeRegTfbsClusteredV3)
    gmWrite(gm, path="gm_csv")
    
    gencon <- gm$context[,list(count=length(chr)),by=c("pattern","call")]
    gencon$call <- factor(gencon$call,levels=c("promoter","exon","intron","3' end","intergenic"))
    gencon <- gencon[,list(call=call,count=count,total=sum(count),
                           fraction=count/sum(count)),by="pattern"]
    gencon
    
    ggplot(gencon,aes(x=call,y=fraction,fill=pattern)) + geom_bar(stat="identity",
                                                                  position="dodge") + ggnice() + labs(title="Genomic Context of DMRs")
    
    
    featcon <- gm$context[,list(CPGisland=sum(cpgIsland_per>0)/length(chr),
                                CPGshore=sum(cpgShore_per>0)/length(chr),
                                CPGshelf=sum(cpgShelf_per>0)/length(chr),
                                TFBS=sum(wgEncodeRegTfbsClusteredV3_per>0)/length(chr),
                                DNaseI=sum(wgEncodeRegDnaseClusteredV3_per>0)/length(chr),
                                Enhancers=sum(enhancers_per>0)/length(chr)),
                          by=c("pattern")]
    featcon <- melt(featcon,id.vars=c("pattern"))
    setnames(featcon,c("variable","value"),c("call","percent"))
    featcon
    ggplot(featcon,aes(x=call,y=percent,fill=pattern)) + geom_bar(stat="identity",
                                                                  position="dodge") + ggnice() + labs(title="Feature Context of DMRs")
    
    etf <- getFeatures("wgEncodeRegTfbsClusteredV3",
                       genome=genome,cachedir=cachedir)[[1]]
    length(etf)
    length(unique(etf$name))
    
    con <- goldmine(etf,genome=genome,cachedir=cachedir,contextonly=TRUE)
    
    
    agg <- con[,list(n=length(chr)),by=c("name","call")]
    agg <- agg[,list(call=call,n=n,frac=n/sum(n)),by="name"]
    head(agg)
    
    levs <- agg[call=="promoter",][order(frac,decreasing=T),]$name
    agg$name <- factor(agg$name,levels=levs)
    agg <- agg[order(agg$name),]
    agg$call <- factor(agg$call,
                       levels=c("promoter","exon","intron","3' end","intergenic"))
    agg <- agg[order(agg$call),]
    ggplot(agg,aes(x=name,y=frac,fill=call)) + geom_bar(stat="identity",
                                                        width=1.3) + ggnice() + coord_flip() + theme(
                                                            axis.ticks.y=element_blank(),
                                                            axis.text.y=element_blank()) + scale_fill_manual(values=
                                                                                                                 c("promoter"="#e41a1c","exon"="#4daf4a","intron"="#377eb8",
                                                                                                                   "3' end"="#ff7f00","intergenic"="#984ea3")) + theme(legend.position="bottom")
    
    tab <- getUCSCTable(table="wgEncodeAwgTfbsUwHct116CtcfUniPk",
                        genome=genome, cachedir=cachedir)
}

# input.bed.file <- "~/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_2.bed"
# input.bed.file <- "/home/aiminyan/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_2_in_mm9.bed"
# useChromHMM(input.bed.file)

useChromHMM <- function(input.bed.file) {
    
    # write.table(peak.all,file ="~/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_2.bed",col.names = F,row.names = F, sep="\t",quote = F)
    
    system("/home/aiminyan/kentUtils/bin/linux.x86_64/liftOver /home/aiminyan/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_2.bed ~/Downloads/mm10ToMm9.over.chain.gz /home/aiminyan/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_2_in_mm9.bed /home/aiminyan/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_in_unMapped.bed")
    
    peak.all.in.mm9 <- read.table(file ="/home/aiminyan/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_2_in_mm9.bed")
    
    colnames(peak.all.in.mm9) <- colnames(peak.all)
    
    #
    #   mkdir ~/Dropbox/Test_chromhmm
    #   /home/aiminyan/miniconda3/share/chromhmm-1.14-/download_chromhmm_data.sh ~/Dropbox/Test_chromhmm
    #
    #   cd ~/Dropbox/Test_chromhmm
    #   java -mx1200M -jar /home/aiminyan/miniconda3/share/chromhmm-1.14-/ChromHMM.jar LearnModel SAMPLEDATA_HG18 OUTPUTSAMPLE 10 hg18
    #   ls -lrth OUTPUTSAMPLE/
    #     total 4.9M
    #   -rw-rw-r-- 1 aiminyan aiminyan 1.9K Mar 23 14:52 transitions_10.txt
    #   -rw-rw-r-- 1 aiminyan aiminyan 2.2K Mar 23 14:52 emissions_10.txt
    #   -rw-rw-r-- 1 aiminyan aiminyan  70K Mar 23 14:52 emissions_10.svg
    #   -rw-rw-r-- 1 aiminyan aiminyan  13K Mar 23 14:52 emissions_10.png
    #   -rw-rw-r-- 1 aiminyan aiminyan  12K Mar 23 14:52 transitions_10.png
    #   -rw-rw-r-- 1 aiminyan aiminyan  53K Mar 23 14:52 transitions_10.svg
    #   -rw-rw-r-- 1 aiminyan aiminyan  14K Mar 23 14:52 model_10.txt
    #   -rw-rw-r-- 1 aiminyan aiminyan 681K Mar 23 14:52 K562_10_segments.bed
    #   -rw-rw-r-- 1 aiminyan aiminyan 599K Mar 23 14:52 GM12878_10_segments.bed
    #   -rw-rw-r-- 1 aiminyan aiminyan 1.3M Mar 23 14:52 GM12878_10_dense.bed
    #   -rw-rw-r-- 1 aiminyan aiminyan 304K Mar 23 14:52 GM12878_10_expanded.bed
    #   -rw-rw-r-- 1 aiminyan aiminyan 1.5M Mar 23 14:52 K562_10_dense.bed
    #   -rw-rw-r-- 1 aiminyan aiminyan 1.5K Mar 23 14:52 webpage_10.html
    #   -rw-rw-r-- 1 aiminyan aiminyan 343K Mar 23 14:52 K562_10_expanded.bed
    #
    #   cd ~/enhancer-snakemake-demo/
    #
    #   snakemake -npr
    #   /home/aiminyan/miniconda3/share/chromhmm-1.14-/COORDS/mm9/RefSeqTSS.mm9.bed.gz compare/links/RefSeqTSS.mm9.bed.gz
    #   /home/aiminyan/miniconda3/share/chromhmm-1.14-/COORDS/mm9/RefSeqTES.mm9.bed.gz compare/links/RefSeqTES.mm9.bed.gz
    #   Building DAG of jobs...
    #   MissingInputException in line 73 of /home/aiminyan/enhancer-snakemake-demo/Snakefile:
    #     Missing input files for rule to_compare:
    #     compare/links/RefSeqTES.mm9.bed.gz
    #   compare/links/RefSeqTSS.mm9.bed.gz
    #
    #   ln -sfn ~/Dropbox/Test_chromhmm/COORDS/mm9/RefSeqTSS.mm9.bed.gz RefSeqTSS.mm9.bed.gz
    #   ln -sfn ~/Dropbox/Test_chromhmm/COORDS/mm9/RefSeqTES.mm9.bed.gz RefSeqTES.mm9.bed.gz
    #
    #
    #
    #   sort -k 1,1 -k2,2n /home/aiminyan/enhancer-snakemake-demo/output/5-state/embryonic-liver_5_segments.bed > /home/aiminyan/enhancer-snakemake-demo/output/5-state/embryonic-liver_5_segments_sorted.bed
    #
    #   sort -k 1,1 -k2,2n /H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll.bed > /H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_sorted.bed
    #
    #   bedtools fisher -a /home/aiminyan/enhancer-snakemake-demo/output/5-state/embryonic-liver_5_segments_sorted.bed -b ~/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_sorted.bed -g mm9
    #
    #   bedtools fisher -a ~/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_sorted.bed -b /home/aiminyan/enhancer-snakemake-demo/output/5-state/embryonic-liver_5_segments_sorted.bed -g mm9
    #
    #   system("intersectBed -a /home/aiminyan/enhancer-snakemake-demo/output/5-state/embryonic-liver_5_segments.bed -b ~/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_2.bed > ~/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_anno.bed")
    
    #  system("intersectBed -a /home/aiminyan/chromatin_states_chromHMM_mm9/spleen_cStates_HMM.bed -b ~/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_2.bed > ~/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_anno.bed")
    
    system("intersectBed -a /home/aiminyan/chromatin_states_chromHMM_mm9/spleen_cStates_HMM.bed -b /home/aiminyan/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_2_in_mm9.bed > ~/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_anno_in_mm9.bed")
    
    peakAll.anno <- read.table("~/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_anno.bed")
    
    peakAll.anno <- read.table("~/H_driver/Aimin_project/ATAC-Seq/CountUseAllPeaks/peakAll_anno_in_mm9.bed")
    
    head(peakAll.anno)
    
    ID <- paste(as.character(peakAll.anno[,1]),peakAll.anno[,2],peakAll.anno[,3],sep = "_")
    
    peakAll.anno.2 <- cbind.data.frame(peakAll.anno,ID=ID)
    
    ID2 <- paste(as.character(peak.all.in.mm9[,1]),peak.all.in.mm9[,2],peak.all.in.mm9[,3],sep = "_")
    
    
    peakAll.anno.3 <- cbind.data.frame(peak.all.in.mm9,ID=ID2)
    
    peakAll.anno.4 <- merge.data.frame(peakAll.anno.2,peakAll.anno.3,by="ID")
    
    getState <- function(g){
        
        y <- lapply(g, function(u){
            
            if(!is.na(str_locate(u,"_")[1])){
                x <- str_sub(u,1,str_locate(u,"_")[1]-1)
            }else
            {
                x <- u
            }
            x
        })
        
        gg <- unlist(y)
        
        ggg <- gg[!is.na(gg)]
        ggg
        
    }
    
    getAnno <- function(g){
        
        y <- lapply(g, function(u){
            
            if(!is.na(str_locate(u,"_")[1])){
                x <- str_sub(u,str_locate(u,"_")[1]+1,str_length(u))
            }else
            {
                x <- u
            }
            x
        })
        
        gg <- unlist(y)
        
        ggg <- gg[!is.na(gg)]
        ggg
        
    }
    
    state <- getState(peakAll.anno.4$V4)
    anno <- getAnno(peakAll.anno.4$V4)
    
    anno.6 <- anno
    
    anno.6[which(anno.6 %in% c("Active_Promoter","Poised_Promoter"))] <- "Promoter"
    
    anno.6[which(anno.6 %in% c("Strong_Enhancer","Poised_Enhancer"))] <- "Enhancer"
    
    anno.6[which(anno.6 %in% c("Txn_Transition","Txn_Elongation","Weak_Txn"))] <- "Transcribed"
    
    peakAll.anno.5 <- cbind.data.frame(peakAll.anno.4,state=state,anno=anno,anno.6 = anno.6)
    
    mytable <- table(peakAll.anno.5$anno.6)
    
    lbls <- paste(names(mytable),mytable, sep="-")
    pie(mytable, labels = lbls,main="Annotation")
    
    mytable2 <- as.data.frame(mytable)
    colnames(mytable2)= c("Function_Annotation","counts")
    
    table_lables <- mytable2 %>%
        mutate(Function_Annotation=factor(Function_Annotation,levels=Function_Annotation[length(Function_Annotation):1]),
               cumulative=cumsum(counts),
               midpoint= cumulative-(counts/2),
               labels=paste0(round((counts/sum(counts))*100,2),"%"," (",counts,") "))
    
    ggplot(table_lables,aes(x="",y=counts,fill=Function_Annotation))+
        geom_bar(width = 1,stat="identity") +
        coord_polar(theta="y",start = 0,direction = 1) +
        scale_fill_manual(values = c("Lightblue","#AD7366","Lightgreen","Orange","Coral","Yellow"))+
        labs(x="",y="",title="Peak functional annotations\n",fill="Function_Annotation")+
        geom_text(aes(x=1.2,y=midpoint,label=labels),color="black",fontface="bold",size=3.3) +
        theme(plot.title=element_text(hjust=0.5),
              legend.title = element_text(hjust = 0.5,face="bold",size=10))
    
    ggplot(table_lables,aes(x="",y=counts,fill=Function_Annotation))+
        geom_bar(width = 1,stat="identity") +
        scale_fill_manual(values = c("Lightblue","#AD7366","Lightgreen","Orange","Coral","Yellow"))+
        labs(x="",y="",title="Peak functional annotations\n",fill="Function_Annotation")+
        geom_text(aes(x=c(1.2,1.2,1,1.2,1.3,1.2),y=midpoint,label=labels),color="black",fontface="bold",size=3.3) +
        theme(plot.title=element_text(hjust=0.5),
              legend.title = element_text(hjust = 0.5,face="bold",size=10))+ theme(axis.ticks = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background= element_blank())
    
    # +  theme_bw()
    # + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    #                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    
}


# input.bam.dir <- "~/Dropbox/Alejandro_AtacSeq_uploaded_2"
# bam.file.pattern <- "*nodup.bam$"
# output.dir <- "~/H_driver/Aimin_project/ATAC-Seq/ShiftedBam"
# generateShiftedBam(input.bam.dir,bam.file.pattern,output.dir)

generateShiftedBam <- function(input.bam.dir,bam.file.pattern,output.dir){
    
    
    library(ATACseqQC)
    library(rtracklayer)
    file.1 <- list.files(input.bam.dir,pattern=bam.file.pattern,all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    
    
    if (!dir.exists(output.dir)){dir.create(output.dir, recursive = TRUE)}
    
    
    temp.Data <- data.frame(ID=paste(basename(dirname(dirname(dirname(file.1)))),basename(dirname(file.1)),sep="_"),condition=basename(dirname(dirname(dirname(file.1)))),rep=basename(dirname(file.1)),fileName=file.1)
    
    
    temp.Data.1 <- temp.Data[-which(temp.Data$condition %in% c("IL-2vsPBS")),]
    
    
    null <- apply(temp.Data.1, 1, function(u,output.dir){
        cat(u[1],u[4],"\n")
        
        bamfile <- u[4]
        tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
        
        gal <- readBamFile(bamfile, tag=tags,asMates=TRUE)
        gal1 <- shiftGAlignmentsList(gal)
        outPath <- output.dir
        shiftedBamfile <- file.path(outPath, paste0(u[1],"_shifted.bam"))
        export(gal1, shiftedBamfile)
        
    },output.dir)
    
    
    #  bamfile <- "/home/aiminyan/Dropbox/Alejandro_AtacSeq_uploaded_2/FP-at-16/align/rep1/2017-10-26-03_S6_R1.PE2SE.nodup.bam"
    #  bamfile.labels <- "FP-at-16-rep1"
    #fragSize <- fragSizeDist(bamfile, bamfile.labels)
    
    #  gal <- readBamFile(bamfile, tag=tags,asMates=TRUE)
    #  gal1 <- shiftGAlignmentsList(gal)
    #  outPath <- "~/Dropbox/Atac-Seq/mm10"
    #  shiftedBamfile <- file.path(outPath, "FP-at-16-rep1_shifted.bam")
    #  export(gal1, shiftedBamfile)
    
}

# input.bam.dir <- "~/Dropbox/Alejandro_AtacSeq_uploaded_2"
# bam.file.pattern <- "*nodup.bam$"
# output.dir <- "~/Dropbox/Atac-Seq/ShiftedBam"

# generateShiftedBam2(input.bam.dir,bam.file.pattern,output.dir)

generateShiftedBam2 <- function(input.bam.dir,bam.file.pattern,output.dir){
    
    
    library(ATACseqQC)
    library(rtracklayer)
    file.1 <- list.files(input.bam.dir,pattern=bam.file.pattern,all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    if (!dir.exists(output.dir)){dir.create(output.dir, recursive = TRUE)}
    
    
    temp.Data <- data.frame(ID=paste(basename(dirname(dirname(dirname(file.1)))),basename(dirname(file.1)),sep="_"),condition=basename(dirname(dirname(dirname(file.1)))),rep=basename(dirname(file.1)),fileName=file.1)
    
    
    temp.Data.1 <- temp.Data[-which(temp.Data$condition %in% c("IL-2vsPBS")),]
    
    
    null <- apply(temp.Data.1, 1, function(u,output.dir){
        #cat(u[1],u[4],"\n")
        
        outPath <- output.dir
        shiftedBamfile <- file.path(outPath, paste0(u[1],"_shifted"))
        cmd = paste("perl ~/ATAC_seq_mm10/auyar/ATAC_BAM_shifter_gappedAlign.pl",u[4],shiftedBamfile,sep=" ")
        
        cat(cmd,"\n\n")
        
        system(cmd)
        
        
    },output.dir)
    
    
}

# source("https://bioconductor.org/biocLite.R")
# biocLite("ChIPseeker")
# biocLite("MotifDb")
# biocLite("MotIV")
# biocLite("ATACseqQC")
# biocLite("BSgenome.Mmusculus.UCSC.mm10")


library(MotifDb)
library(ATACseqQC)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomeInfoDb)
library(rtracklayer)
library(ChIPpeakAnno)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)


input.bam.dir <- "~/Dropbox/Atac-Seq/ShiftedBam"
output.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif"

sortBam <- function(input.bam.dir) {
    shiftedBamfile <- list.files(input.bam.dir,pattern="*shifted.bam$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    null <- lapply(shiftedBamfile, function(u){
        system(paste("sambamba sort",u,sep= " "))
    })
}


library(MotifDb)
library(ATACseqQC)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomeInfoDb)
library(rtracklayer)
library(ChIPpeakAnno)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(stringr)

# input.bam.dir <- "~/Dropbox/Atac-Seq/ShiftedBam"
# output.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif"
# findMotif(output.dir, input.bam.dir, MotifDb, Mmusculus)
#
# findMotif <- function(output.dir, input.bam.dir, MotifDb, Mmusculus) {
#
#   if (!dir.exists(output.dir)){dir.create(output.dir, recursive = TRUE)}
#
#   shiftedBamfile.sorted <- list.files(input.bam.dir,pattern="*.sorted.bam$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
#
#   CTCF <- query(MotifDb, c("CTCF"))
#   CTCF <- as.list(CTCF)
#   #print(CTCF[[1]], digits=2)
#
#   genome <- Mmusculus
#   # For all chr
#
#     null <- lapply(shiftedBamfile.sorted, function(u,output.dir,MotifDb, Mmusculus){
#
#       x <- str_sub(basename(u),1,str_locate(basename(u),"_shifted")[1]-1)
#
#       f.Fig <- file.path(output.dir,paste0(x,"_CTCF.pdf"))
#       f.text <- file.path(output.dir,paste0(x,"_CTCF.txt"))
#
#       cat(f.Fig,"\n")
#       cat(f.text,"\n")
#
#       pdf(f.Fig)
#       sigs <- factorFootprints(u, pfm=CTCF[[1]],
#                                 genome=genome,
#                                 min.score="90%",seqlev = paste0("chr", c(1:19, "X","Y")),
#                                 upstream=100, downstream=100)
#       dev.off()
#       write.table(data.frame(sigs$bindingSites),f.text,sep="\t",row.names = F,quote = F,col.names = T)
#
#     },output.dir,MotifDb, Mmusculus)
#
#   # pdf(file.path(output.dir,"FP-at-0.75_rep1_CTCF.pdf"))
#   # sigs <- factorFootprints(shiftedBamfile.sorted[1], pfm=CTCF[[1]],
#   #                          genome=genome,
#   #                          min.score="90%",seqlev = paste0("chr", c(1:19, "X","Y")),
#   #                          upstream=100, downstream=100)
#   # dev.off()
#   # write.table(data.frame(sigs$bindingSites),"~/Dropbox/Atac-Seq/FP-at-16-rep1_CTCT_all_chrs.txt",sep="\t",row.names = F,quote = F,col.names = T)
#
#
#   #STAT5 <- query(MotifDb, c("STAT5"))
#   #STAT5 <- as.list(STAT5)
#   #print(STAT5[[1]], digits=2)
# }


# input.bam.dir <- "~/Dropbox/Atac-Seq/ShiftedBam"
# output.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif"

# output.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_2"

# output.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_3"
# which.chr <- paste0("chr", c(1:19, "X","Y")

# findMotif2(output.dir, input.bam.dir, MotifDb, Mmusculus,"CTCF",3)


# output.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_stat5"
# which.chr <- paste0("chr", c(4))

# output.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_stat5_chr11"

# which.chr <- paste0("chr", c(9))
# findMotif2(output.dir, input.bam.dir, MotifDb, Mmusculus,"STAT5",2,which.chr)

# which.chr <- paste0("chr", c(8))
# findMotif2(output.dir, input.bam.dir, MotifDb, Mmusculus,"STAT5",2,which.chr)

# which.chr <- paste0("chr", c(7))
# findMotif2(output.dir, input.bam.dir, MotifDb, Mmusculus,"STAT5",2,which.chr)

# which.chr <- paste0("chr", c(1))
# findMotif2(output.dir, input.bam.dir, MotifDb, Mmusculus,"STAT5",2,which.chr)

# which.chr <- paste0("chr", c(10))
# findMotif2(output.dir, input.bam.dir, MotifDb, Mmusculus,"STAT5",2,which.chr)

# which.chr <- paste0("chr", c(11))
# findMotif2(output.dir, input.bam.dir, MotifDb, Mmusculus,"STAT5",2,15,which.chr)

# which.chr <- paste0("chr", c(11))
# findMotif2(output.dir, input.bam.dir, MotifDb, Mmusculus,"STAT5",2,15,which.chr,anchor="fragment center")

# CTCF
# output.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_ctcf"
# which.chr <- paste0("chr", c(1))
# findMotif2(output.dir, input.bam.dir, MotifDb, Mmusculus,"CTCF",3,100,which.chr,anchor="fragment center")

# CTCF
# output.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_ctcf_cut_site"
# which.chr <- paste0("chr", c(1))
# findMotif2(output.dir, input.bam.dir, MotifDb, Mmusculus,"CTCF",3,100,which.chr,anchor="cut site")


findMotif2 <- function(output.dir, input.bam.dir, MotifDb, Mmusculus,motif,which.motif,dd,which.chr,...) {
    
    if (!dir.exists(output.dir)){dir.create(output.dir, recursive = TRUE)}
    
    shiftedBamfile.sorted <- list.files(input.bam.dir,pattern="*.sorted.bam$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    CTCF <- query(MotifDb, c(motif))
    CTCF <- as.list(CTCF)
    #print(CTCF[[1]], digits=2)
    
    genome <- Mmusculus
    # For all chr
    
    null <- lapply(shiftedBamfile.sorted, function(u,output.dir,CTCF,genome,which.motif,which.chr){
        
        x <- str_sub(basename(u),1,str_locate(basename(u),"_shifted")[1]-1)
        
        f.Fig <- file.path(output.dir,paste0(x,"_",motif,"_",which.chr,".pdf"))
        f.text <- file.path(output.dir,paste0(x,"_",motif,"_",which.chr,".txt"))
        
        cat(u,"\t",f.Fig,"\n")
        cat(u,"\t",f.text,"\n")
        
        pdf(f.Fig)
        sigs <- factorFootprints(u, pfm=CTCF[[which.motif]],
                                 genome=genome,
                                 min.score="90%",seqlev = which.chr,
                                 upstream=dd, downstream=dd,legTitle = x,...)
        dev.off()
        write.table(data.frame(sigs$bindingSites),f.text,sep="\t",row.names = F,quote = F,col.names = T)
        
    },output.dir,CTCF,genome,which.motif,which.chr)
    
    
    #system("pdftk /home/aiminyan/Dropbox/Atac-Seq/motif/*.pdf cat output 18_2.pdf")
    
    
    # pdf(file.path(output.dir,"FP-at-0.75_rep1_CTCF.pdf"))
    # sigs <- factorFootprints(shiftedBamfile.sorted[1], pfm=CTCF[[1]],
    #                          genome=genome,
    #                          min.score="90%",seqlev = paste0("chr", c(1:19, "X","Y")),
    #                          upstream=100, downstream=100)
    # dev.off()
    # write.table(data.frame(sigs$bindingSites),"~/Dropbox/Atac-Seq/FP-at-16-rep1_CTCT_all_chrs.txt",sep="\t",row.names = F,quote = F,col.names = T)
    
    
    #STAT5 <- query(MotifDb, c("STAT5"))
    #STAT5 <- as.list(STAT5)
    #print(STAT5[[1]], digits=2)
}

# Ex1:
# input.file.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_3"
# input.file.pattern <- "*.pdf"
# output.file.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_3"
# output.file.name <- "CTCF_4_18_samples.pdf"
# combinePdf(input.file.dir,input.file.pattern,output.file.dir,output.file.name)

# Ex2:
# input.file.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_stat5"
# input.file.pattern <- "*STAT5_chr10.pdf"
# output.file.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_stat5"
# output.file.name <- "stat5_4_18_samples_chr10.pdf"
# combinePdf(input.file.dir,input.file.pattern,output.file.dir,output.file.name)

# For chr11:
# input.file.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_ctcf_cut_site"
# input.file.pattern <- "*STAT5_chr11.pdf"
# output.file.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_ctcf_cut_site"
# output.file.name <- "stat5_4_18_samples_chr11.pdf"
# combinePdf(input.file.dir,input.file.pattern,output.file.dir,output.file.name)

# For CTCF in chr1:
# input.file.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_ctcf_cut_site"
# input.file.pattern <- "*CTCF_chr1.pdf"
# output.file.dir <- "/home/aiminyan/Dropbox/Atac-Seq/motif_ctcf_cut_site"
# output.file.name <- "ctcf_4_18_samples_chr1.pdf"
# combinePdf(input.file.dir,input.file.pattern,output.file.dir,output.file.name)

combinePdf <- function(input.file.dir,input.file.pattern,output.file.dir,output.file.name){
    
    if (!dir.exists(output.file.dir)){dir.create(output.file.dir, recursive = TRUE)}
    
    cmd=paste("pdftk",file.path(input.file.dir,input.file.pattern),"cat","output",file.path(output.file.dir,output.file.name))
    
    cat(cmd,"\n")
    
    system(cmd)
    
}

# input.bam.dir <- "~/Dropbox/Atac-Seq/ShiftedBam"
# option <- "--binSize 10 --normalizeTo1x 2652783500"
# output.file.dir <- "~/Dropbox/Atac-Seq/Bw4ShiftedBam_bin_10"
# convertBam2Bw(input.bam.dir,option,output.file.dir)

convertBam2Bw <- function(input.bam.dir,option,output.file.dir){
    
    shiftedBamfile.sorted <- list.files(input.bam.dir,pattern="*.sorted.bam$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
    
    if (!dir.exists(output.file.dir)){dir.create(output.file.dir, recursive = TRUE)}
    
    null <- lapply(shiftedBamfile.sorted, function(u,output.file.dir){
        
        xx <- u
        x <- str_sub(basename(u),1,str_locate(basename(u),"_shifted")[1]-1)
        
        
        f.bw <- file.path(output.file.dir,paste0(x,".bw"))
        
        cmd1 <- paste("bamCoverage --bam",xx,sep=" ")
        cmd2 <- "-o"
        
        
        cmd3 <- paste(cmd1,cmd2,f.bw,option,sep=" ")
        
        cat(cmd3,"\n")
        system(cmd3)
        
    },output.file.dir)
    
}

# system("find ~/Dropbox/Atac-Seq/Bw4ShiftedBam_bin_10/ -type f -name '*.bw' | xargs /home/aiminyan/Dropbox/Aimin_project/AtacSeq/inst/UCSC-custom-track-generator/make-tracks.py -url https://miamiedu-my.sharepoint.com/:f:/g/personal/axy148_miami_edu/EismF-w1t91AoJ4fDoulY8YBeNGECdp6pkvTQi2LIe6yng/ -o ~/Dropbox/Atac-Seq/Bw4ShiftedBam_bin_10/my_custom_tracks.txt")


# set a account at https://user.cyverse.org/register
# login in at https://de.cyverse.org/de
#
#
# Install iCommands by following
# https://pods.iplantcollaborative.org/wiki/display/DS/Setting+Up+iCommands#SettingUpiCommands-co
#
# Use icommands by following https://docs.irods.org/master/icommands/user/
#
# iinit
# ils -A
# imkdir AtacSeq
# iput -bf ~/Dropbox/Atac-Seq/Bw4ShiftedBam_bin_10/*.bw AtacSeq
#
# After upload to your space
# go to your space, select "Name", share-> View in Genome Browser to get the following urls for makeUcscTrack function
#

url=c("https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/FP-at-0.75_rep1.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/FP-at-0.75_rep2.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/FP-at-0.75_rep3.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/FP-at-16_rep1.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/FP-at-16_rep2.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/FP-at-16_rep3.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/IL-2-at-16_rep1.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/IL-2-at-16_rep2.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/IL-2-at-16_rep3.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/IL2-at-0.75_rep1.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/IL2-at-0.75_rep2.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/IL2-at-0.75_rep3.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/PBS-at-0.75_rep1.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/PBS-at-0.75_rep2.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/PBS-at-0.75_rep3.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/PBS-at-16_rep1.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/PBS-at-16_rep2.bw
https://de.cyverse.org/anon-files/iplant/home/axy148/AtacSeq/PBS-at-16_rep3.bw")

output.dir <- "~/Dropbox/Atac-Seq/Bw4ShiftedBam_bin_10"
output.file <- "18_AtacSeq_custom_tracks.txt"

makeUcscTrack(url,output.dir,output.file)

makeUcscTrack <- function(url,output.dir,output.file,smoothingWindow="4",color="123,100,50",autoScale="on",viewLimits="1:200",visibility="full",windowingFunction="maximum"){
    
    if (!dir.exists(output.dir)){dir.create(output.dir, recursive = TRUE)}
    
    ucsc.track <- cbind.data.frame(name=paste0("name=",tools::file_path_sans_ext(basename(unlist(strsplit(url,"\n"))))),bigDataUrl=paste0("bigDataUrl=",unlist(strsplit(url,"\n"))))
    
    n <- dim(ucsc.track)[1]
    
    track <- rep("track",n)
    type <- rep("type=bigWig",n)
    sW <- rep(paste0("smoothingWindow=",smoothingWindow),n)
    color <- rep(paste0("color=",color),n)
    autoScale <- rep(paste0("autoScale=",autoScale),n)
    viewLimits <- rep(paste0("viewLimits=",viewLimits),n)
    visibility <- rep(paste0("visibility=",visibility),n)
    windowingFunction <- rep(paste0("windowingFunction=",windowingFunction),n)
    
    ucsc.track.1 <- cbind.data.frame(track,type,name=ucsc.track[,1],sW,color,autoScale,viewLimits,visibility,windowingFunction,bigDataUrl=ucsc.track[,2])
    
    write.table(ucsc.track.1,file = file.path(output.dir,output.file),quote = F, sep = " ",row.names = FALSE,
                col.names = FALSE)
    
}
