#source("http://www.bioconductor.org/biocLite.R")
#biocLite(c("GenomicRanges"))
library("GenomicRanges")

dir<-"/media/2T-Disk/DMS/comb-p-and-data5"
setwd(dir)

sink("output.txt")
#suppressWarnings()
for (i in 1:50) { 
    #i <- 3
    # 1.1 read in chromsome and locus info
    train_file_name <- paste0(dir, "/Train_", toString(i), "_anno.hg19.bed")
    if (file.size(train_file_name) == 0) next
    
    train_file <- read.table(train_file_name, header=T, row.names=NULL)
    colnames(train_file) <- c(colnames(train_file)[-1],"x")
    train_file$x <- NULL
    
    train_file[,2] <- as.character(train_file[,2])
    train_file[,3] <- as.character(train_file[,3])
    
    df1 <- makeGRangesFromDataFrame(train_file, keep.extra.columns = TRUE, start.field = "start", end.field = "end")
    #df1 <- df1[c(1,2,3)]
     
    #test
    test_file_name <- paste0(dir, "/Test_", toString(i), "_anno.hg19.bed")
    if (file.size(test_file_name) == 0) next
   
    
    test_file <- read.table(test_file_name, header=T, row.names=NULL)
    colnames(test_file) <- c(colnames(test_file)[-1],"x")
    test_file$x <- NULL
    
    test_file[,2] <- as.character(test_file[,2])
    test_file[,3] <- as.character(test_file[,3])
    
    df2 <- makeGRangesFromDataFrame(test_file, keep.extra.columns = TRUE, start.field = "start", end.field = "end")
    
    #df2 <- df2[c(1,2,3)]
    #suppressWarnings()
    hits <- suppressWarnings(findOverlaps(df1, df2))
    hits.df = data.frame(df1[queryHits(hits),c(1,2,3)], df2[subjectHits(hits),c(1,2,3)]) # important
    
    if (length(hits)!=0) {
        statement <- paste0("The cycle No. ", toString(i), " has matched locus in training and test splits")
        print(statement)
        print(hits.df)
        #print()
    }
}

sink()

