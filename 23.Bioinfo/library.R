install.packages("Rtsne")
install.packages("RCurl") # library("RCurl") #sudo apt-get install libcurl4-gnutls-dev
install.packages("scran")
install.packages("pheatmap")
install.packages("UpSetR")
install.packages("VennDiagram")

library(UpSetR)

if (!require("BiocManager", quietly = TRUE))  install.packages("BiocManager")

# There's a trick to this where one needs to add biocViews: to the package Description. 
# That's the only solution I've ever seen to allowing automatic installation of bioconductor dependencies.
BiocManager::install("GenomicRanges")
BiocManager::install("SummarizedExperiment")
BiocManager::install("DESeq2")
BiocManager::install("SRAdb")

library(SRAdb)

help(SRAdb) # BiocManager::install("getSRAfile") #https://www.biostars.org/p/93494/


# https://www.biostars.org/p/93494/
library(remotes)
install_version("foobarbaz", "0.1.2")
#An alternative is to install from the GitHub CRAN mirror.

library(remotes)
install_github("cran/foobarbaz")

library("GenomicRanges","SummarizedExperiment","DESeq2")
library("getSRAfile")

# SRP133642 -O /home/gao/Desktop/Code/Bioinfo/scRNA/
getSRAfile(in_acc = c("SRP133642"), sra_con = sra_con,
           destDir = "/home/gao/Desktop/Code/Bioinfo/scRNA/", fileType = 'sra', srcType='ftp')