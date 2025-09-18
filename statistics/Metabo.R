# https://www.metaboanalyst.ca/docs/RTutorial.xhtml
metanr_packages <- function(){
    metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea", "devtools", "crmn")
    list_installed <- installed.packages()
    new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
    if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        BiocManager::install(new_pkgs)
        print(c(new_pkgs, " packages added..."))
    }
    
    if((length(new_pkgs)<1)){
        print("No new packages added...")
    }
}
metanr_packages()
install.packages("pacman")

library(pacman)

pacman::p_load(c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea"))


#'RBGL', 'fgsea', 'impute', 'pcaMethods', 'MSnbase', 'siggenes'

BiocManager::install('RBGL')
BiocManager::install('fgsea')
BiocManager::install('impute')
BiocManager::install('pcaMethods')
BiocManager::install('MSnbase')
BiocManager::install('siggenes')
BiocManager::install('OptiLCMS')

# Step 1: Install devtools
install.packages("devtools")
install.packages("OptiLCMS")
library(devtools)

# Step 2: Install MetaboAnalystR without documentation
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)
devtools::install_github("xia-lab/OptiLCMS", build = TRUE, build_vignettes = FALSE, build_manual =TRUE)

library(googledrive);
temp <- tempfile(fileext = ".csv")
dl2 <- drive_download(
    as_id("1KaBnSNRrirVPvpRxIubGCqpjX8asNeVA"), path = temp, overwrite = TRUE)



library("MetaboAnalystR")

# Load MetaboAnalystR
library(MetaboAnalystR)
# Load OptiLCMS
library(OptiLCMS)
