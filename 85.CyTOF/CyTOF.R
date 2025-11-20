# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("diffcyt")
# 
# BiocManager::install("HDCytoData")
# BiocManager::install("CATALYST")

# https://bioconductor.org/packages/devel/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html

library(cytolib)
library("HDCytoData")

d_flowSet <- Bodenmiller_BCR_XL_flowSet()
suppressPackageStartupMessages(library(flowCore))


d_flowSet
pData(d_flowSet)
fsApply(d_flowSet, nrow)
dim(exprs(d_flowSet[[1]]))
exprs(d_flowSet[[1]])[1:6, 1:5]

# Alternatively: load data from '.fcs' files
# files <- list.files(
#     path = "path/to/files", pattern = "\\.fcs$", full.names = TRUE
# )
# d_flowSet <- read.flowSet(
#     files, transformation = FALSE, truncate_max_range = FALSE
# )

# Meta-data: experiment information

# check sample order
filenames <- as.character(pData(d_flowSet)$name)

# sample information
sample_id <- gsub("^PBMC8_30min_", "", gsub("\\.fcs$", "", filenames))
group_id <- factor(
    gsub("^patient[0-9]+_", "", sample_id), levels = c("Reference", "BCR-XL")
)
patient_id <- factor(gsub("_.*$", "", sample_id))

experiment_info <- data.frame(
    group_id, patient_id, sample_id, stringsAsFactors = FALSE
)
experiment_info


# Meta-data: marker information

# source: Bruggner et al. (2014), Table 1

# column indices of all markers, lineage markers, and functional markers
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)

# channel and marker names
channel_name <- colnames(d_flowSet)
marker_name <- gsub("\\(.*$", "", channel_name)

# marker classes
# note: using lineage markers for 'cell type', and functional markers for 
# 'cell state'
marker_class <- rep("none", ncol(d_flowSet[[1]]))
marker_class[cols_lineage] <- "type"
marker_class[cols_func] <- "state"
marker_class <- factor(marker_class, levels = c("type", "state", "none"))

marker_info <- data.frame(
    channel_name, marker_name, marker_class, stringsAsFactors = FALSE
)
marker_info
marker_info$channel_name



