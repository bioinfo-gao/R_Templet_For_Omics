# Mutation Annotation Format (MAF) 
# This package attempts to summarize, analyze, annotate and visualize MAF files any studies in MAF format, including TCGA 
# maftools functions can be categorized into mainly Visualization and Analysis modules.
# Install from Bioconductor repository 
BiocManager::install("maftools")

# install.packages("mclust")

#Install from GitHub repository
# BiocManager::install("PoisonAlien/maftools")

library(maftools)

# path to TCGA LAML MAF file
# AML as the abbreviation for acute myeloid leukemia
# LAML used within TCGA for unkown reason 
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations

laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 

#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

laml = read.maf(maf = laml.maf, clinicalData = laml.clin)

laml

#For VCF files or simple tabular files, easy option is to use vcf2maf utility 
# ANNOVAR for variant annotations, maftools has a handy function annovarToMaf