# "koenrich.xls" 
# The term "KO enrichment" refers to KEGG Orthology enrichment, a common bioinformatics analysis
# 
# A koenrich.xls file is a standard output from bioinformatics tools like DAVID or those based on the clusterProfiler
#  representing the results of a KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway enrichment analysis. 
#  The spreadsheet show which KEGG pathways are statistically over-represented (or "enriched") in a given gene list. 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")  

browseVignettes("clusterProfiler")





# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
# 
library(clusterProfiler)
search_kegg_organism('ece', by='kegg_code')

ecoli <- search_kegg_organism('Escherichia coli', by='scientific_name')


?search_kegg_organism

# ecoli <- search_kegg_organism('Homo sapiens', by='scientific_name')
# dim(ecoli)
# head(ecoli)

# 7.2 KEGG pathway over-representation analysis

data(geneList, package="DOSE")

str(geneList)


gene <- names(geneList)[abs(geneList) > 2]

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

setwd("C:/Users/zhen-/Documents/Code/Rcode/R_Templet_For_Omics/plotting/Publish_quanlity_plots/")
write.csv(kk, "koenrich.csv", row.names = FALSE)

#install.packages("writexl")
library("writexl")

write.csv(kk, "koenrich.xls", row.names = FALSE)










# Perform KEGG enrichment analysis
kegg_results <- enrichKEGG(gene = my_genes, 
                           organism = 'hsa', # 'hsa' for human. See search_kegg_organism() for others.
                           keyType = 'kegg')

# View the results
head(kegg_results)



# 7.3 KEGG pathway gene set enrichment analysis
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)


# 7.4 KEGG module over-representation analysis
# KEGG Module is a collection of manually defined function units. In some situation, KEGG Modules have a more straightforward interpretation.

mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
head(mkk)                   


# 7.5 KEGG module gene set enrichment analysis
mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'hsa',
                 pvalueCutoff = 1)
head(mkk2)

# 7.6 Visualize enriched KEGG pathways
browseKEGG(kk, 'hsa04110')

# another way to view
#BiocManager::install("pathview")  
library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))
# the plot was produced in wd ==>> Working in directory C:/Users/zhen-/Documents/Code/Rcode/R_Templet_For_Omics