#bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# WARNING: Rtools is required to build R packages but is not currently installed.
# Please download and install the appropriate version of Rtools before proceeding:
BiocManager::install("ggtree")
#BiocManager::install("ggtreeExtra")
BiocManager::install("SummarizedExperiment")
BiocManager::install("GEOquery")
# Alternatively, install just pillar:
install.packages('lifecycle')
install.packages('vctrs')
install.packages("pillar")

# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("r-lib/pillar")



getGEO()

install.packages("rtools")
install.packages("devtools")
install.packages("ggplot2")
install.packages("devtools")
install.packages("roxygen2")
install.packages("RCurl")
install.packages("randomForest")

devtools::install_github("r-lib/devtools")

# https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html
#require("rtools")
require("devtools")
require("ggplot2")
require("devtools")
require("roxygen2")

devtools::create("myfirstpackage")
# ✔ Creating 'myfirstpackage/'
# ✔ Setting active project to 'C:/Users/zhen-/Documents/myfirstpackage'
# ✔ Creating 'R/'
# ✔ Writing 'DESCRIPTION'
# Package: myfirstpackage
# Title: What the Package Does (One Line, Title Case)
# Version: 0.0.0.9000
# Authors@R (parsed):
#     * First Last <first.last@example.com> [aut, cre] (YOUR-ORCID-ID)
# Description: What the package does (one paragraph).
# License: `use_mit_license()`, `use_gpl3_license()` or friends to
# pick a license
# Encoding: UTF-8
# Roxygen: list(markdown = TRUE)
# RoxygenNote: 7.3.2
# ✔ Writing 'NAMESPACE'
# ✔ Writing 'myfirstpackage.Rproj'
# ✔ Adding '^myfirstpackage\\.Rproj$' to '.Rbuildignore'
# ✔ Adding '.Rproj.user' to '.gitignore'
# ✔ Adding '^\\.Rproj\\.user$' to '.Rbuildignore'
# ✔ Setting active project to '<no active project>'
# # https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html
#
# This ends up creating a folder with the same name as your package name with 4 files inside the folder:
#
# DESCRIPTION: This is where all the meta-data about your package goes. Rather than try to explain the contents,
# I will refer you to Hadley’s detailed explanation on the contents of this file.
# === myfirstpackage.Rproj: This is a RStudio specific file. As I do not use RStudio,
# I will not comment on this file as I never use it.
# === NAMESPACE: In short, this file indicates what needs to be exposed to users for your R package. From my experience,
# I’ve never edited this file as devtools takes care of the changes as you’ll see below.
# === R: This is where all your R code goes for your package.
# You now have the bare bones of your first R package. First start by filling out the details in the DESCRIPTION file.
# When that is done, we can start adding some functions!


# https://stackoverflow.com/questions/73990995/r-cmd-check-warning-non-standard-license-specification
# use_mit_license()
#
# #　https://stackoverflow.com/questions/54419777/devtoolsuse-data-internal-true
#
# devtools::use_data has now been deprecated and transferred to the usethis package, hence the message; the functions remains the same.
#
# If you set internal = TRUE the data will only be accessible by functions of the package internally while if FALSE the data will be exported, users will be able to load it with data("aaa").
#
# There is no need for multiple files, you can save multiple variables in a single file, i.e.: usethis::use_data(cars, mtcars).

# usethat::use_vignette("introduction")
devtools::install_github("yourusername/myfirstpackage")

source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite("Rsamtools")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.19")
#To install core packages, type the following in an R command window:
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
library("GenomicFeatures", "AnnotationDbi")
BiocManager::install("simpleSingleCell")
install.packages("foreign", "MASS", "nlme", "survival")

require("foreign")
require("MASS")
require("nlme")
require("survival")

# It was a permission issue for me. First, I identified where
# the packages were installed using installed.packages()[, c("Package", "LibPath")]. This outputs a long 2 column matrix with the names and locations of the packages. Then you will see where the offending packages are. In my case, they were at /usr/lib/R/site-library and /usr/lib/R/library. Then I changed the permission of these folders by chmod (I used chmod -R 777 on the main R folder, this is my personal computer, so security is not a big concern here I think).

help(package = "dplyr")
help(package = "data.table")
help(package = "MASS")
help(package = "devtools")
# https://github.com/r-lib/devtools/issues/2145
name = "Introduction_to_myfirstpackag"
usethis::use_vignette(name, title = name)
?usethis::use_vignette

BiocManager::install("coMethDMR")

library("SummarizedExperiment")
library("GEOquery")
