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



library("SummarizedExperiment")
library("GEOquery")


getGEO()

install.packages("rtools")
install.packages("devtools")
install.packages("ggplot2")
install.packages("devtools")
install.packages("roxygen2")
install.packages("data.table")


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

# DESCRIPTION file:
# Notice how I also specified the version of the data.table. Basically I am saying that this package I am building requires the data.table package and specifically a version of it that is >= 1.9.4.
# You can indicate multiple external dependencies by just adding them in the next line:
# ===========
# Imports:
#     data.table (>= 1.9.4),
#     dplyr
# ===========
# Notice how I didn’t specify any version for dplyr which simply indicates that the package requires some version of dplyr.
# Also remember the comma between each dependency. I’ve been burned a few times by that!

# You can additionally add packages to the
# Depends section of the DESCRIPTION file instead of the Imports section.
# What’s the difference? The only difference is that packages in the Depends section are loading and attached while packages in the Imports section are only loaded. For more details on this, I refer you to the namespace section of R packages.

devtools::document() # This will generate the load_mat.Rd file in the man folder:
devtools::check()

x <- c(1:10)
devtools::use_data(x)

