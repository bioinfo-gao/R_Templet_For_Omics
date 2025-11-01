install.packages("rsthemes")
install.packages("cli")
# install.packages(
#   "rsthemes",
#   repos = c(gadenbuie = 'https://gadenbuie.r-universe.dev', getOption("repos"))
# )

# Or you can install rsthemes from GitHub with:
  
install.packages("devtools")
# or
devtools::install_github("gadenbuie/rsthemes") # pretty theme 

## larger collecton, but I don't need 
# devtools::install_github("max-alletsee/rstudio-themes") 
# Then, install the included, hand-crafted themes with:

library(rsthemes)
  
rsthemes::install_rsthemes()

rsthemes::list_rsthemes() #  to list installed themes
rsthemes::try_rsthemes() #   to try all installed themesnn

rstudioapi::applyTheme("Oceanic Plus {rsthemes}")     # 1
rstudioapi::applyTheme("One Dark {rsthemes}")         # 2
rstudioapi::applyTheme("Serendipity Dark {rsthemes}") # 3
rstudioapi::applyTheme("a11y-dark {rsthemes}")        # 4       
rstudioapi::applyTheme("Fairyfloss {rsthemes}")       # 5

rstudioapi::applyTheme("Material Darker {rsthemes}")  # 
rstudioapi::applyTheme("Serendipity Light {rsthemes}")
rstudioapi::applyTheme("Yule RStudio {rsthemes}")
rstudioapi::applyTheme("Material Ocean {rsthemes}")
