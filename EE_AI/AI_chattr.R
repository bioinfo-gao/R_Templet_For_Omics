# https://www.youtube.com/watch?v=dxBR191Vc8A
installed.packages("chattr")
library(chattr)

##  Setting Up chattr in RStudio: A Complete Guide
# https://rstudio-pubs-static.s3.amazonaws.com/1360735_64ad9b3ba74d4f60bff2f36e94e2ab35.html
# https://cloud.r-project.org/web/packages/chattr/readme/README.html

## Get an Open AI PI key
# 1) create an account
# 2) setup payment method with limit
# 3) create an API key


Sys.getenv("OPENAI_API_KEY")
chattr_app()
