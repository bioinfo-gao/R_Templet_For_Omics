# https://www.youtube.com/watch?v=dxBR191Vc8A
#install.packages("curl", dependencies = T)
install.packages("chattr", dependencies = T)
uninstall.packages("chattr")
# At least in Windows, there is a popup show that "do you want to install the library need compilation" , I choose NO
# there is a list of compilateion needed version and sereral not show in Rstudio
#install.packages("chattr")
library(chattr)

##  Setting Up chattr in RStudio: A Complete Guide
# https://rstudio-pubs-static.s3.amazonaws.com/1360735_64ad9b3ba74d4f60bff2f36e94e2ab35.html
# https://cloud.r-project.org/web/packages/chattr/readme/README.html

## Get an Open AI PI key
# 1) create an account
# 2) setup payment method with limit
# 3) create an API key

Sys.getenv("OPENAI_API_KEY")

install.packages("usethis")
# need click the popup in windows 
usethis::edit_r_environ()
.rs.restartR()
# ☐ Modify C:/Users/zhen-/.Renviron.
# ☐ Restart R for changes to take effect.
# https://aistudio.google.com/app/api-keys

chattr_use("gpt41")
# Open chat interface
chattr()



chattr("How do I create a ggplot2 scatter plot?")

# Change model on the fly
chattr_use("gpt-4")  # Switch to GPT-4
chattr_use("claude-3")  # Switch to Claude 3
chattr_use("gemini-pro")  # Switch to Gemini Pro

# View available models
chattr_available()

# Configure defaults
chattr_defaults(
    model = "gpt-4",
    max_tokens = 1000,
    temperature = 0.7
)

?chattr_use
chattr_app()
# Check current settings
chattr_defaults()

# Test the connection
chattr_test()


# Option 1: Provide a prompt directly
#  chattr("How do I make a scatter plot with ggplot2?")
# Option 2: Use the Interactive App (Recommended)
library(chattr)
chattr_app()


chattr()
