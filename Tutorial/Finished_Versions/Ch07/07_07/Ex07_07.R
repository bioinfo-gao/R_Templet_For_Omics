# R Statistics Essential Training
# Ex07_07
# Creating crosstabs for categorical variables

# Load data
?Titanic
str(Titanic)
Titanic
ftable(Titanic)  # Makes "flat" table

# Convert table to 1 table
as.data.frame.table(Titanic)
as.data.frame.table(Titanic)$Freq
?rep
rep(1:4, 2)
lapply(as.data.frame.table(Titanic), function(x)rep(x, as.data.frame.table(Titanic)$Freq))

as.data.frame(lapply(as.data.frame.table(Titanic), function(x)rep(x, as.data.frame.table(Titanic)$Freq)))
as.data.frame(lapply(as.data.frame.table(Titanic), function(x)rep(x, as.data.frame.table(Titanic)$Freq)))[, -5]

# Convert table to data frame with one row per observation
tdf <- as.data.frame(lapply(as.data.frame.table(Titanic), function(x)rep(x, as.data.frame.table(Titanic)$Freq)))[, -5]
tdf[1:5, ]  # Check first five rows of d
tdf[]  # Check first five rows of 



# Create contingency table
ttab <- table(tdf$Class, tdf$Survived)
ttab

ttab1 <- table(tdf$Sex, tdf$Survived)
ttab1

ttab2 <- table(tdf$Age, tdf$Survived)
ttab2



# Call also get cell, row, and column %
# With rounding to get just 2 decimal places
# Multiplied by 100 to make %
ttab
round(prop.table(ttab, 1), 2) * 100 # row %, here 1 means row, 2 means 2 digit for percent
round(prop.table(ttab, 2), 2) * 100 # column %
round(prop.table(ttab), 2) * 100    # cell %

# Chi-squared test
tchi <- chisq.test(ttab)
tchi

# Additional tables
tchi$observed   # Observed frequencies (same as ttab)
tchi$expected   # Expected frequencies
tchi$residuals  # Pearson's residual
tchi$stdres     # Standardized residual

rm(list = ls())  # Clean up