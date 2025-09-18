# R Statistics Essential Training
# Ex04_01
# Examining outliers

# Categorical data
# Outlier is < 10%
# Worldwide shipments of smartphone OS
# in millions for 2013 Q1
OS <- read.csv("~/Desktop/R/OS.csv", header = TRUE)
View(OS)
OS

# Outlier has proportion < .10
# Either combine into "other" (if homogeneous) or delete
OS.hi <- subset(OS, Proportion > 0.1)
OS.hi

# Quantitative data
# See outliers in boxplots
require("datasets")
?rivers
rivers
data(rivers)  # Lengths of Major North American Rivers
hist(rivers)
boxplot(rivers, horizontal = TRUE)
boxplot.stats(rivers)
rivers.low  <- rivers[rivers < 1210]  # Remove outliers
boxplot(rivers.low, horizontal = TRUE)  # Has new outliers
boxplot.stats(rivers.low)
rivers.low2  <- rivers[rivers < 1055]  # Remove again
boxplot(rivers.low2)  # Still one outlier

rm(list = ls())  # Clean up

c1=c(1,2,3,4,5)
c2=c(1,2,3,4,5)
cp=c1*c2
cp
