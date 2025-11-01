
df <- data.frame(
    fruit = c("Banana", "Banana", "Orange", "Orange", "Banana", "Orange","Banana", "Banana", "Orange", "Orange", "Banana", "Orange"),
    color = c("Yellow", "Yellow", "Orange", "Orange", "Yellow", "Yellow","Yellow", "Yellow", "Orange", "Orange", "Yellow", "Orange")
)

df 

# Frequency of multiple columns (contingency table)
freq_table_multi <- table(df$fruit, df$color)
print(freq_table_multi)

chisq.test(freq_table_multi)     

M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
dimnames(M) <- list(gender = c("F", "M"),
                    party = c("Democrat","Independent", "Republican"))
M

(Xsq <- chisq.test(M))  # Prints test summary
Xsq$observed   # observed counts (same as M)
Xsq$expected   # expected counts under the null
Xsq$residuals  # Pearson residuals
Xsq$stdres     # standardized residuals
