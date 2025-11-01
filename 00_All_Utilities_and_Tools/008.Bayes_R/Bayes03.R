n                           <- 500
diet_effect                 <- 0.1        # the diet effect 
hidden_factor_effect        <- c(0, 0.95) # the hidden factor is a gene gene effect 
names(hidden_factor_effect) <- c("FA", "FB")

set.seed(123) 

hidden_factor_chance        <- runif(n)
hidden_factor               <- ifelse( hidden_factor_chance < 0.9, "FA", "FB")
table(hidden_factor)

random_data <- runif(n)
random_data
group <- ifelse( random_data < 0.5, "control", "drug")

# group_and_data <- cbind(group, random_data )
# group_and_data
# head(group_and_data)

diet.chance <- runif(n)
drug.chance <- runif(n)

hidden_factor
hidden_factor_effect 
hidden_factor_effect[hidden_factor] 

# diet          <- 0.1
hidden_facto_effect 
outcome <- ((drug.chance < diet ) | (drug.chance <hidden_factor_effect[hidden_factor] * (group == "drug")))
outcome
table(outcome)
summary(outcome)#summary(factor(outcome))

trail <- data.frame(group = group, Hidden_factor=hidden_factor, treatment= outcome)
trail
summary(trail)
head(trail)
str(trail)

with(trail[group == "control",], table(Hidden_factor, treatment))
with(trail[group ==    "drug",], table(Hidden_factor, treatment))


trail2 <- data.frame(group = as.factor(group),  treatment= as.factor(outcome))
trail2

freq_table_multi <- table(trail2$group, trail2$treatment)
print(freq_table_multi)

chisq.test(freq_table_multi)
