# 安装和加载 brms
#install.packages("brms")
library(brms)

# 数据：美国和全球 iPhone 销量（百万部，转换为整数）
data <- data.frame(
    year = c(2022, 2023, 2024),
    us_sales = round(c(72.28, 75, 72)),  # round() 将销量数据四舍五入为整数 
    global_sales = round(c(231.8, 231.8, 232.5))  # 四舍五入为整数：232, 232, 233
)

# 计算比例（仅用于参考）
data$prop <- data$us_sales / data$global_sales

# 贝叶斯二项式模型：估计美国销量占比
fit <- brm(
    us_sales | trials(global_sales) ~ 1,  # 二项式模型：us_sales 成功次数，global_sales 总试验次数
    data = data,
    family = binomial(link = "logit"), # need integer
    prior = prior(normal(0, 1), class = "Intercept"),  # 弱信息先验
    chains = 4, iter = 2000, warmup = 1000, seed = 42
)

# 输出结果
summary(fit)

# 可视化后验分布
library(bayesplot)
posterior <- as_draws_df(fit)
mcmc_areas(posterior, pars = "b_Intercept", prob = 0.95) +
    labs(title = "后验分布：美国 iPhone 销量占比 (2022-2024)")

