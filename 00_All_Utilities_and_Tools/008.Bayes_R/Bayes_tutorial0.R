# 安装必要包（如果未安装）
install.packages(c("brms", "rstan", "ggplot2"))

# 加载包
vignette("brms")
library(brms)
library(ggplot2)

# 生成模拟数据：简单线性回归示例（y = a + b*x + error）
set.seed(42)
n <- 100
x <- rnorm(n)
y <- 2 + 3 * x + rnorm(n, sd = 1)

# 创建数据框
data <- data.frame(x = x, y = y)
head(data)
dim(data)

# 拟合贝叶斯线性回归模型
# 使用默认先验（Normal(0,10)），4条链，
# ==== 2000次迭代
fit <- brm(y ~ x, data = data, family = gaussian(),
           chains = 4, iter = 2000, warmup = 1000, thin = 1,
           seed = 42)

# 查看模型摘要（参数估计、后验分布）
summary(fit)

# 绘制后验分布图（例如参数 b 的分布）
plot(fit, pars = "^b_")

# 提取后验样本
posterior <- as_draws_df(fit)

# 绘制参数 b 的密度图
ggplot(posterior, aes(x = b_x)) +
    geom_density(fill = "blue", alpha = 0.5) +
    labs(title = "Posterior Distribution of Slope (b_x)", x = "b_x", y = "Density") +
    theme_minimal()

# 预测新数据
new_data <- data.frame(x = seq(-3, 3, length.out = 50))
pred <- predict(fit, newdata = new_data, summary = TRUE)

head(new_data)
dim(new_data)
summary(new_data)

head(data)
dim(data)
summary(data)

head(pred)
dim(pred)
summary(pred) # $ operator is invalid for atomic vectors
typeof(pred) # $ operator is invalid for atomic vectors
class(pred) # $ operator is invalid for atomic vectors Here is matrix

pred=as.data.frame(pred)

# 绘制预测间obooks
ggplot() +
    geom_point(aes(x = x, y = y), data = data) +
    geom_line(aes(x = new_data$x, y = pred$Estimate), color = "red") +
    geom_ribbon(aes(x = new_data$x, ymin = pred$`Q2.5`, ymax = pred$`Q97.5`), alpha = 0.2, fill = "red") +
    labs(title = "Bayesian Linear Regression Prediction", x = "x", y = "y") +
    theme_minimal()
