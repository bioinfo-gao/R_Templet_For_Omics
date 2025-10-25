# brms: Bayesian Regression Models using 'Stan'
# Stan 的名称来源于 Stanisław (Stan) Ulam，一位波兰裔美国数学家（1909-1984），他与 Nicholas Metropolis 共同开发了蒙特卡洛方法（Monte Carlo method），这是 Stan 平台核心算法（如 Hamiltonian Monte Carlo，HMC）的基础。命名“Stan”是为了向 Ulam 的贡献致敬，反映了平台在统计推断和概率建模中的技术根源。
library(rstan)

# Stan 模型代码
stan_code <- "
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  y ~ normal(alpha + beta * x, sigma);
}
"

# 数据
data <- list(N = 100, x = rnorm(100), y = 2 + 3 * rnorm(100) + rnorm(100))
data
str(data)

# 拟合模型
fit <- stan(model_code = stan_code, data = data, iter = 1000, chains = 2)
print(fit)
