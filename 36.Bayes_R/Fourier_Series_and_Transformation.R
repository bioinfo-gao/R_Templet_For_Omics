# 周期方波的傅立叶级数
T <- 2  # 周期
n_terms <- 10  # 级数项数
t <- seq(-T, T, length.out = 1000)
f <- function(t) ifelse(t %% T < T/2, 1, -1)  # 方波

# 计算傅立叶级数
fourier_series <- function(t, n_terms, T) {
    result <- rep(0, length(t))
    for (n in 1:n_terms) {
        bn <- 4 / (pi * n) * sin(n * pi / 2)  # 方波的 b_n 系数
        result <- result + bn * sin(2 * pi * n * t / T)
    }
    result
}

# 计算并绘制
y <- fourier_series(t, n_terms, T)
plot(t, f(t), type = "l", col = "black", lwd = 2, main = "Fourier Series of Square Wave")
lines(t, y, col = "red", lwd = 1)
legend("topright", legend = c("Original", "Fourier Approx"), col = c("black", "red"), lwd = 2)

###### ==============================================
# 傅立叶变换（非周期高斯脉冲）
library(ggplot2)
t <- seq(-10, 10, length.out = 1000)
f <- exp(-t^2)  # 高斯脉冲

# 计算傅立叶变换（使用 fft）
n <- length(t)
dt <- t[2] - t[1]
freq <- seq(-0.5/dt, 0.5/dt, length.out = n)
fft_result <- fft(f) / n
magnitude <- abs(fft_result)

# 绘制
ggplot(data.frame(freq = freq, mag = magnitude), aes(x = freq, y = mag)) +
    geom_line(color = "#3B82F6") +
    labs(title = "Fourier Transform of Gaussian Pulse", x = "Frequency", y = "Magnitude") +
    theme_minimal()