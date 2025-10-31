# ----------------------------------------------------
# 确保所有必需的并行包都被加载
library(future)       
library(future.apply) # 🌟 关键：显式加载 future.apply
# ----------------------------------------------------

# 设置并行策略
plan(multisession)

# 临时测试 future_lapply (确认现在可以找到函数)
test_list <- list(1, 2, 3)
test_result <- future_lapply(test_list, FUN = function(x) x * 2)
print(test_result)
