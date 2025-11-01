# 目标：
# 1. 计算一个向量中所有“正数”的累积和。
# 2. 返回一个列表，包含最终的总和以及每一步的累积和向量。

process_data <- function(data_vector) {
    
    cumulative_sum <- 0
    # 我们预先分配一个向量来存储结果
    processed_vector <- numeric(length(data_vector)) 
    
    # 在这里设置一个断点
    for (i in 1:length(data_vector)) {
        
        current_val <- data_vector[i]
        
        # 🐛 Bug 1 隐藏在这里
        if (current_val > 0) {
            cumulative_sum <- cumulative_sum + current_val
        }
        
        # 🐛 Bug 2 隐藏在这里 (逻辑错误)
        processed_vector[i] <- current_val 
    }
    
    return(list(total = cumulative_sum, vec = processed_vector))
}


# --- 运行我们的代码 ---

# 准备一个包含“脏数据”的向量
# 它包含一个负数和一个会导致崩溃的字符串
#input_data <- c(10, 20, -5, 15, "oops", 30) 
input_data <- c(10, 20, -5, 15, 30) 

# 运行函数
result <- process_data(input_data)

print(result)

