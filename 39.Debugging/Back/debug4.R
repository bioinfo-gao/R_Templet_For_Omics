# 定义一个需要调试的函数
my_function <- function(x) {
    intermediate <- x * 2
    browser()  # 在这里进入调试模式
    result <- intermediate + 10
    return(result)
}

# 调用函数
my_function(5)