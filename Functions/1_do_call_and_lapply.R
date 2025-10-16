# lapply：想象成 "List APPLY" - 对列表的每个元素应用函数，作用在list的每个元素上
 

# do.call：想象成 "DO a function CALL" - 执行一个函数调用, 作用在list整体上


# 总结：lapply 用于数据转换，do.call 用于函数调用。理解这个核心区别就能正确选择使用哪个函数了。

# 创建一个列表
my_list <- list(1, 2, 3)

# 使用 lapply 将每个元素乘以2
result <- lapply(my_list, function(x) x * 2)
# 结果为: [[1]] 2, [[2]] 4, [[3]] 6



lapply(iris, class)
Map(class, iris)
 
# 创建一个列表作为参数
args_list <- list(1, 2)

# 使用 do.call 调用函数
result <- do.call(sum, args_list)
# 结果为: 3
# 
# 
x <- lapply(iris, class)
y= do.call(c, x)

y
