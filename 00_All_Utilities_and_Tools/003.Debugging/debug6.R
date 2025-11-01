# https://www.kimi.com/chat/cui4ls1djjppps7pg2gg
# 1. 定义一个故意留 bug 的函数
add_ratio <- function(a, b) {
    #browser()                       # <-- 调试入口
    num <- a + b
    den <- b - a                    # 若 a == b 会除 0
    ratio <- num / den
    return(ratio)
}

# 2. 调用函数（会触发 browser）
add_ratio(3, 3)


debugonce(add_ratio)   # 仅下次调用进入调试
add_ratio(5, 2)

# 总结口诀
# “先 browser 或断点，
# 再 Source；
# F10 逐行看变量，(n)
# Shift+F4 进子函数；
# Shift+F5 继续跑，
# Shift+F8 随时停。”
# q or Q quit debugger
Q