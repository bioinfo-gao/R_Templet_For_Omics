library(ggplot2)
library(dplyr)

# 1. 美学 (Aesthetics) Basic 
# color	离散（分类）或连续	点/线的颜色
# fill	离散（分类）或连续	区域的填充颜色（如条形图）
# size	连续	点的大小或线的粗细
# alpha	连续（0到1）	图形的透明度
# shape	离散（分类）	点的形状（0-25）
# linetype	离散（分类）	线的类型（实线、虚线、点线等）
# 使用内置的 mtcars 数据
data(mtcars)
mtcars$cyl_f <- factor(mtcars$cyl) # 将气缸数       转为因子（离散变量）
mtcars$am <- factor(mtcars$am)     # 将transmission 转为因子（离散变量） 
mtcars

ggplot(data = mtcars, 
       aes(x = wt, y = mpg,
           color = cyl_f,   # 映射：颜色由气缸数决定
           size = hp,       # 映射：点的大小由马力决定
           shape = am)) +   # 映射：点的形状由变速箱类型（am）决定 # error prone
    geom_point(alpha = 0.7)

# 设置要点： 你可以在 ggplot() 中设置全局映射，也可以在特定的 geom_ 图层中设置
# 局部映射或固定值。
# 局部映射： geom_point(aes(shape = gear)) - 仅对点图层生效。

# 固定值： geom_point(color = "blue") - 所有的点都固定为蓝色，不需要在 aes() 内部。

# 2. 坐标轴标度 (scale_*()) 的高级定制scale_ 函数用于控制如何将数据值转换为视觉效果。
# 它是定制图例、颜色、坐标轴刻度的核心。
# 2.1. 坐标轴标度 (scale_x_*, scale_y_*)函数描述关键参数
# scale_x_continuous() / 
# scale_y_continuous() 调整连续型轴name（轴标题），
# breaks（刻度位置），
# labels（刻度标签），limits（轴范围）
# scale_x_log10() / scale_y_log10()将轴转换为对数刻度breaks, labels
# 
ggplot(mtcars, aes(x = wt, y = mpg)) +
    geom_point() +
    scale_y_continuous(
        name = "每加仑英里数 (MPG)",         # 新标题
        limits = c(0, 35),                   # 强制 Y 轴范围
        breaks = seq(10, 30, by = 5),        # 自定义刻度线
        labels = paste0(seq(10, 30, by = 5), "L") # 自定义刻度标签
    )

#  颜色和填充标度 (scale_color_*, scale_fill_*)这是让图表美观的关键。
# 你需要根据你映射的数据类型（离散或连续）选择不同的 scale。
# 函数描述关键参数
# scale_color_manual()   /scale_fill_manual()手动指定每个类别的颜色values (颜色向量)，
# scale_color_discrete() /scale_fill_discrete()默认的离散颜色name, 
# scale_color_brewer()   / scale_fill_brewer()使用ColorBrewer 调色板（推荐）palette (调色板名称，如 "Set1", "Dark2")
# scale_color_gradient() / scale_fill_gradient()连续型数据颜色
# 渐变low, high (起始色和结束色)
# 
# 
# 3.1 使用 ColorBrewer 调色板
ggplot(mtcars, aes(x = wt, y = mpg, color = cyl_f)) + # 
    geom_point(size = 3) +
    scale_color_brewer(                                 # color 
        palette = "Set1", # 使用 Set1 调色板，非常适合离散数据
        name = "气缸数"
    )


# 4. 主题 (theme()) 的完全定制theme() 函数负责图表中所有不依赖于数据的元素：
# 字体、背景、网格线、图例位置、标题格式等。
# 3.1. 预设主题首先，你可以使用预设主题快速设置风格：
# theme_minimal(): # 干净、简洁
# theme_bw(): 黑白背景
# theme_classic(): 经典风格（无网格线）
# 
# 3.2. theme() 核心设置元素theme() 函数的参数结构是 element_ 函数。
# theme()参数           描述              对应的 element_ 函数
# plot.title         图表主标题           element_text()
# axis.text.x        X 轴刻度标签         element_text()
# panel.grid.major   主要网格线           element_line()
# panel.background   绘图区域背景         element_rect()
# legend.position    图例位置             字符串 ("top", "bottom", "none") 或坐标向量 c(x, y)
# 
# 

final_plot <- ggplot(mtcars, aes(x = wt, y = mpg, color = cyl_f)) +
    geom_point(size = 3) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = "汽车重量与油耗分析",
         subtitle = "按气缸数着色",
         caption = "数据来源：mtcars") +
    
    # --- 自定义主题开始 ---
    theme_minimal(base_size = 14) + # 从简洁主题开始，并设置基础字体大小
    theme(
        # 1. 标题设置
        plot.title = element_text(
            size = 20,              # 字体大小
            face = "bold",          # 字体加粗
            hjust = 0               # 靠左对齐 (0=左，0.5=中，1=右)
        ),
        # 2. 坐标轴文本设置
        axis.text.x = element_text(
            angle = 45,             # X 轴标签旋转45度
            vjust = 1,              # 垂直对齐底部
            hjust = 1               # 水平对齐右侧
        ),
        # 3. 网格线设置
        panel.grid.minor = element_blank(), # 移除次要网格线
        panel.grid.major = element_line(
            color = "grey85",
            linewidth = 0.5
        ),
        # 4. 图例位置
        legend.position = "bottom", # 将图例移到下方
        legend.title = element_blank(), # 隐藏图例标题
        # 5. 边框
        panel.border = element_rect(
            colour = "black", 
            fill = NA, 
            linewidth = 1
        )
    )
# --- 自定义主题结束 ---

print(final_plot)


# 5 图例（Legend）定制
# 图例的定制主要通过 theme() 和 guides() 函数来完成。
# 
# theme(legend.position = ...): 控制图例的整体位置。
# 
# guides(): 对特定的美学（如 color, size）进行精细调整。
# 
final_plot + 
    guides(
        color = guide_legend(
            #title.position = "bottom", # 图例标题放在顶部
            title.position = "right", # 图例标题放在顶部
            keywidth = 20,           # 键（Key）的宽度
            keyheight = 1,          # 键的高度
            nrow = 1                # 强制图例为一行
        ),
        size = "none" # 隐藏 size 的图例
    )

print(final_plot)

# 4. 强大的辅助函数最后，一些常用的辅助函数可以极大地提高效率：函数描述labs()快速设置标题 (title, subtitle, x, y, caption)facet_wrap() / facet_grid()分面：根据一个或多个分类变量，将图表拆分为多个子图coord_flip()交换 X 轴和 Y 轴（常用于条形图）scale_y_reverse()反转 Y 轴方向
# 增加一个分类变量，根据气缸数和变速箱类型分面
ggplot(mtcars, aes(x = hp, y = qsec)) +
    geom_point() +
    # 按 cyl_f 分面，分成两行，自由缩放 Y 轴 (scales = "free_y")
    facet_wrap(~ cyl_f + am, 
               ncol = 2, 
               scales = "free_y") + 
    labs(title = "分面图：不同气缸/变速箱下的性能")

# 5 对于连续变量，更合适的视觉映射是 颜色 (color) 或 大小 (size)，而不是形状。
# 尝试使用 size 来映射连续变量（如 hp 马力）
ggplot(data = mtcars, aes(x = wt, y = mpg, size = hp)) +
    geom_point()

# 或者使用 color 来映射连续变量（如 hp 马力），这会得到一个颜色渐变图例
ggplot(data = mtcars, aes(x = wt, y = mpg, color = hp)) +
    geom_point()
mtcars
# 或者使用 color 来映射连续变量（如 hp 马力），这会得到一个颜色渐变图例
ggplot(data = mtcars, aes(x = wt, y = mpg, size = hp, color = disp)) +
    geom_point()


# 方案三：使用 scale_shape_binned()（高级/分箱）
# 如果您坚持要使用形状来映射连续变量，您必须首先将这个连续变量分箱 (binning)，也就是把它分成几个离散的区间。ggplot2 为此提供了专门的函数。
# 
# R

# 使用 scale_shape_binned() 自动将连续变量 hp 分成几个区间，并映射到形状
ggplot(data = mtcars, aes(x = wt, y = mpg, shape = hp)) +
    geom_point() +
    scale_shape_binned() # <--- 必须添加这个标度函数
