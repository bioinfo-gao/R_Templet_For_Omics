# 1. 安装包（如果这是你第一次使用它们）
# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("dplyr")
# install.packages("gapminder")

# 2. 加载这些包
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gapminder)

# 1. 查看 gapminder 数据的结构
str(gapminder)
# 'data.frame':    1704 obs. of  6 variables:
# $ country  : Factor w/ 142 levels "Afghanistan",..: 1 2 3 4 5 6 7 8 9 10 ...
# $ continent: Factor w/ 5 levels "Africa","Americas",..: 3 3 3 3 3 3 3 3 3 3 ...
# $ year     : int  1952 1957 1962 1967 1972 1977 1982 1987 1992 1997 ...
# $ lifeExp  : num  28.8 30.3 32 34 36.1 ...
# $ pop      : int  8425333 9240934 10267083 11537966 13079460 14880372 12881816 13867957 16317921 22227415 ...
# $ gdpPercap: num  779 821 853 836 740 ...

# 2. 我们只对2007年的数据感兴趣
# 使用 dplyr 包中的 filter() 函数来筛选数据
data_2007 <- filter(gapminder, year == 2007)

# 3. 查看筛选后的数据
head(data_2007)
#       country continent year lifeExp      pop gdpPercap
# 1 Afghanistan      Asia 2007  43.828 31889923   974.5803
# 2     Albania    Europe 2007  76.423  3600523  5937.0295
# 3     Algeria    Africa 2007  72.301 33333216  6223.3675
# 4      Angola    Africa 2007  42.731 12420476  4797.2313
# 5   Argentina  Americas 2007  75.320 40301927 12779.3796
# 6   Australia   Oceania 2007  81.235 20434176 34435.3674
# 
# 
# # 1. 设置画布：
#    data = data_2007 (我们刚创建的数据)
#    aes() 映射：x 轴是 gdpPercap, y 轴是 lifeExp
p1 <- ggplot(data = data_2007, aes(x = gdpPercap, y = lifeExp)) +
    # 2. 添加图层：   #    geom_point() 告诉 ggplot 画点 <<======= 
    geom_point()

# 3. 显示图表
print(p1)

p_with_bad_labels <- ggplot(data = data_2007, aes(x = gdpPercap, y = lifeExp, label = country)) +
    geom_point() +
    geom_text()     # 警告：这将是一场灾难！

print(p_with_bad_labels)


# 1. 映射保持不变 (x, y, label)
p_with_good_labels <- ggplot(data = data_2007, aes(x = gdpPercap, y = lifeExp, label = country)) +
    geom_point() +
    geom_text_repel()

# 3. 显示图表
print(p_with_good_labels)

getwd()
setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/30.Plotting/ggplot2_and_ggrepel")

p_themed <- ggplot(data = data_2007, 
                   aes(x = gdpPercap, y = lifeExp, label = country)) +
    geom_point() +
    geom_text_repel() +
    scale_x_log10() +                      # 1. 将 X 轴转换为 Log10 刻度
    labs(                                  # 2. 添加标题、副标题和坐标轴标签
        title = "2007年各国健康与财富",
        subtitle = "人均GDP vs. 预期寿命",
        x = "人均GDP (对数刻度)",
        y = "预期寿命 (岁)",
        caption = "数据来源: Gapminder"
    )

print(p_themed)
# earlier for contrast print(p_with_good_labels)


# 一张好的图表会用多种视觉元素传递信息。让我们用 continent（大洲）来控制颜色，用 pop（人口）来控制点的大小。
p_colorful <- ggplot(data = data_2007, 
                     aes(x = gdpPercap, 
                         y = lifeExp, 
                         label = country, 
                         color = continent,  # <-- 新增：颜色
                         size = pop)) +       # <-- 新增：大小
    geom_point(alpha = 0.6) + # (alpha = 0.6 增加透明度，防止点重叠)
    geom_text_repel() +
    scale_x_log10() +
    #scale_size(name = "人口",                       # 新增：自定义大小图例 ===>> 右下的legend
    #           labels = scales::comma) + # 使用 scales 包格式化数字（例如 1,000,000）
    scale_color_discrete(name = "大洲") +     # 新增：自定义颜色图例
    labs(
        title = "2007年各国健康与财富",
        #subtitle = "人均GDP vs. 预期寿命",
        x = "人均GDP (对数刻度)",
        y = "预期寿命 (岁)",
        caption = "数据来源: Gapminder"
    )

print(p_colorful)

## 风格统一

p_repel_custom <- ggplot(data = data_2007, 
                         aes(x = gdpPercap, y = lifeExp, 
                             label = country, color = continent, size = pop)) +
    
    geom_point(alpha = 0.6) +
    
    # --- 定制 ggrepel ---
    geom_text_repel(
        size = 2.5,                 # 调小标签字体大小 << ======= MOST IMPORTANT 
        color = "black",            # 强制标签为黑色（否则会跟随 color=continent 映射）
        max.overlaps = 15,          # 提高重叠的“容忍度”（允许更多标签显示）
        min.segment.length = 0,     # 总是绘制线段（即使很短）
        box.padding = 0.1,          # 标签文本框的内边距
        segment.color = 'grey50',   # 线段颜色
        segment.alpha = 0.5         # 线段透明度
    ) +
    # --- 结束定制 ---
    
    scale_x_log10(labels = scales::dollar) + # 使用 scales 包将 x 轴格式化为美元 $
    scale_size(name = "人口", labels = scales::comma) +
    scale_color_discrete(name = "大洲") +
    
    # 换一个更干净的主题
    theme_minimal(base_size = 14) + # (base_size 调整所有字体大小)
    
    # 调整主题：把图例移到下方
    theme(legend.position = "bottom") + 
    
    labs(
        title = "2007年各国健康与财富",
        subtitle = "人均GDP vs. 预期寿命",
        x = "人均GDP (对数刻度)",
        y = "预期寿命 (岁)",
        caption = "数据来源: Gapminder"
    )

print(p_repel_custom)


# 第5步：高级技巧 - 只标记特定点
# 在 142 个国家中，我们可能并不想标记所有国家，也许只想标记“金砖四国”（巴西、俄罗斯、印度、中国）和美国、日本。
# 
# ggplot2 的核心技巧：你可以为不同的 geom 层指定不同的数据！
# 
# # 1. 准备一个我们感兴趣的国家的列表
highlight_countries <- c("China", "India", "United States", "Japan", "Brazil", "Russia")

# 2. 从 data_2007 中筛选出这些国家的数据
highlight_data <- filter(data_2007, country %in% highlight_countries)

# 3. 绘制！
p_highlight <- ggplot(data = data_2007, # <-- 基础图层使用“所有数据”
                      aes(x = gdpPercap, y = lifeExp, 
                          color = continent, size = pop)) +
    
    # 图层 1：绘制“所有”的点（设为半透明）
    geom_point(alpha = 0.3) +
    
    # 图层 2：只为“高亮数据”绘制标签
    geom_text_repel(
        data = highlight_data,      # <-- 关键！为这层指定“高亮数据”
        aes(label = country),       # 标签也来自高亮数据
        color = "black",            # 强制黑色
        size = 4,                   # 字体大一点
        segment.color = "black"
    ) +
    
    # 图层 3：为“高亮数据”再画一层点（更不透明），让它们凸显出来
    geom_point(data = highlight_data, 
               shape = 21,          # 使用带边框的形状
               color = "black",     # 边框为黑色
               size = 5,            # 大一点
               stroke = 1)          # 边框粗细

# ... (其他 scale 和 theme 设置保持不变) ...
scale_x_log10(labels = scales::dollar) + 
    scale_size(name = "人口", labels = scales::comma, range = c(1, 15)) + # (range 调整点大小范围)
    scale_color_discrete(name = "大洲") +
    theme_minimal(base_size = 14) + 
    theme(legend.position = "bottom") + 
    guides(color = guide_legend(override.aes = list(size=5)), # 增大图例中的点
           size = "none") + # 隐藏人口图例，因为它现在会误导
    labs(
        title = "2007年各国健康与财富",
        subtitle = "重点突出显示金砖四国、美国和日本",
        x = "人均GDP (对数刻G刻度)",
        y = "预期寿命 (岁)",
        caption = "数据来源: Gapminder"
    )

print(p_highlight)
