# 1. 强制移除包
remove.packages(c("future", "future.apply"))

# 2. 重新安装
# 务必在有网络连接的情况下运行
install.packages(c("future", "future.apply"))

library(future)    # 用于并行计算
library(future.apply)    # 用于并行计算
library(limma)     # 用于数据处理，特别是 avereps() 函数
library(ggplot2)   # 用于绘图
library(ggpubr)    # 用于在图中添加相关性统计 (stat_cor)
library(ggExtra)   # 尽管未在代码中使用，但保留

gene="ARL15"                # 指定用于计算相关性的目标基因
corFilter=0.3               # 相关性系数的筛选阈值
pFilter=0.05                # P 值的筛选阈值

# 设置数据目录路径。注意：在 R 中使用正斜杠 (/) 是最佳实践。
data_dir = "C:/Users/zhen-/Code/R_code/R_For_DS_Omics/01.New_TCGA"

# 拼接完整的表达文件路径
expFile = file.path(data_dir, "combined_RNAseq_FPKM.txt")

# 检查文件路径是否正确
print(expFile) 

# 读取表达矩阵文件
# 假设第一列是基因名，第二列（可能）是基因 ID，后续是样本的 FPKM 值
rt1 = read.table(expFile, header=T, sep="\t", check.names=F, row.names = NULL)
# 将数据转换为矩阵格式
rt = as.matrix(rt1)

# 打印前两行和前五列进行检查
rt[1:2, 1:5]


# 设置行名为唯一基因名，假设基因名在第一列。
# 注意：如果文件第一列是基因名，直接用 row.names=1 更好，但为了保持原代码逻辑，我们使用 make.unique。
# 此处假定 rt[,1] 是基因名
rownames(rt) = make.unique(rt[,2]) 

# 提取表达数据。原代码注释提到 "not 2"，通常基因名是第1列。
# 故从第2列（或第3列，取决于文件结构）开始取表达值。这里沿用原代码的 rt[,3:ncol(rt)]。
# 如果第2列是基因 ID，则从第3列开始；如果第2列就是第一个样本，则应该从 rt[,2:ncol(rt)]。
exp = rt[,4:ncol(rt)] 
exp[1:2, 1:5]

# 转换为数值型矩阵
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)

# 使用 limma 包中的 avereps() 函数处理重复的行名（如果 make.unique 导致了重复行名，或者原始数据中存在重复）
data = avereps(data)

# 过滤掉表达量过低的基因（在所有样本中平均表达量大于 1）
data = data[rowMeans(data)>1,]
data[1:2, 1:6]

# **样本分组处理 (TCGA 特有)**
# 根据 TCGA 样本 ID (如 TCGA-XX-XXXX-01) 的第4段分割，提取第4段的第1个字符 (0/1/2)。
group = sapply(strsplit(colnames(data),"\\-"), "[", 4)
group = sapply(strsplit(group,""), "[", 1)

# 将数字 "2" 替换为 "1"。
# 0: 肿瘤 (Tumor), 1: 正常/癌旁 (Normal/Solid Tissue Normal), 2: 复发肿瘤（Rare）
group = gsub("2", "1", group)

# **关键逻辑错误修复：**
# 原代码 data=data[,group==0] 只保留了 '0' 组（肿瘤）的样本。
# 如果要计算基因间相关性，通常需要使用所有样本，或者仅使用肿瘤样本。
# **此处保留原代码逻辑：只保留肿瘤样本 (group == "0")**
data = data[,group=="0"] 

# 对表达数据进行 Log2(x+1) 转换，以使其更接近正态分布
data = log2(data+1)
data[1:2, 1:6]

# --- 相关性分析核心部分 ---

# 提取目标基因 (gene) 的表达向量
x = as.numeric(data[gene,])

# **重要修复：** 初始化用于存储结果的数据框，而不是空的 res3。
# 结果列表 `res_list` 应该在 `future_lapply` 之前初始化
# `future_lapply` 返回的是一个列表，而不是直接的 data.frame
res_list <- list() 

# 初始化用于存储相关性结果的 outTab 数据框
# 必须在循环或并行计算外部初始化
outTab = data.frame(Query=character(), Gene=character(), cor=numeric(), pvalue=numeric())


# 设置并行策略：使用多核/多会话
plan(multisession) 
# launch a number of background R sessions equal to the number of available cores detected by availableCores()

# 计时并使用 future_lapply 进行并行计算
system.time({
    # **关键修复：** future_lapply 的正确语法是：
    # future_lapply(X, FUN, ...)
    # X 是要迭代的向量/列表。这里是 data 的行名，即所有基因。
    # FUN 是应用于每个元素的函数。
    res_list <- future_lapply(rownames(data), FUN = function(j){
        # print(j) # 并行运行时，print 信息可能不会立即或按顺序显示

        # 跳过目标基因本身
        if(gene==j){
            return(NULL) # 并行计算中，跳过并返回 NULL
        }

        # 提取当前基因的表达向量
        y = as.numeric(data[j,])

        # 计算皮尔逊相关性
        corT = cor.test(x, y, method = 'pearson')
        cor = corT$estimate
        pvalue = corT$p.value

        # 构建当前基因的结果行。注意：在并行计算中，不应使用 rbind(outTab, ...)，
        # 因为 outTab 在各个会话中是独立的，应在函数内部返回单个结果。
        result_row = data.frame(
            Query = gene,
            Gene = j,
            cor = cor,
            pvalue = pvalue
        )

        # 绘制并保存相关性图
        # **重要：** 在并行环境中，文件写入（如 pdf/dev.off）可能导致问题或覆盖。
        # 推荐在非并行或后续串行步骤中进行绘图，但为保持原意，保留绘图代码。
        if((abs(cor)>corFilter) & (pvalue<pFilter)){
            df1 = as.data.frame(cbind(x,y))
            p1 = ggplot(df1, aes(x, y)) +
                xlab(paste0(gene, " expression"))+ ylab(paste0(j, " expression"))+
                geom_point()+
                geom_smooth(method="lm", formula=y~x) +
                theme_bw()+
                stat_cor(method = 'pearson', aes(x =x, y =y))

            # 使用 tryCatch 确保一个绘图错误不会中断整个并行进程
            tryCatch({
                pdf(file=paste0("cor.", j, ".pdf"), width=5, height=4.6)
                print(p1)
                dev.off()
            }, error = function(e) {
                cat("Error drawing plot for gene:", j, "\n")
            })
        }

        return(result_row)
    } )
} )


# # ... (代码前半部分不变，保持 plan(multisession)) ...
# 
# # 计时并使用 future_lapply 进行并行计算
# system.time({ 
#     res_list <- future_lapply(rownames(data), FUN = function(j){ 
#         
#         # 🌟 关键修复：确保在每个子进程中加载所需的包
#         library(limma)
#         library(ggplot2)
#         library(ggpubr)
#         
#         # ... (函数体内部其余代码不变) ...
#         
#         # 跳过目标基因本身
#         if(gene==j){
#             return(NULL)
#         } 
#         
#         # ... (后续计算和绘图代码) ...
#         
#         return(result_row)
#     } )
# } )

# ... (代码后半部分不变) ...


# 将所有并行计算返回的结果 (列表) 合并成一个数据框
# **重要修复：** 使用 res_list 替换原代码中的 res3
res3 = do.call(rbind, res_list)

# 将结果从列表转换为数据框后，才能进行后续的筛选和保存操作
# **重要修复：** 检查 res3 是否包含数据
if (!is.null(res3) && nrow(res3) > 0) {
    
    # 将最终结果写入文件
    write.table(file="corResult.txt", res3, sep="\t", quote=F, row.names=F)
    
    # 根据阈值筛选显著相关的基因
    outTab_sig = res3[abs(as.numeric(res3$cor)) > corFilter & as.numeric(res3$pvalue) < pFilter,]
    
    # 将显著相关的结果写入文件
    write.table(file="corSig.txt", outTab_sig, sep="\t", quote=F, row.names=F)
} else {
    print("Warning: No correlation results were generated.")
}

# 恢复默认的串行计划
plan(sequential)

