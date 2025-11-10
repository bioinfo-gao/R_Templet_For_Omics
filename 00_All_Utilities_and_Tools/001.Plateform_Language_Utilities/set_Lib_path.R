# 1. 设置新的库路径变量
new_lib_path <- normalizePath("~", mustWork = TRUE) # 获取用户主目录
new_lib <- file.path(new_lib_path, "R_Library") 
new_lib
# 2. 创建目录 (如果不存在)
if (!dir.exists(new_lib)) {
    dir.create(new_lib, recursive = TRUE)
}

# 3. 将新路径写入 .Rprofile 文件，使设置永久生效
dotRprofile <- file.path(new_lib_path, ".Rprofile")

# 检查 .Rprofile 是否存在，如果不存在则创建
if (!file.exists(dotRprofile)) {
    file.create(dotRprofile)
}

# 写入新的路径配置到 .Rprofile (注意使用 'append = TRUE' 或覆盖)
# 我们使用 readLines 和 writeLines 来确保我们只添加，不破坏现有内容
current_lines <- readLines(dotRprofile)
new_line <- paste0('.libPaths(c("', new_lib, '", .libPaths()))')

# 检查是否已存在，防止重复添加
if (!any(grepl(".libPaths", current_lines))) {
    write(new_line, file = dotRprofile, append = TRUE)
}

message(paste("新的R包安装路径已设置为:", new_lib))
