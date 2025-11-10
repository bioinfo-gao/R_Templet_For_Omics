# 1. 定义一个新的、全英文的R主目录
new_home_dir <- "C:/R_User_Home" 

# 2. 创建该目录
if (!dir.exists(new_home_dir)) {
    dir.create(new_home_dir, recursive = TRUE)
}

# 3. 设置R的环境变量 R_USER，将其指向新的主目录
Sys.setenv(R_USER = new_home_dir)
Sys.setenv(R_LIBS_USER = file.path(new_home_dir, "R_Library")) # 确保包安装在新目录

message(paste("R 的用户主目录已设置为:", new_home_dir))

# 4. 再次设置编译路径（确保完整性）
Sys.setenv(TMPDIR = new_home_dir)
Sys.setenv(TEMP = new_home_dir)

