# 永久更改（设置默认启动目录）：
# 进入 RStudio 的顶部菜单，选择 工具（Tools） > 全局选项（Global Options）。
# 在弹出的窗口中，找到“默认工作目录（Default working directory）”选项。
# 点击 浏览（Browse），选择一个你希望 RStudio 每次启动时都作为主目录的文件夹。
# 点击“确定”保存设置。下次启动 RStudio 时，它就会自动使用新的默认目录

# C:/Users/zhen-/AppData/Local/R/win-library/4.5/maftools


.libPaths()

# [1] "C:/Users/zhen-/AppData/Local/R/win-library/4.5"
# [2] "C:/Program Files/R/R-4.5.1/library"     

sessionInfo()

find.package("maftools")
# [1] "C:/Users/zhen-/AppData/Local/R/win-library/4.5/maftools"


# =============================  设置用户家目录

# 这一步是找出您的用户主目录，通常是 C:\Users\zhen-
normalizePath("~")
# 运行以下命令，创建一个 .Rprofile 文件并打开它进行编辑
file.edit("~/.Rprofile")
