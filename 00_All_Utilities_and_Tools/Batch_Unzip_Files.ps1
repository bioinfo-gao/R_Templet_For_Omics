# --- PowerShell 批量解压脚本 ---

# 设置要搜索的文件类型
$ZipExtension = "*.zip"

# Get-ChildItem -Recurse 会遍历当前目录及其所有子目录
Get-ChildItem -Path . -Recurse -Filter $ZipExtension | ForEach-Object {
    
    # 定义解压的目标路径：在 zip 文件所在的目录创建同名文件夹（不含 .zip 扩展名）
    # $_.BaseName 是文件名（不带扩展名）
    $Destination = Join-Path -Path $_.Directory.FullName -ChildPath $_.BaseName

    # --- 安全检查 ---
    # 检查目标文件夹是否存在，如果不存在则创建
    if (-not (Test-Path -Path $Destination)) {
        Write-Host "创建目录: $Destination" -ForegroundColor Green
        New-Item -Path $Destination -ItemType Directory | Out-Null
    }

    Write-Host "开始解压: $($_.Name) 到 $Destination" -ForegroundColor Yellow

    # --- 核心解压逻辑 ---
    # 使用 Windows Shell.Application COM 对象进行解压（内置方法）
    try {
        $shell = new-object -com shell.application
        $zipFile = $shell.NameSpace($_.FullName)
        $destinationFolder = $shell.NameSpace($Destination)
        
        # 复制文件（即执行解压操作）。参数 0 表示静默操作。
        $destinationFolder.CopyHere($zipFile.Items(), 0)

        Write-Host "解压完成: $($_.Name)" -ForegroundColor Cyan
    } catch {
        Write-Host "!!! 解压失败: $($_.Name). 错误: $($_.Exception.Message)" -ForegroundColor Red
    }
}

Write-Host "`n--- 所有 ZIP 文件的批量解压任务已完成 ---" -ForegroundColor White