# 这个脚本会列出仓库中前 10 个最大的文件，包括其在历史记录中的大小
git rev-list --all --objects | \
  git cat-file --batch-check='%(size) %(objectname) %(rest)' | \
  sed -n 's/^[^ ]* //p' | \
  sort -n -r | \
  head -n 10


  # 适用于旧版本 Git 寻找大文件的命令
git rev-list --objects --all | grep "$(git verify-pack -v .git/objects/pack/*.idx | sort -k 3 -n | tail -5 | awk '{print $1}')"