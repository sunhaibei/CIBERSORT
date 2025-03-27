# 读取数据
data <- read.table("symbols.txt", header = TRUE, stringsAsFactors = FALSE)

# 1. 删除 NA 值
# 检查是否有 NA 值
print("NA 值检查：")
print(colSums(is.na(data)))

# 删除 geneID 列为 NA 的行
data <- data[!is.na(data$geneID), ]

# 2. 删除 geneID 重复列
# 检查重复的 geneID
duplicated_genes <- duplicated(data$geneID)
print("重复的 geneID：")
print(data$geneID[duplicated_genes])

# 保留第一个出现的 geneID，删除后续重复的行
data <- data[!duplicated(data$geneID), ]

# 查看处理后的数据
print("处理后的数据：")
print(head(data))
# 3. 将 geneID 列的小写字母转换为大写字母
data$geneID <- toupper(data$geneID)

# 4. 保存处理后的数据到新的 .txt 文件
write.table(data, "processed_symbol.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 5. 查看处理后的数据
print("处理后的数据：")
print(head(data))
