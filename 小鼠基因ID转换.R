library(org.Mm.eg.db)

# 读取数据（文件路径为 D:\GSE264139_genes.readcount.txt）
data <- read.table("GSE264139_genes.readcount.txt", header = TRUE, stringsAsFactors = FALSE)

# 提取第一列的 Ensembl ID
ensembl_ids <- data[, 1]

# 使用 org.Mm.eg.db 将 Ensembl ID 转换为 Gene Symbol
gene_symbols <- mapIds(
  org.Mm.eg.db,              # 使用的注释包
  keys = ensembl_ids,        # 输入的 Ensembl ID
  keytype = "ENSEMBL",       # 输入 ID 的类型
  column = "SYMBOL",         # 输出的注释信息（Gene Symbol）
  multiVals = "first"        # 如果有多个匹配，取第一个
)

# 将第一列替换为 Gene Symbol
data[, 1] <- gene_symbols

# 检查是否有未匹配的 Ensembl ID（显示为 NA）
if (any(is.na(gene_symbols))) {
  warning("部分 Ensembl ID 未找到对应的 Gene Symbol，已标记为 NA。")
}

# 查看转换后的数据
head(data)

# 保存结果到新文件
write.table(data, file = "D:\\GSE264139_genes_with_symbols.txt", sep = "\t", quote = FALSE, row.names = FALSE)
