#install.packages("Matrix")
#install.packages("tidyverse")
library(Matrix)
library(tidyverse)
# 读取 features 文件
features <- read.table(gzfile("features.tsv"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(features) <- c("Ensembl_ID", "Gene_Symbol", "Gene_Type")

# 读取 barcodes 文件
barcodes <- read.table(gzfile("barcodes.tsv.gz"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(barcodes) <- "Barcode"

# 读取 matrix 文件
matrix_data <- readMM(gzfile("matrix.mtx"))
# 将稀疏矩阵转换为稠密矩阵
expression_matrix <- as.matrix(matrix_data)
# 直接使用稀疏矩阵
expression_matrix <- matrix_data
# 设置行名和列名
rownames(expression_matrix) <- features$Gene_Symbol
colnames(expression_matrix) <- barcodes$Barcode

# 去除重复基因（如果有）
expression_matrix <- expression_matrix[!duplicated(features$Gene_Symbol), ]
# 保存为制表符分隔的文本文件
write.table(expression_matrix, "gene_expression_matrix.txt", sep = "\t", quote = FALSE)
getwd()
head(expression_matrix)

# 计算每个基因的总表达量
gene_sums <- rowSums(expression_matrix)

# 过滤掉低表达基因（例如总表达量小于 10）
keep_genes <- gene_sums >= 10
expression_matrix <- expression_matrix[keep_genes, ]

# 计算每个细胞的总表达量
cell_sums <- colSums(expression_matrix)

# 过滤掉低质量细胞（例如总表达量小于 100）
keep_cells <- cell_sums >= 100
expression_matrix <- expression_matrix[, keep_cells]


# 保存过滤后的稀疏矩阵
writeMM(expression_matrix, "filtered_expression_matrix.mtx")

# 转换为稠密矩阵
#dense_matrix <- as.matrix(expression_matrix)
rm(dense_matrix)
# 保存为文本文件
#write.table(dense_matrix, "filtered_expression_matrix.txt", sep = "\t", quote = FALSE)