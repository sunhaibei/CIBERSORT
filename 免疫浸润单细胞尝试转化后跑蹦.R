library(Matrix)
library(data.table)

load_sample_data <- function(sample_id) {
  # 构造文件路径
  barcodes_file <- paste0(sample_id, "_barcodes.tsv.gz")
  features_file <- paste0(sample_id, "_features.tsv.gz")
  matrix_file <- paste0(sample_id, "_matrix.mtx.gz")
  
  # 读取 barcodes 和 features
  barcodes <- fread(barcodes_file, header = FALSE)$V1
  features <- fread(features_file, header = FALSE)$V1
  
  # 读取稀疏矩阵
  matrix <- readMM(matrix_file)
  
  # 设置行名和列名
  rownames(matrix) <- features
  colnames(matrix) <- barcodes
  
  return(matrix)
}

samples <- c('GSM5333784_HC1', 'GSM5333785_HC2', 'GSM5333786_Sep1', 'GSM5333787_Sep2', 
             'GSM5333788_S20', 'GSM5333789_S30', 'GSM5333790_H20', 'GSM5333791_H27', 'GSM5333792_H30')

# 加载所有样本的数据
data <- lapply(samples, load_sample_data)
names(data) <- samples

# 合并所有样本的矩阵
combined_matrix <- do.call(cbind, data)

# 设置列名（样本名称 + 细胞条形码）
colnames(combined_matrix) <- paste(rep(names(data), sapply(data, ncol)), unlist(lapply(data, colnames)), sep = "_")


# 保存合并后的矩阵
writeMM(combined_matrix, 'combined_expression_matrix.mtx')
# 保存基因名称和样本名称
write.table(rownames(combined_matrix), 'gene_names.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(combined_matrix), 'sample_names.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)
library(Matrix)
library(data.table)

# 加载稀疏矩阵
matrix <- readMM('combined_expression_matrix.mtx')

# 加载基因名称和样本名称
gene_names <- read.table('gene_names.tsv', header = FALSE, stringsAsFactors = FALSE)[[1]]
sample_names <- read.table('sample_names.tsv', header = FALSE, stringsAsFactors = FALSE)[[1]]

# 设置行名和列名
rownames(matrix) <- gene_names
colnames(matrix) <- sample_names

# 定义分块大小（例如每次处理 1000 行）
chunk_size <- 1000
num_genes <- nrow(matrix)

# 创建一个空的输出文件，并写入表头
header <- c("Gene", sample_names)
write.table(t(header), 'gene_expression_matrix.csv', sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 分块处理并写入文件(后面跑崩了，别试)
for (start_row in seq(1, num_genes, chunk_size)) {
  end_row <- min(start_row + chunk_size - 1, num_genes)
  chunk <- matrix[start_row:end_row, ]
  
  # 提取非零值
  non_zero_values <- summary(chunk)
  
  # 转换为数据框
  chunk_df <- data.frame(
    Gene = gene_names[start_row:end_row][non_zero_values$i],
    Sample = sample_names[non_zero_values$j],
    Value = non_zero_values$x
  )
  
  # 将数据框转换为宽格式
  chunk_wide <- dcast(chunk_df, Gene ~ Sample, value.var = "Value", fill = 0)
  
  # 追加写入文件
  write.table(chunk_wide, 'gene_expression_matrix.csv', sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}
write.table(chunk_wide, 'gene_expression_matrix.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
