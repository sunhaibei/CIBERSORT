

library(Matrix)
library(readr)

load_mtx_data <- function(sample_id) {
  # 构造文件路径
  barcodes_file <- paste0(sample_id, "_barcodes.tsv.gz")
  features_file <- paste0(sample_id, "_features.tsv.gz")
  matrix_file <- paste0(sample_id, "_matrix.mtx.gz")
  
  # 读取 barcodes 和 features
  barcodes <- read_tsv(barcodes_file, col_names = FALSE, col_types = cols())
  features <- read_tsv(features_file, col_names = FALSE, col_types = cols())
  
  # 读取稀疏矩阵
  matrix <- readMM(matrix_file)
  
  # 设置行名和列名
  rownames(matrix) <- features$X1  # 行名为基因名称
  colnames(matrix) <- barcodes$X1  # 列名为条形码
  
  return(matrix)
}
samples <- c('GSM5333784_HC1', 'GSM5333785_HC2', 'GSM5333786_Sep1', 'GSM5333787_Sep2', 
             'GSM5333788_S20', 'GSM5333789_S30', 'GSM5333790_H20', 'GSM5333791_H27', 'GSM5333792_H30')

# 加载所有样本的数据
data <- lapply(samples, load_mtx_data)
names(data) <- samples
# 合并所有样本的稀疏矩阵
combined_matrix <- do.call(cbind, data)

# 设置列名
colnames(combined_matrix) <- paste(rep(names(data), sapply(data, ncol)), unlist(lapply(data, colnames)), sep = "_")

# 保存为 .mtx 文件
writeMM(combined_matrix, 'combined_expression_matrix.mtx')

# 保存基因名称和样本名称
write.table(rownames(combined_matrix), 'gene_names.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(combined_matrix), 'sample_names.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)

library(Matrix)
library(data.table)

# 加载稀疏矩阵、基因名称和样本名称
matrix <- readMM('combined_expression_matrix.mtx')
gene_names <- read.table('gene_names.tsv', header = FALSE, stringsAsFactors = FALSE)[[1]]
sample_names <- read.table('sample_names.tsv', header = FALSE, stringsAsFactors = FALSE)[[1]]

# 设置行名和列名
rownames(matrix) <- gene_names
colnames(matrix) <- sample_names

# 定义分块大小（例如每次处理 100 列）
chunk_size <- 100
num_samples <- ncol(matrix)

# 创建一个空的输出文件，并写入表头
header <- c("Gene", sample_names)
write.table(t(header), 'cibersortx_input.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 分块处理并写入文件
for (start_col in seq(1, num_samples, chunk_size)) {
  end_col <- min(start_col + chunk_size - 1, num_samples)
  chunk <- matrix[, start_col:end_col]
  
  # 转换为数据框并添加基因名称列
  chunk_df <- as.data.frame(as.matrix(chunk))
  chunk_df <- cbind(Gene = rownames(chunk_df), chunk_df)
  
  # 追加写入文件
  write.table(chunk_df, 'cibersortx_input.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}












