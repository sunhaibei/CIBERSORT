# 加载库
library(Seurat)
library(Matrix)
library(dplyr)
library(data.table)
library(preprocessCore)

# 设置你的文件路径（修改成你的数据路径）
data_dir <- "D:/GSE175453_RAW (1)/GSM533787"# 你的数据文件夹路径

# 读取 10X Genomics 数据
seurat_obj <- Read10X(data.dir = data_dir)

# 转换成 Seurat 对象
seurat <- CreateSeuratObject(counts = seurat_obj)

# 计算基因表达矩阵（取 logTPM 以适配 CIBERSORT）
bulk_expr <- as.data.frame(as.matrix(seurat@assays$RNA@counts))
bulk_expr <- log1p(bulk_expr)  # 取 log(TPM+1)
slotNames(seurat@assays$RNA)
bulk_expr <- as.data.frame(as.matrix(seurat@assays$RNA@layers$counts))
names(seurat@assays$RNA@layers)
bulk_expr <- as.data.frame(as.matrix(seurat@assays$RNA@layers[["counts.Gene Expression"]]))
library(Matrix)
# 检查 Assay5 中的所有插槽
slotNames(seurat@assays$RNA)
# 查看 layers 中的所有内容
names(seurat@assays$RNA$layers)

# 将 counts 转换为稀疏矩阵
counts_matrix <- seurat@assays$RNA@counts
counts_sparse <- as(counts_matrix, "CsparseMatrix")

# 然后可以将稀疏矩阵转换为数据框
bulk_expr <- as.data.frame(as.matrix(counts_sparse))



# 转换成 CIBERSORT 需要的格式
bulk_expr$GeneSymbol <- rownames(bulk_expr)
bulk_expr <- bulk_expr %>% select(GeneSymbol, everything())

# 保存数据以供 CIBERSORT 使用
write.table(bulk_expr, "bulk_expr_for_CIBERSORT.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 运行 CIBERSORT
source("CIBERSORT.R")  # 你需要下载 CIBERSORT.R 脚本并放在当前目录
lm22_file <- "LM22.txt"  # 免疫细胞参考基因集（CIBERSORT 官方提供）

# 运行 CIBERSORT
cibersort_res <- CIBERSORT(lm22_file, "bulk_expr_for_CIBERSORT.txt", perm=1000)

# 保存 CIBERSORT 结果
write.table(cibersort_res, "CIBERSORT_results.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# 画图可视化结果
library(ggplot2)
cibersort_res <- read.table("CIBERSORT_results.txt", header = TRUE, sep = "\t", row.names = 1)

# 转换为长格式
cibersort_long <- melt(cibersort_res, id.vars = "Mixture")

# 画图
ggplot(cibersort_long, aes(x = variable, y = value, fill = Mixture)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "免疫细胞类型", y = "细胞比例", title = "免疫细胞浸润分析")