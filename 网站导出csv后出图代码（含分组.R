# 安装并加载必要的库
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(pheatmap)) install.packages("pheatmap")

library(tidyverse)
library(ggplot2)
library(pheatmap)
# 读取 CIBERSORTx 结果
cibersort_results <- read.csv("CIBERSORTx_Job6_Results.csv", check.names = FALSE)


# 添加分组信息
cibersort_results <- cibersort_results %>%
  mutate(Group = case_when(
    Mixture %in% c("A1", "A2", "A3") ~ "Group A",
    Mixture %in% c("D1", "D2", "D3") ~ "Group B",
    Mixture %in% c("C1", "C2", "C3") ~ "Group C",
    TRUE ~ NA_character_  # 如果样本不属于任何组，标记为 NA
  ))

# 查看添加分组信息后的数据
print("添加分组信息后的数据：")
print(head(cibersort_results))
# 查看数据
print("CIBERSORTx 结果的前几行：")
print(head(cibersort_results))
# 提取免疫细胞比例
cell_proportions <- cibersort_results %>%
  select(-c("P-value", "Correlation", "RMSE"))

# 将数据转换为长格式以便 ggplot2 绘图
cell_proportions_long <- cell_proportions %>%
  pivot_longer(cols = -Mixture, names_to = "CellType", values_to = "Proportion")

# 查看整理后的数据
print("整理后的数据：")
print(head(cell_proportions_long))
# 绘制堆叠柱状图
ggplot(cell_proportions_long, aes(x = Mixture, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "CIBERSORTx 免疫细胞比例",
       x = "样本",
       y = "细胞比例") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存图表
ggsave("cibersortx_stacked_barplot.png", width = 10, height = 6)
# 将数据转换为矩阵格式
cell_proportions_matrix <- cell_proportions %>%
  column_to_rownames("Mixture") %>%
  as.matrix()

# 绘制热图
pheatmap(cell_proportions_matrix,
         scale = "none",  # 不对数据进行标准化
         clustering_distance_rows = "euclidean",  # 样本聚类距离
         clustering_distance_cols = "euclidean",  # 细胞类型聚类距离
         clustering_method = "complete",  # 聚类方法
         color = colorRampPalette(c("white", "blue"))(100),  # 颜色渐变
         main = "CIBERSORTx 免疫细胞比例热图")

# 保存热图
png("cibersortx_heatmap.png", width = 800, height = 600)
pheatmap(cell_proportions_matrix, ...)
dev.off()

# 绘制小提琴图
ggplot(cell_proportions_long, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_violin(trim = FALSE, scale = "width") +  # 小提琴图
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # 添加箱线图
  labs(title = "CIBERSORTx 免疫细胞比例分布",
       x = "免疫细胞类型",
       y = "细胞比例") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存图表
ggsave("cibersortx_violin_plot.png", width = 12, height = 6)

# 分组后提取免疫细胞比例
cell_proportions <- cibersort_results %>%
  select(-c("P-value", "Correlation", "RMSE"))

# 将数据转换为长格式以便 ggplot2 绘图
cell_proportions_long <- cell_proportions %>%
  pivot_longer(cols = -c(Mixture, Group), names_to = "CellType", values_to = "Proportion")

# 查看整理后的数据
print("整理后的数据：")
print(head(cell_proportions_long))
# 绘制分组小提琴图
ggplot(cell_proportions_long, aes(x = CellType, y = Proportion, fill = Group)) +
  geom_violin(trim = FALSE, scale = "width", position = position_dodge(0.8)) +  # 分组小提琴图
  geom_boxplot(width = 0.1, position = position_dodge(0.8), outlier.shape = NA) +  # 添加分组箱线图
  labs(title = "CIBERSORTx 免疫细胞比例分布（按组别）",
       x = "免疫细胞类型",
       y = "细胞比例",
       fill = "组别") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存图表
ggsave("cibersortx_grouped_violin_plot.png", width = 12, height = 6)


# 提取 T 细胞相关的列
t_cell_columns <- c("T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting", 
                    "T cells CD4 memory activated", "T cells follicular helper", 
                    "T cells regulatory (Tregs)", "T cells gamma delta")

t_cell_data <- cell_proportions_long %>%
  filter(CellType %in% t_cell_columns)

# 绘制 T 细胞的小提琴图
ggplot(t_cell_data, aes(x = CellType, y = Proportion, fill = Group)) +
  geom_violin(trim = FALSE, scale = "width", position = position_dodge(0.8)) +
  geom_boxplot(width = 0.1, position = position_dodge(0.8), outlier.shape = NA) +
  labs(title = "CIBERSORTx T 细胞比例分布（按组别）",
       x = "T 细胞类型",
       y = "细胞比例",
       fill = "组别") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存图表
ggsave("cibersortx_t_cell_violin_plot.png", width = 10, height = 6)

