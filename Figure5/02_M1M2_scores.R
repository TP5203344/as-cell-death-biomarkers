# ============================================================================
# Figure 5 补丁: M1/M2 巨噬极化评分 (无需重跑 UMAP)
# ============================================================================
suppressMessages({
  library(Seurat)
  library(ggplot2)
})

setwd("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 5 单细胞定位")
dir.create("Results", showWarnings = FALSE)

# 直接读表达矩阵, 快速 Normalize 后算 M1/M2 分数
cat("读取表达矩阵...\n")
expr <- read.csv("C:/Users/19267/.qclaw/workspace-agent-b69b3ebd/as_bioinfo/_gse115469.csv.gz",
                 row.names = 1, check.names = FALSE)
so <- CreateSeuratObject(counts = expr, project = "GSE115469", min.cells = 3, min.features = 200)
so <- NormalizeData(so)

# 直接对全部细胞算 M1/M2 score (不需要聚类)
m1_genes <- intersect(c("NOS2", "IL1B", "TNF", "CXCL9", "CXCL10", "CXCL11", "CD80", "CD86"), rownames(so))
m2_genes <- intersect(c("MRC1", "CD163", "MSR1", "IL10", "TGFB1", "ARG1", "CCL18", "CCL22"), rownames(so))
cat("M1 marker 命中:", paste(m1_genes, collapse=","), "\n")
cat("M2 marker 命中:", paste(m2_genes, collapse=","), "\n")

data_mat <- GetAssayData(so, layer = "data")
so$M1_score <- colMeans(data_mat[m1_genes, , drop = FALSE])
so$M2_score <- colMeans(data_mat[m2_genes, , drop = FALSE])

# 保存 M1/M2 分数表
scores <- data.frame(cell = colnames(so), M1_score = so$M1_score, M2_score = so$M2_score)
write.csv(scores, "Results/M1_M2_scores.csv", row.names = FALSE)

# M1 vs M2 相关性与分布图 (不用 UMAP)
png("Results/10_M1_M2_distribution.png", width = 2400, height = 1200, res = 300)
par(mfrow = c(1, 2))
hist(so$M1_score, main = "M1 Score Distribution", xlab = "M1 score", col = "steelblue")
hist(so$M2_score, main = "M2 Score Distribution", xlab = "M2 score", col = "firebrick")
dev.off()

# M1/M2 ratio
so$M1M2_ratio <- so$M1_score / (so$M2_score + 0.01)
cat("\nM1/M2 ratio 统计:\n")
print(summary(so$M1M2_ratio))
cat("M1>M2 细胞占比:", round(mean(so$M1_score > so$M2_score) * 100, 1), "%\n")

saveRDS(data.frame(M1=so$M1_score, M2=so$M2_score, ratio=so$M1M2_ratio),
        "Results/M1_M2_scores.rds")
cat("\n✅ M1/M2 极化评分完成\n")
