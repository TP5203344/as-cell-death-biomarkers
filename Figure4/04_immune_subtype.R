# ============================================================================
# 方案A Step⑤: 免疫分型 (ConsensusClusterPlus)
# 基于 GSE100927 免疫浸润比例 (23 种免疫细胞) 对样本聚类
# 输出: 2-4 个免疫亚型 + 亚型间免疫细胞差异 + 与 AS/Control 关联
# ============================================================================
suppressMessages({
  library(ConsensusClusterPlus)
  library(pheatmap)
})

outdir <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 4 免疫浸润与诊断模型/ImmuneSubtype"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- 加载免疫浸润结果 ----
imm <- read.csv("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 4 免疫浸润与诊断模型/GSE100927_immune_infiltration.csv", row.names = 1, check.names = FALSE)
group_df <- read.csv("C:/Users/19267/.qclaw/workspace-agent-b69b3ebd/as_bioinfo/_gse100927_group.csv", check.names = FALSE)
colnames(group_df) <- c("sample", "title", "group")
group_map <- setNames(group_df$group, group_df$sample)

# 去掉 P-value 列
imm_mat <- imm[, !colnames(imm) %in% "P-value"]
cat("免疫浸润矩阵:", dim(imm_mat), "\n")

# 样本顺序与 group 匹配
samples <- rownames(imm_mat)
grp <- group_map[samples]
cat("AS:", sum(grp == "Atherosclerotic"), "Control:", sum(grp == "Control"), "\n")

# ---- ConsensusClusterPlus ----
set.seed(123)
# 转置: 行=样本, 列=细胞
mat <- t(imm_mat)
# 标准化
mat <- scale(mat)
mat[is.na(mat)] <- 0

# 聚类 2-4 类
res <- ConsensusClusterPlus(
  mat, maxK = 4, reps = 500, pItem = 0.8, pFeature = 1,
  clusterAlg = "km", distance = "euclidean", seed = 123,
  plot = "png", title = outdir, verbose = FALSE
)

# 保存聚类结果
clusters <- data.frame(sample = samples, group = grp)
for (k in 2:4) {
  cl <- res[[k]]$consensusClass
  clusters[[paste0("K", k)]] <- cl[clusters$sample]
}
write.csv(clusters, file.path(outdir, "immune_subtype_clusters.csv"), row.names = FALSE)

# ---- K=3 详细分析 ----
res3 <- res[[3]]
if (!is.null(res3)) {
  k <- 3
  cl3 <- res3$consensusClass
  cat("\n=== K=3 免疫分型 ===\n")
  print(table(cl3))
  
  # 亚型 vs 疾病状态
  ct <- table(cl3[samples], grp)
  print(ct)
  write.csv(ct, file.path(outdir, "subtype_vs_disease_K3.csv"))
  
  # 各亚型免疫细胞均值热图
  imm_annot <- imm_mat
  imm_annot$Subtype <- cl3[samples]
  subtype_means <- aggregate(. ~ Subtype, data = imm_annot, mean)
  rownames(subtype_means) <- paste0("Subtype", subtype_means$Subtype)
  subtype_means$Subtype <- NULL
  
  png(file.path(outdir, "subtype_immune_heatmap_K3.png"), width = 2000, height = 1600, res = 300)
  pheatmap(as.matrix(subtype_means), scale = "row",
           main = "Immune Cell Abundance by Subtype (K=3)",
           fontsize = 9, cluster_cols = TRUE, cluster_rows = TRUE,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
  dev.off()
  
  svg(file.path(outdir, "subtype_immune_heatmap_K3.svg"), width = 8, height = 6.5)
  pheatmap(as.matrix(subtype_means), scale = "row",
           main = "Immune Cell Abundance by Subtype (K=3)",
           fontsize = 9, cluster_cols = TRUE, cluster_rows = TRUE,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
  dev.off()
  
  # 亚型间差异检验 (Kruskal-Wallis)
  kw <- data.frame()
  for (cell in colnames(imm_mat)) {
    vals <- split(imm_mat[[cell]], cl3[samples])
    if (length(unique(cl3)) >= 2 && all(sapply(vals, length) >= 2)) {
      p <- kruskal.test(vals)$p.value
      kw <- rbind(kw, data.frame(Cell = cell, Kruskal_P = p))
    }
  }
  kw$FDR <- p.adjust(kw$Kruskal_P, method = "BH")
  write.csv(kw, file.path(outdir, "subtype_immune_diff_K3.csv"), row.names = FALSE)
  cat("\n亚型间差异免疫细胞 (FDR<0.05):", sum(kw$FDR < 0.05), "/", nrow(kw), "\n")
}

cat("\n✅ 免疫分型完成, 输出:", outdir, "\n")
