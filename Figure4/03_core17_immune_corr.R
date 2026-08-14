# ============================================================================
# Figure 4 补充: 17个细胞死亡核心基因 x 22免疫细胞 相关性热图
# 数据: GSE100927_immune_infiltration.csv (104样本 x 22细胞)
#        _gse100927_expr_full.csv (34964探针 x 104样本)
#        _probe_symbol_map.csv (50724探针->Symbol)
# ============================================================================
suppressMessages({
  library(pheatmap)
})

outdir <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 4 免疫浸润与诊断模型"
workdir <- "C:/Users/19267/.qclaw/workspace-agent-b69b3ebd/as_bioinfo"

# ---- 1. 17 个核心基因 ----
core_genes <- c("AIM2","CASP1","NLRP3","RIPK3","PYCARD","IL18",
                "ITGAM","MPO","NCF2","CTSG","CTSD","PGD",
                "ZEB1","PRKAA2","BMF","ITGA5","CYBA")

# ---- 2. 免疫浸润矩阵 ----
imm <- read.csv(file.path(outdir, "GSE100927_immune_infiltration.csv"), row.names = 1, check.names = FALSE)
imm <- as.matrix(imm[, 1:22])
cat("免疫浸润矩阵:", nrow(imm), "样本 x", ncol(imm), "细胞\n")

# ---- 3. 表达矩阵 + 探针映射 ----
cat("读取表达矩阵 (59MB, 约1分钟)...\n")
expr <- as.matrix(read.csv(file.path(workdir, "_gse100927_expr_full.csv"), row.names = 1, check.names = FALSE))
pmap <- read.csv(file.path(workdir, "_probe_symbol_map.csv"), stringsAsFactors = FALSE)
probe2sym <- pmap$symbol
names(probe2sym) <- pmap$probe
cat("表达矩阵:", nrow(expr), "探针 x", ncol(expr), "样本\n")

# ---- 4. 匹配 17 个核心基因 (多探针取均值) ----
gene_expr <- matrix(NA, nrow = length(core_genes), ncol = ncol(expr))
rownames(gene_expr) <- core_genes
colnames(gene_expr) <- colnames(expr)

for (g in core_genes) {
  probes <- names(probe2sym)[probe2sym == g]
  probes <- probes[probes %in% rownames(expr)]
  if (length(probes) == 0) next
  gene_expr[g, ] <- colMeans(expr[probes, , drop = FALSE])
}
hit <- rownames(gene_expr)[!is.na(gene_expr[, 1])]
cat("核心基因命中:", length(hit), "/", length(core_genes), ":", paste(hit, collapse = ","), "\n")
gene_expr <- gene_expr[!is.na(gene_expr[, 1]), , drop = FALSE]

# ---- 5. 样本对齐 ----
common <- intersect(colnames(gene_expr), rownames(imm))
cat("共同样本:", length(common), "\n")
gene_expr <- gene_expr[, common]
imm <- imm[common, ]

# ---- 6. Spearman 相关 ----
rho_mat <- matrix(NA, nrow = nrow(gene_expr), ncol = ncol(imm),
                  dimnames = list(rownames(gene_expr), colnames(imm)))
p_mat <- rho_mat
for (i in 1:nrow(gene_expr)) {
  for (j in 1:ncol(imm)) {
    ct <- suppressWarnings(cor.test(gene_expr[i, ], imm[, j], method = "spearman"))
    rho_mat[i, j] <- ct$estimate
    p_mat[i, j] <- ct$p.value
  }
}
star_mat <- ifelse(p_mat < 0.001, "***", ifelse(p_mat < 0.01, "**", ifelse(p_mat < 0.05, "*", "")))

# ---- 7. 保存 ----
write.csv(rho_mat, file.path(outdir, "core17_immune_correlation.csv"))
write.csv(p_mat, file.path(outdir, "core17_immune_correlation_pvalues.csv"))
saveRDS(rho_mat, file.path(outdir, "core17_immune_correlation.rds"))

# ---- 8. 热图 (PNG + SVG) ----
png(file.path(outdir, "core17_immune_corr_heatmap.png"), width = 3000, height = 2000, res = 300)
pheatmap(rho_mat,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(-1, 1, length.out = 101),
         display_numbers = star_mat,
         number_color = "grey30",
         fontsize_number = 6,
         cluster_cols = TRUE, cluster_rows = TRUE,
         main = "17 Cell Death Core Genes vs 22 Immune Cells (Spearman)",
         fontsize_row = 9, fontsize_col = 8)
dev.off()

svg(file.path(outdir, "core17_immune_corr_heatmap.svg"), width = 12, height = 8)
pheatmap(rho_mat,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(-1, 1, length.out = 101),
         display_numbers = star_mat,
         number_color = "grey30",
         fontsize_number = 6,
         cluster_cols = TRUE, cluster_rows = TRUE,
         main = "17 Cell Death Core Genes vs 22 Immune Cells (Spearman)",
         fontsize_row = 9, fontsize_col = 8)
dev.off()

cat("\n✅ 17 基因 x 免疫细胞相关性热图完成\n")
cat("显著相关 (P<0.05) 数量:", sum(p_mat < 0.05), "/", length(p_mat), "\n")
cat("输出:", file.path(outdir, "core17_immune_corr_heatmap.png"), "\n")
