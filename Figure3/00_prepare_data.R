# ============================================================
# Figure 3 机器学习降维锁定核心基因
# step0_prepare_data.R — 跨平台数据准备（统一 Entrez ID）
# 候选基因池: GSE100927 DEGs ∩ GSE43292 WGCNA模块基因
# 训练集: GSE43292 (64样本) | 验证集: GSE28829 (29样本)
# 原则: 跨平台交集 + 单平台训练，唯一标识符 = Entrez ID
# ============================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(hgu133plus2.db)
})

setwd("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 3 机器学习降维锁定核心基因")
dir.create("Data", showWarnings = FALSE)
dir.create("Results", showWarnings = FALSE)

# ------------------------------------------------------------
# 1. GSE100927 DEGs（发现数据集）
# ------------------------------------------------------------
deg <- read.csv("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 1 研究设计与转录组差异分析/火山图GSE100927/GSE100927_DE_results.csv",
                stringsAsFactors = FALSE)
cat("GSE100927 探针总数:", nrow(deg), "\n")

# 差异基因筛选: |logFC| > 1 & FDR < 0.05
deg_sig <- deg[abs(deg$logFC) > 1 & deg$FDR < 0.05, ]
cat("显著 DEGs:", nrow(deg_sig), "\n")

# Symbol -> Entrez ID（唯一标识符）
sym <- unique(deg_sig$external_gene_name)
sym <- sym[!is.na(sym) & sym != ""]
map <- select(org.Hs.eg.db, keys = sym, columns = "ENTREZID", keytype = "SYMBOL")
map <- map[!is.na(map$ENTREZID), ]
# 一个 symbol 对应多个 Entrez 时取第一个
map <- map[!duplicated(map$SYMBOL), ]
deg_entrez <- unique(map$ENTREZID)
cat("DEGs 映射到 Entrez ID:", length(deg_entrez), "\n")

# ------------------------------------------------------------
# 2. GSE43292 WGCNA 模块基因（训练集）
# ------------------------------------------------------------
mod_dir <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 2 WGCNA鉴定斑块演进模块/06_Data_Output"
mod_files <- list.files(mod_dir, pattern = "^05_Module_.*_Genes\\.txt$", full.names = TRUE)
cat("\nWGCNA 模块文件:\n")
print(basename(mod_files))

module_genes <- list()
for (f in mod_files) {
  mod_name <- gsub("^05_Module_|_Genes\\.txt$", "", basename(f))
  module_genes[[mod_name]] <- readLines(f)
}
# 模块基因 -> Entrez
all_mod_sym <- unique(unlist(module_genes))
mod_map <- select(org.Hs.eg.db, keys = all_mod_sym, columns = "ENTREZID", keytype = "SYMBOL")
mod_map <- mod_map[!is.na(mod_map$ENTREZID), ]
mod_map <- mod_map[!duplicated(mod_map$SYMBOL), ]
mod_entrez <- unique(mod_map$ENTREZID)
cat("WGCNA 模块基因映射到 Entrez ID:", length(mod_entrez), "\n")

# ------------------------------------------------------------
# 3. 候选基因池 = DEGs ∩ WGCNA模块
# ------------------------------------------------------------
candidate_entrez <- intersect(deg_entrez, mod_entrez)
cat("\n候选基因池 (DEGs ∩ WGCNA):", length(candidate_entrez), "\n")

# 保存候选基因（带 symbol 便于阅读）
cand_map <- mod_map[mod_map$ENTREZID %in% candidate_entrez, ]
write.csv(cand_map, "Data/candidate_genes_entrez.csv", row.names = FALSE)

# ------------------------------------------------------------
# 4. GSE43292 训练集表达矩阵（行 = Entrez ID）
# ------------------------------------------------------------
train_mat <- read.csv("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 2 WGCNA鉴定斑块演进模块/Sample Type Matrix.csv",
                      row.names = 1, check.names = FALSE)
cat("\nGSE43292 原始矩阵:", nrow(train_mat), "基因 x", ncol(train_mat), "样本\n")

# 行名 symbol -> Entrez
train_sym <- rownames(train_mat)
train_map <- select(org.Hs.eg.db, keys = train_sym, columns = "ENTREZID", keytype = "SYMBOL")
train_map <- train_map[!is.na(train_map$ENTREZID), ]
train_map <- train_map[!duplicated(train_map$SYMBOL), ]
# 匹配候选基因
train_cand <- train_map[train_map$ENTREZID %in% candidate_entrez, ]
cat("训练集匹配候选基因:", nrow(train_cand), "\n")

# 构建 Entrez 行名矩阵（多探针取均值）
train_expr <- train_mat[train_cand$SYMBOL, , drop = FALSE]
rownames(train_expr) <- train_cand$ENTREZID
# 同一 Entrez 多行时取均值
if (any(duplicated(rownames(train_expr)))) {
  train_expr <- rowsum(train_expr, group = rownames(train_expr), reorder = FALSE) / 
    as.numeric(table(rownames(train_expr)))
}
cat("训练集最终矩阵:", nrow(train_expr), "基因 x", ncol(train_expr), "样本\n")

# 样本分组: _con = Control(0), _tre = Disease(1)
train_group <- ifelse(grepl("_con$", colnames(train_expr)), 0, 1)
names(train_group) <- colnames(train_expr)
cat("训练集: Control =", sum(train_group == 0), ", Disease =", sum(train_group == 1), "\n")

saveRDS(train_expr, "Data/GSE43292_train_expr_entrez.rds")
saveRDS(train_group, "Data/GSE43292_train_group.rds")
write.csv(data.frame(Sample = colnames(train_expr), Group = train_group),
          "Data/GSE43292_train_group.csv", row.names = FALSE)

# ------------------------------------------------------------
# 5. GSE28829 验证集表达矩阵（GPL570 探针 -> Entrez）
# ------------------------------------------------------------
gz <- gzfile("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/数据文件/GSE28829_series_matrix.txt.gz", "rt")
lines <- readLines(gz)
close(gz)
ds <- grep("!series_matrix_table_begin", lines)
hdr <- strsplit(lines[ds + 1], "\t")[[1]]
data_lines <- lines[(ds + 1):(grep("!series_matrix_table_end", lines) - 1)]
val_mat <- read.table(text = data_lines, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cat("\nGSE28829 矩阵:", nrow(val_mat), "探针 x", ncol(val_mat) - 1, "样本\n")

# 分组: phenotype advanced vs early
ph_line <- lines[41]
pheno <- strsplit(sub("^!Sample_characteristics_ch1\\s*", "", ph_line), "\t")[[1]]
pheno <- gsub('^"phenotype: |"$', "", pheno)
val_group <- ifelse(pheno == "advanced atherosclerotic plaque", 1, 0)
names(val_group) <- colnames(val_mat)[-1]
cat("验证集: early =", sum(val_group == 0), ", advanced =", sum(val_group == 1), "\n")

# 探针 -> Entrez
probe_ids <- as.character(val_mat$ID_REF)
probe_map <- select(hgu133plus2.db, keys = probe_ids, columns = "ENTREZID", keytype = "PROBEID")
probe_map <- probe_map[!is.na(probe_map$ENTREZID), ]
probe_map <- probe_map[!duplicated(probe_map$PROBEID), ]
cat("GSE28829 探针映射到 Entrez:", length(unique(probe_map$ENTREZID)), "\n")

# 匹配候选基因
val_cand <- probe_map[probe_map$ENTREZID %in% candidate_entrez, ]
cat("验证集匹配候选基因:", nrow(val_cand), "\n")

val_expr <- as.matrix(val_mat[match(val_cand$PROBEID, val_mat$ID_REF), -1])
rownames(val_expr) <- val_cand$ENTREZID
if (any(duplicated(rownames(val_expr)))) {
  val_expr <- rowsum(val_expr, group = rownames(val_expr), reorder = FALSE) / 
    as.numeric(table(rownames(val_expr)))
}
cat("验证集最终矩阵:", nrow(val_expr), "基因 x", ncol(val_expr), "样本\n")

saveRDS(val_expr, "Data/GSE28829_val_expr_entrez.rds")
saveRDS(val_group, "Data/GSE28829_val_group.rds")
write.csv(data.frame(Sample = names(val_group), Group = val_group),
          "Data/GSE28829_val_group.csv", row.names = FALSE)

# ------------------------------------------------------------
# 6. 交集总结
# ------------------------------------------------------------
cat("\n========== 候选基因池 ==========\n")
cat("DEGs Entrez:", length(deg_entrez), "\n")
cat("WGCNA模块 Entrez:", length(mod_entrez), "\n")
cat("交集候选基因:", length(candidate_entrez), "\n")
cat("候选基因列表 (Entrez | Symbol):\n")
print(cand_map[, c("ENTREZID", "SYMBOL")])
write.csv(cand_map[, c("ENTREZID", "SYMBOL")], "Results/candidate_genes_final.csv", row.names = FALSE)

cat("\nstep0 完成！\n")
