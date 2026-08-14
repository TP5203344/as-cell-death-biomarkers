# ============================================================================
# Figure 3 机器学习降维锁定核心基因 (细胞死亡版)
# step2_validation_cell_death.R — GSE28829 外部验证 + ROC + Nomogram
# ============================================================================
suppressMessages({
  library(pROC)
  library(org.Hs.eg.db)
})

setwd("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 3 机器学习降维锁定核心基因")
outdir <- "Results_CellDeath"
dir.create(outdir, showWarnings = FALSE)

# ---- 1. 加载训练/验证数据 ----
train_expr <- readRDS("Data/GSE43292_train_expr_entrez.rds")
train_group <- readRDS("Data/GSE43292_train_group.rds")
val_expr <- readRDS("Data/GSE28829_val_expr_entrez.rds")
val_group <- readRDS("Data/GSE28829_val_group.rds")
cat("训练集:", nrow(train_expr), "x", ncol(train_expr), "\n")
cat("验证集:", nrow(val_expr), "x", ncol(val_expr), "\n")

# 核心基因 (LASSO 17 个)
core_entrez <- readRDS(file.path(outdir, "core_genes.rds"))
cat("核心基因:", length(core_entrez), "\n")
print(core_entrez)

# ---- 2. 提取核心基因表达 ----
train_common <- intersect(rownames(train_expr), core_entrez)
val_common <- intersect(rownames(val_expr), core_entrez)
cat("训练集命中:", length(train_common), "/", length(core_entrez), "\n")
cat("验证集命中:", length(val_common), "/", length(core_entrez), "\n")
missing_train <- setdiff(core_entrez, train_common)
missing_val <- setdiff(core_entrez, val_common)
if (length(missing_train) > 0) {
  cat("训练集缺失:", paste(missing_train, collapse=", "), "\n")
}
if (length(missing_val) > 0) {
  cat("验证集缺失:", paste(missing_val, collapse=", "), "\n")
}

# 用共有基因
use_genes <- intersect(train_common, val_common)
cat("两个数据集共有核心基因:", length(use_genes), "\n")
print(use_genes)

if (length(use_genes) < 3) {
  stop("共有基因不足 3 个, 无法建模!")
}

X_train <- t(train_expr[use_genes, ])
X_val <- t(val_expr[use_genes, ])
y_train <- factor(train_group, levels = c(0,1), labels = c("Control","Disease"))
y_val <- factor(val_group, levels = c(0,1), labels = c("Control","Disease"))

# ---- 3. 逻辑回归 (训练集拟合) ----
set.seed(123)
logit_model <- glm(y_train ~ ., data = data.frame(X_train), family = binomial)
train_pred <- predict(logit_model, newdata = data.frame(X_train), type = "response")
val_pred <- predict(logit_model, newdata = data.frame(X_val), type = "response")
saveRDS(train_pred, file.path(outdir, "train_pred.rds"))
saveRDS(val_pred, file.path(outdir, "val_pred.rds"))

roc_train <- roc(as.numeric(y_train) - 1, train_pred, quiet = TRUE)
roc_val <- roc(as.numeric(y_val) - 1, val_pred, quiet = TRUE)
cat("\n训练集 AUC:", round(auc(roc_train), 4), "\n")
cat("验证集 AUC:", round(auc(roc_val), 4), "\n")
saveRDS(roc_train, file.path(outdir, "roc_train.rds"))
saveRDS(roc_val, file.path(outdir, "roc_val.rds"))

# ROC 汇总表
roc_summary <- data.frame(
  Dataset = c("Training (GSE43292)", "Validation (GSE28829)"),
  AUC = c(auc(roc_train), auc(roc_val)),
  CI_lower = c(ci.auc(roc_train)[1], ci.auc(roc_val)[1]),
  CI_upper = c(ci.auc(roc_train)[3], ci.auc(roc_val)[3]),
  N = c(length(train_pred), length(val_pred))
)
write.csv(roc_summary, file.path(outdir, "ROC_summary.csv"), row.names = FALSE)
print(roc_summary)

# ---- 4. ROC 图 - PNG + SVG ----
png(file.path(outdir, "04_ROC_train_val.png"), width = 2000, height = 1800, res = 300)
plot(roc_train, col = "#C00000", lwd = 2.5, main = "ROC Curves - Cell Death Signature")
plot(roc_val, col = "#1F4E79", lwd = 2.5, add = TRUE)
legend("bottomright",
       legend = c(paste0("Training AUC=", round(auc(roc_train),3)),
                  paste0("Validation AUC=", round(auc(roc_val),3))),
       col = c("#C00000","#1F4E79"), lwd = 2.5, bty = "n")
dev.off()

svg(file.path(outdir, "04_ROC_train_val.svg"), width = 7, height = 6)
plot(roc_train, col = "#C00000", lwd = 2.5, main = "ROC Curves - Cell Death Signature")
plot(roc_val, col = "#1F4E79", lwd = 2.5, add = TRUE)
legend("bottomright",
       legend = c(paste0("Training AUC=", round(auc(roc_train),3)),
                  paste0("Validation AUC=", round(auc(roc_val),3))),
       col = c("#C00000","#1F4E79"), lwd = 2.5, bty = "n")
dev.off()

# ---- 5. Bootstrap 验证集 CI ----
set.seed(999)
boot_auc <- numeric(1000)
n_val <- length(val_pred)
for (i in 1:1000) {
  idx <- sample(1:n_val, n_val, replace = TRUE)
  boot_auc[i] <- auc(roc(as.numeric(y_val)[idx] - 1, val_pred[idx], quiet = TRUE))
}
cat("\n验证集 Bootstrap AUC 95% CI:", round(quantile(boot_auc, c(0.025, 0.975)), 4), "\n")
saveRDS(boot_auc, file.path(outdir, "val_boot_aucs.rds"))

# ---- 6. 基因 SYMBOL 映射 ----
mapped <- select(org.Hs.eg.db, keys = use_genes, columns = "SYMBOL", keytype = "ENTREZID")
mapped <- mapped[!duplicated(mapped$ENTREZID), ]
entrez2sym <- setNames(mapped$SYMBOL, mapped$ENTREZID)
core_sym <- entrez2sym[use_genes]

# 系数表
coef_df <- data.frame(
  Entrez = use_genes,
  Symbol = unname(core_sym),
  Coef = round(coef(logit_model)[-1], 4),
  OR = round(exp(coef(logit_model)[-1]), 4),
  P = summary(logit_model)$coefficients[-1, 4]
)
write.csv(coef_df, file.path(outdir, "logistic_coefficients.csv"), row.names = FALSE)
cat("\n===== 逻辑回归系数 =====\n")
print(coef_df)

cat("\n✅ 细胞死亡基因验证完成\n")
