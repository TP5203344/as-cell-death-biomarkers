# ============================================================================
# Figure 3 机器学习降维锁定核心基因 (方案A: 细胞死亡基因专项)
# step1_machine_learning_cell_death.R
# 输入: GSE43292_train_expr_entrez.rds (759基因 x 64样本, Entrez ID)
#       death_candidate_genes.csv (40个细胞死亡基因, SYMBOL)
# 流程: SYMBOL -> Entrez -> 过滤表达矩阵 -> LASSO + RF-RFE
# ============================================================================
suppressMessages({
  library(glmnet)
  library(randomForest)
  library(caret)
  library(org.Hs.eg.db)
})

setwd("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 3 机器学习降维锁定核心基因")
dir.create("Results_CellDeath", showWarnings = FALSE)

# ---- 1. SYMBOL -> Entrez 映射 ----
death_genes_sym <- read.csv("Data/death_candidate_genes.csv", stringsAsFactors = FALSE)$gene
cat("细胞死亡候选基因:", length(death_genes_sym), "个\n")
print(death_genes_sym)

mapped <- select(org.Hs.eg.db, keys = death_genes_sym, columns = "ENTREZID", keytype = "SYMBOL")
mapped <- mapped[!is.na(mapped$ENTREZID), ]
# 去重: 取第一个 Entrez
mapped <- mapped[!duplicated(mapped$SYMBOL), ]
cat("成功映射:", nrow(mapped), "/", length(death_genes_sym), "个\n")
print(mapped)
sym2entrez <- setNames(mapped$ENTREZID, mapped$SYMBOL)

# ---- 2. 加载训练数据 ----
train_expr <- readRDS("Data/GSE43292_train_expr_entrez.rds")
train_group <- readRDS("Data/GSE43292_train_group.rds")
cat("\n原始训练集:", nrow(train_expr), "基因 x", ncol(train_expr), "样本\n")

# ---- 3. 过滤到细胞死亡基因 ----
death_entrez <- as.character(mapped$ENTREZID)
common_entrez <- intersect(rownames(train_expr), death_entrez)
cat("交集基因:", length(common_entrez), "个\n")
# 保留顺序
common_entrez <- common_entrez[order(match(common_entrez, death_entrez))]

# 若交集太少, 放宽到所有 759 探针中含 death gene symbol 的
if (length(common_entrez) < 10) {
  cat("交集不足, 使用全部训练数据\n")
  X_all <- t(train_expr)
  y_all <- factor(train_group, levels = c(0,1), labels = c("Control","Disease"))
} else {
  train_death <- train_expr[common_entrez, ]
  cat("过滤后:", nrow(train_death), "基因 x", ncol(train_death), "样本\n")
  X_all <- t(train_death)
  y_all <- factor(train_group, levels = c(0,1), labels = c("Control","Disease"))
}
cat("最终特征数:", ncol(X_all), "\n")

# 标准化
X_scaled <- scale(X_all)
colnames(X_scaled) <- colnames(X_all)

# ---- 4. LASSO ----
set.seed(123)
cv_fit <- cv.glmnet(X_scaled, y_all, family = "binomial", alpha = 1, nfolds = 10, type.measure = "class")
lasso_coef_min <- as.matrix(coef(cv_fit, s = "lambda.min"))
lasso_coef_1se <- as.matrix(coef(cv_fit, s = "lambda.1se"))
lasso_genes_min <- rownames(lasso_coef_min)[lasso_coef_min[,1] != 0]
lasso_genes_min <- lasso_genes_min[lasso_genes_min != "(Intercept)"]
lasso_genes_1se <- rownames(lasso_coef_1se)[lasso_coef_1se[,1] != 0]
lasso_genes_1se <- lasso_genes_1se[lasso_genes_1se != "(Intercept)"]
cat("\nLASSO (lambda.min):", length(lasso_genes_min), "个基因\n")
cat(lasso_genes_min, "\n")
cat("LASSO (lambda.1se):", length(lasso_genes_1se), "个基因\n")
cat(lasso_genes_1se, "\n")

# 保存系数
lasso_df <- data.frame(Entrez = rownames(lasso_coef_min)[-1],
                       Coef = as.numeric(lasso_coef_min[-1,1]))
lasso_df <- lasso_df[order(-abs(lasso_df$Coef)), ]
write.csv(lasso_df, "Results_CellDeath/LASSO_coefficients.csv", row.names = FALSE)

# CV 曲线 — PNG + SVG
png("Results_CellDeath/01_LASSO_CV.png", width = 2400, height = 1800, res = 300)
par(mar = c(5,5,4,2))
plot(cv_fit, main = "LASSO Cross-Validation (Cell Death Genes, 10-fold)")
abline(v = log(cv_fit$lambda.min), col = "red", lty = 2)
abline(v = log(cv_fit$lambda.1se), col = "blue", lty = 2)
legend("topright", legend = c("lambda.min","lambda.1se"), col = c("red","blue"), lty = 2, bty = "n")
dev.off()

svg("Results_CellDeath/01_LASSO_CV.svg", width = 8, height = 6)
par(mar = c(5,5,4,2))
plot(cv_fit, main = "LASSO Cross-Validation (Cell Death Genes, 10-fold)")
abline(v = log(cv_fit$lambda.min), col = "red", lty = 2)
abline(v = log(cv_fit$lambda.1se), col = "blue", lty = 2)
legend("topright", legend = c("lambda.min","lambda.1se"), col = c("red","blue"), lty = 2, bty = "n")
dev.off()

# 系数路径图 — PNG + SVG
png("Results_CellDeath/02_LASSO_Coefficient_Path.png", width = 2400, height = 1800, res = 300)
par(mar = c(5,5,4,2))
plot(cv_fit$glmnet.fit, xvar = "lambda", label = FALSE, main = "LASSO Coefficient Path")
abline(v = log(cv_fit$lambda.min), col = "red", lty = 2)
abline(v = log(cv_fit$lambda.1se), col = "blue", lty = 2)
dev.off()

svg("Results_CellDeath/02_LASSO_Coefficient_Path.svg", width = 8, height = 6)
par(mar = c(5,5,4,2))
plot(cv_fit$glmnet.fit, xvar = "lambda", label = FALSE, main = "LASSO Coefficient Path")
abline(v = log(cv_fit$lambda.min), col = "red", lty = 2)
abline(v = log(cv_fit$lambda.1se), col = "blue", lty = 2)
dev.off()

# ---- 5. RF-RFE ----
if (length(lasso_genes_min) >= 5) {
  X_lasso <- X_scaled[, lasso_genes_min, drop = FALSE]
  set.seed(456)
  ctrl <- rfeControl(functions = rfFuncs, method = "cv", number = 5, verbose = FALSE)
  sizes <- seq(3, min(20, length(lasso_genes_min)), by = 1)
  rfe_result <- rfe(X_lasso, y_all, sizes = sizes, rfeControl = ctrl)
  cat("\nRFE 最优特征数:", length(rfe_result$optVariables), "\n")
  cat(rfe_result$optVariables, "\n")
  saveRDS(rfe_result, "Results_CellDeath/rfe_result.rds")
  
  png("Results_CellDeath/03_RFE_Accuracy.png", width = 2400, height = 1800, res = 300)
  plot(rfe_result, type = c("o","g"), main = "RF-RFE Feature Selection")
  dev.off()
  
  svg("Results_CellDeath/03_RFE_Accuracy.svg", width = 8, height = 6)
  plot(rfe_result, type = c("o","g"), main = "RF-RFE Feature Selection")
  dev.off()
} else {
  cat("\nLASSO 特征不足5个，跳过 RFE\n")
  rfe_result <- NULL
}

# ---- 6. 逻辑回归建模 + 预测 ----
if (length(lasso_genes_min) > 0) {
  X_final <- X_scaled[, lasso_genes_min, drop = FALSE]
  logit_model <- glm(y_all ~ X_final, family = binomial(link = "logit"))
  train_pred <- predict(logit_model, type = "response")
  names(train_pred) <- rownames(X_all)
  saveRDS(train_pred, "Results_CellDeath/train_pred.rds")
  
  # AUC
  suppressMessages(library(pROC))
  roc_obj <- roc(as.numeric(y_all) - 1, train_pred, quiet = TRUE)
  cat("\n训练集 AUC:", round(auc(roc_obj), 4), "\n")
  saveRDS(roc_obj, "Results_CellDeath/roc_train.rds")
  
  png("Results_CellDeath/04_LASSO_ROC.png", width = 1800, height = 1600, res = 300)
  plot(roc_obj, col = "#C00000", lwd = 2.5, main = paste0("LASSO ROC (AUC=", round(auc(roc_obj),3), ")"))
  dev.off()
  
  svg("Results_CellDeath/04_LASSO_ROC.svg", width = 6, height = 5.5)
  plot(roc_obj, col = "#C00000", lwd = 2.5, main = paste0("LASSO ROC (AUC=", round(auc(roc_obj),3), ")"))
  dev.off()
  
  # 核心基因列表
  core_genes_entrez <- lasso_genes_min
  # SYMBOL 映射
  entrez2sym <- setNames(mapped$SYMBOL, mapped$ENTREZID)
  core_genes_sym <- entrez2sym[core_genes_entrez]
  core_df <- data.frame(Entrez = core_genes_entrez, Symbol = unname(core_genes_sym),
                        stringsAsFactors = FALSE)
  write.csv(core_df, "Results_CellDeath/core_genes.csv", row.names = FALSE)
  saveRDS(core_df$Entrez, "Results_CellDeath/core_genes.rds")
  
  cat("\n===== 核心基因 (LASSO 选中的细胞死亡基因) =====\n")
  print(core_df)
  cat("\n核心基因所属细胞死亡类型:\n")
  death_types <- c(Ferroptosis="CYBB,DPP4,FTL,HILPDA,HMOX1,PRKAA2,SLC40A1,STEAP3",
                   Pyroptosis="AIM2,CASP1,CASP8,IL18,IL1B,NLRC4,NLRP3,PYCARD,TNF",
                   PANoptosis="AIM2,CASP1,CASP8,NLRP3,RIPK3,ZBP1",
                   Necroptosis="CASP8,RIPK3,TNF",Autophagy="CTSB,CTSD",
                   NETosis="CTSG,CYBA,CYBB,ITGAM,ITGB2,MPO,NCF1,NCF2,NCF4,RAC2,TLR2")
  for (dt in names(death_types)) {
    hits <- intersect(core_df$Symbol, strsplit(death_types[[dt]], ",")[[1]])
    if (length(hits) > 0) cat(" -", dt, ":", paste(hits, collapse=", "), "\n")
  }
} else {
  cat("LASSO 未选中任何基因!\n")
}

cat("\n✅ 细胞死亡基因机器学习完成\n")
cat("结果目录: Results_CellDeath/\n")
