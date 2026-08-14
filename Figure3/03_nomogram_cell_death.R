# ============================================================================
# Figure 3 细胞死亡版: Nomogram + Calibration + DCA
# ============================================================================
suppressMessages({
  library(pROC)
  library(rms)
})

setwd("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 3 机器学习降维锁定核心基因")
outdir <- "Results_CellDeath"

# ---- 1. 数据 ----
train_expr <- readRDS("Data/GSE43292_train_expr_entrez.rds")
train_group <- readRDS("Data/GSE43292_train_group.rds")
val_expr <- readRDS("Data/GSE28829_val_expr_entrez.rds")
val_group <- readRDS("Data/GSE28829_val_group.rds")
core_entrez <- readRDS(file.path(outdir, "core_genes.rds"))
use_genes <- intersect(intersect(rownames(train_expr), rownames(val_expr)), core_entrez)

X_train <- t(train_expr[use_genes, ])
X_val <- t(val_expr[use_genes, ])
y_train <- as.numeric(train_group)
y_val <- as.numeric(val_group)

# ---- 2. Nomogram ----
X_train_df <- as.data.frame(X_train)
X_val_df <- as.data.frame(X_val)

dd <- datadist(X_train_df)
options(datadist = "dd")
logit_model <- lrm(y_train ~ ., data = X_train_df, x = TRUE, y = TRUE)
# 训练集预测
train_pred <- predict(logit_model, newdata = X_train_df, type = "fitted")
val_pred <- predict(logit_model, newdata = X_val_df, type = "fitted")

# Nomogram 图 (用全部 17 基因) — PNG + SVG
# 修复：减少概率刻度(3个) + 加宽画布 + 加大下边距，避免 Risk 刻度数字重叠
png(file.path(outdir, "05_Nomogram.png"), width = 3600, height = 2600, res = 300)
par(mar = c(5, 2, 3, 2))
nom <- nomogram(logit_model, fun = plogis, fun.at = c(0.05, 0.5, 0.95),
                funlabel = "Risk of AS", lp = FALSE)
plot(nom, main = "Nomogram for Atherosclerosis Risk (Cell Death Signature)")
dev.off()

svg(file.path(outdir, "05_Nomogram.svg"), width = 16, height = 10)
par(mar = c(5, 2, 3, 2))
plot(nom, main = "Nomogram for Atherosclerosis Risk (Cell Death Signature)")
dev.off()

# ---- 3. Calibration — PNG + SVG ----
cal_train <- calibrate(logit_model, B = 200)
png(file.path(outdir, "06_Calibration.png"), width = 2000, height = 1800, res = 300)
plot(cal_train, xlab = "Predicted Probability", ylab = "Observed Probability",
     main = "Calibration Curve (Training)")
dev.off()

svg(file.path(outdir, "06_Calibration.svg"), width = 8, height = 7)
plot(cal_train, xlab = "Predicted Probability", ylab = "Observed Probability",
     main = "Calibration Curve (Training)")
dev.off()

# ---- 4. DCA (自实现) ----
dca_calc <- function(y, p, thresholds = seq(0.01, 0.99, by = 0.01)) {
  n <- length(y)
  nb_model <- numeric(length(thresholds))
  nb_all <- numeric(length(thresholds))
  for (i in seq_along(thresholds)) {
    pt <- thresholds[i]
    tp <- sum(p > pt & y == 1); fp <- sum(p > pt & y == 0)
    nb_model[i] <- tp/n - fp/n * (pt/(1-pt))
    nb_all[i] <- sum(y==1)/n - sum(y==0)/n * (pt/(1-pt))
  }
  data.frame(threshold = thresholds, model = nb_model, all = nb_all, none = 0)
}

dca_train <- dca_calc(y_train, train_pred)
dca_val <- dca_calc(y_val, val_pred)
write.csv(dca_train, file.path(outdir, "DCA_training.csv"), row.names = FALSE)
write.csv(dca_val, file.path(outdir, "DCA_validation.csv"), row.names = FALSE)

plot_dca <- function(dca_df, title, fname) {
  png(file.path(outdir, fname), width = 2000, height = 1700, res = 300)
  plot(dca_df$threshold, dca_df$model, type = "l", lwd = 2.5, col = "#C00000",
       xlab = "Threshold Probability", ylab = "Net Benefit", main = title,
       ylim = c(min(dca_df$model)-0.05, max(c(dca_df$all, dca_df$model))+0.05))
  lines(dca_df$threshold, dca_df$all, lwd = 2, col = "#4D4D4D", lty = 2)
  lines(dca_df$threshold, dca_df$none, lwd = 2, col = "#4D4D4D", lty = 3)
  legend("topright", legend = c("Model","Treat All","Treat None"),
         col = c("#C00000","#4D4D4D","#4D4D4D"), lwd = c(2.5,2,2), lty = c(1,2,3), bty = "n")
  dev.off()
  # SVG 版
  svg(file.path(outdir, sub("\\.png$", ".svg", fname)), width = 7, height = 6)
  plot(dca_df$threshold, dca_df$model, type = "l", lwd = 2.5, col = "#C00000",
       xlab = "Threshold Probability", ylab = "Net Benefit", main = title,
       ylim = c(min(dca_df$model)-0.05, max(c(dca_df$all, dca_df$model))+0.05))
  lines(dca_df$threshold, dca_df$all, lwd = 2, col = "#4D4D4D", lty = 2)
  lines(dca_df$threshold, dca_df$none, lwd = 2, col = "#4D4D4D", lty = 3)
  legend("topright", legend = c("Model","Treat All","Treat None"),
         col = c("#C00000","#4D4D4D","#4D4D4D"), lwd = c(2.5,2,2), lty = c(1,2,3), bty = "n")
  dev.off()
}
plot_dca(dca_train, "DCA - Training (GSE43292)", "07_DCA_training.png")
plot_dca(dca_val, "DCA - Validation (GSE28829)", "08_DCA_validation.png")

cat("\n✅ Nomogram + Calibration + DCA 完成\n")
cat("输出:", outdir, "\n")
