# ============================================================================
# 方案A: DCA 决策曲线分析 (自实现, 因 rmda 已从 CRAN 下架)
# 基于 Figure 3 已有的 logistic 模型 + 训练/验证集
# ============================================================================
suppressMessages({
  library(pROC)
})

setwd("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 3 机器学习降维锁定核心基因")
outdir <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 4 免疫浸润与诊断模型"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- 自实现 DCA ----
# 输入: 预测概率 p, 真实结局 y (0/1)
# 输出: 不同阈值下的净获益 (Net Benefit)
dca_calc <- function(y, p, thresholds = seq(0.01, 0.99, by = 0.01)) {
  n <- length(y)
  nb_model <- numeric(length(thresholds))
  nb_all <- numeric(length(thresholds))
  nb_none <- numeric(length(thresholds))
  
  for (i in seq_along(thresholds)) {
    pt <- thresholds[i]
    # 模型策略: 预测概率 > pt 则干预
    tp <- sum(p > pt & y == 1)
    fp <- sum(p > pt & y == 0)
    nb_model[i] <- tp / n - fp / n * (pt / (1 - pt))
    # 全部干预
    tpa <- sum(y == 1)
    nb_all[i] <- tpa / n - (n - tpa) / n * (pt / (1 - pt))
    # 不干预
    nb_none[i] <- 0
  }
  data.frame(threshold = thresholds, model = nb_model, all = nb_all, none = nb_none)
}

# ---- 加载模型和数据 ----
train_pred <- readRDS("Results/train_pred.rds")
val_pred <- readRDS("Results/val_pred.rds")
train_group <- readRDS("../Figure 3 机器学习降维锁定核心基因/Data/GSE43292_train_group.rds")
val_group <- readRDS("../Figure 3 机器学习降维锁定核心基因/Data/GSE28829_val_group.rds")

cat("train_pred class:", class(train_pred), "len:", length(train_pred), "\n")
cat("val_pred class:", class(val_pred), "len:", length(val_pred), "\n")
cat("train_group:", length(train_group), "val_group:", length(val_group), "\n")

# 预测概率向量
if (is.list(train_pred)) {
  train_p <- train_pred$pred  # 可能是 pred 列
} else {
  train_p <- train_pred
}
if (is.list(val_pred)) {
  val_p <- val_pred$pred
} else {
  val_p <- val_pred
}

# 确保是概率 (0-1)
cat("train_p range:", range(train_p), "\n")
cat("val_p range:", range(val_p), "\n")

y_train <- as.numeric(train_group)
y_val <- as.numeric(val_group)

# ---- DCA 计算 ----
dca_train <- dca_calc(y_train, as.numeric(train_p))
dca_val <- dca_calc(y_val, as.numeric(val_p))

write.csv(dca_train, file.path(outdir, "DCA_training.csv"), row.names = FALSE)
write.csv(dca_val, file.path(outdir, "DCA_validation.csv"), row.names = FALSE)

# ---- DCA 绘图 ----
plot_dca <- function(dca_df, title, fname) {
  # PNG
  png(file.path(outdir, fname), width = 2400, height = 2000, res = 300)
  plot(dca_df$threshold, dca_df$model, type = "l", lwd = 2.5, col = "#C00000",
       xlab = "Threshold Probability", ylab = "Net Benefit",
       main = title, ylim = c(min(dca_df$model) - 0.05, max(c(dca_df$all, dca_df$model)) + 0.05))
  lines(dca_df$threshold, dca_df$all, lwd = 2, col = "#4D4D4D", lty = 2)
  lines(dca_df$threshold, dca_df$none, lwd = 2, col = "#4D4D4D", lty = 3)
  legend("topright", legend = c("Model", "Treat All", "Treat None"),
         col = c("#C00000", "#4D4D4D", "#4D4D4D"), lwd = c(2.5, 2, 2), lty = c(1, 2, 3), bty = "n",
         x.intersp = 0.8, y.intersp = 1.1, cex = 0.9, text.width = 0.35)
  dev.off()
  # SVG
  svg(sub("\\.png$", ".svg", file.path(outdir, fname)), width = 8, height = 6.5)
  plot(dca_df$threshold, dca_df$model, type = "l", lwd = 2.5, col = "#C00000",
       xlab = "Threshold Probability", ylab = "Net Benefit",
       main = title, ylim = c(min(dca_df$model) - 0.05, max(c(dca_df$all, dca_df$model)) + 0.05))
  lines(dca_df$threshold, dca_df$all, lwd = 2, col = "#4D4D4D", lty = 2)
  lines(dca_df$threshold, dca_df$none, lwd = 2, col = "#4D4D4D", lty = 3)
  legend("topright", legend = c("Model", "Treat All", "Treat None"),
         col = c("#C00000", "#4D4D4D", "#4D4D4D"), lwd = c(2.5, 2, 2), lty = c(1, 2, 3), bty = "n",
         x.intersp = 0.8, y.intersp = 1.1, cex = 0.9, text.width = 0.35)
  dev.off()
}

plot_dca(dca_train, "Decision Curve Analysis - Training Set (GSE43292)", "DCA_training.png")
plot_dca(dca_val, "Decision Curve Analysis - Validation Set (GSE28829)", "DCA_validation.png")

# ---- 合并 ROC 图 (训练+验证) ----
roc_train <- roc(y_train, as.numeric(train_p), quiet = TRUE)
roc_val <- roc(y_val, as.numeric(val_p), quiet = TRUE)

png(file.path(outdir, "ROC_train_val_combined.png"), width = 2400, height = 2000, res = 300)
plot(roc_train, col = "#C00000", lwd = 2.5, main = "ROC Curves - Training & Validation")
plot(roc_val, col = "#1F4E79", lwd = 2.5, add = TRUE)
legend("bottomright",
       legend = c(paste0("Training AUC = ", round(auc(roc_train), 4)),
                  paste0("Validation AUC = ", round(auc(roc_val), 4))),
       col = c("#C00000", "#1F4E79"), lwd = 2.5, bty = "n")
dev.off()

svg(file.path(outdir, "ROC_train_val_combined.svg"), width = 8, height = 6.5)
plot(roc_train, col = "#C00000", lwd = 2.5, main = "ROC Curves - Training & Validation")
plot(roc_val, col = "#1F4E79", lwd = 2.5, add = TRUE)
legend("bottomright",
       legend = c(paste0("Training AUC = ", round(auc(roc_train), 4)),
                  paste0("Validation AUC = ", round(auc(roc_val), 4))),
       col = c("#C00000", "#1F4E79"), lwd = 2.5, bty = "n")
dev.off()

cat("\n===== DCA 完成 =====\n")
cat("训练集 AUC:", round(auc(roc_train), 4), "\n")
cat("验证集 AUC:", round(auc(roc_val), 4), "\n")
cat("输出:", outdir, "\n")
