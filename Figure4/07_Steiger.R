# ============================================================================
# 方案A Step④扩展: Steiger 方向性过滤 (因果方向核查)
# 检验: 基因表达(暴露) -> CAD(结局) 方向是否成立, 排除反向因果(CAD->表达)
# 数据: 复用 05_MR_analysis.R 生成的 mr_dat_*.rds (TwoSampleMR harmonised)
# 输出: MR/Steiger/ 下 per-gene CSV + Steiger_summary_all_genes.csv
# ============================================================================
suppressMessages({
  library(TwoSampleMR)
})

mrdir <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 4 免疫浸润与诊断模型/MR"
outdir <- file.path(mrdir, "Steiger")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

genes <- c("AIM2","BMF","CASP1","CTSD","CTSG","CYBA","ITGA5","ITGAM","MPO","NLRP3")

results <- list()

for (gene in genes) {
  cat("\n===== ", gene, " =====\n", sep = "")
  f <- file.path(mrdir, paste0("mr_dat_", gene, ".rds"))
  if (!file.exists(f)) { cat("  rds 不存在, 跳过\n"); next }
  dat <- readRDS(f)
  dat <- dat[dat$mr_keep, ]
  if (nrow(dat) < 3) { cat("  SNP<3, 跳过\n"); next }

  # 每个 SNP 对暴露/结局的解释方差 R²
  dat$r.exposure <- get_r_from_bsen(b = dat$beta.exposure, se = dat$se.exposure, n = dat$samplesize.exposure)
  dat$r.outcome  <- get_r_from_bsen(b = dat$beta.outcome,  se = dat$se.outcome,  n = dat$samplesize.outcome)

  # 1) 方向性检验
  dir_res <- tryCatch(directionality_test(dat), error = function(e) {
    cat("  directionality_test 错误:", conditionMessage(e), "\n"); NULL
  })
  # 2) Steiger 过滤
  st_res <- tryCatch(steiger_filtering(dat), error = function(e) {
    cat("  steiger_filtering 错误:", conditionMessage(e), "\n"); NULL
  })

  if (is.null(dir_res) && is.null(st_res)) next

  row <- data.frame(gene = gene, nsnp_total = nrow(dat), stringsAsFactors = FALSE)
  if (!is.null(dir_res)) {
    row$correct_causal_direction <- dir_res$correct_causal_direction[1]
    row$steiger_pval <- dir_res$steiger_pval[1]
    row$snp_r2_exposure_mean <- mean(dat$r.exposure, na.rm = TRUE)
    row$snp_r2_outcome_mean  <- mean(dat$r.outcome,  na.rm = TRUE)
  }
  if (!is.null(st_res)) {
    row$nsnp_kept_steiger <- sum(st_res$steiger_dir, na.rm = TRUE)
    row$nsnp_removed_steiger <- sum(!st_res$steiger_dir, na.rm = TRUE)
  }
  results[[gene]] <- row

  if (!is.null(st_res)) {
    write.csv(st_res[, c("SNP","r.exposure","r.outcome","steiger_dir","steiger_pval")],
              file.path(outdir, paste0("steiger_", gene, ".csv")), row.names = FALSE)
  }
}

if (length(results) > 0) {
  all <- do.call(rbind, results)
  write.csv(all, file.path(outdir, "Steiger_summary_all_genes.csv"), row.names = FALSE)
  cat("\n===== Steiger 方向性汇总 =====\n")
  print(all)
}
cat("\n完成。结果目录:", outdir, "\n")
