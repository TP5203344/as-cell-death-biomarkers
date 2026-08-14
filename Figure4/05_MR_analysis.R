# ============================================================================
# 方案A Step④: 孟德尔随机化 (MR) 因果验证 - 17个细胞死亡核心基因
# 数据源: eQTLGen (暴露) -> CAD (结局)
# ⚠️ OpenGWAS API 自 2024-05-01 起强制 JWT token
#    需在 https://api.opengwas.io/ 注册并设置环境变量 OPENGWAS_JWT
# ============================================================================
suppressMessages({
  library(TwoSampleMR)
  library(dplyr)
})

outdir <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 4 免疫浸润与诊断模型/MR"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- 17 个细胞死亡核心基因 (经验证训练/验证 AUC=0.966) ----
# 焦亡 5, 泛凋亡 4, NETosis 5, 自噬 1, 双硫死亡 1, 失巢凋亡 3, 铁死亡 1
core_genes <- c("AIM2","CASP1","NLRP3","RIPK3","PYCARD","IL18",
                "ITGAM","MPO","NCF2","CTSG","CTSD","PGD",
                "ZEB1","PRKAA2","BMF","ITGA5","CYBA")

# ENSG 映射 (org.Hs.eg.db 程序化验证 2026-08-11)
ensg_map <- c(
  AIM2    = "ENSG00000163568",
  CASP1   = "ENSG00000137752",
  NLRP3   = "ENSG00000162711",
  RIPK3   = "ENSG00000129465",
  PYCARD  = "ENSG00000103490",
  IL18    = "ENSG00000150782",
  ITGAM   = "ENSG00000169896",
  MPO     = "ENSG00000005381",
  NCF2    = "ENSG00000116701",
  CTSG    = "ENSG00000100448",
  CTSD    = "ENSG00000117984",
  PGD     = "ENSG00000142657",
  ZEB1    = "ENSG00000148516",
  PRKAA2  = "ENSG00000162409",
  BMF     = "ENSG00000104081",
  ITGA5   = "ENSG00000161638",
  CYBA    = "ENSG00000051523"
)

# 结局: CAD (冠心病) - CARDIoGRAMplusC4D
outcome_id <- "ebi-a-GCST005194"   # Nikpay 2015 CAD GWAS

# 检查 token
token <- Sys.getenv("OPENGWAS_JWT")
if (token == "") {
  cat("⚠️ 未设置 OPENGWAS_JWT 环境变量\n")
  cat("OpenGWAS API 自 2024-05-01 起强制认证\n")
  cat("请在 https://api.opengwas.io/ 注册并获取 token\n")
  cat("设置方法: export OPENGWAS_JWT='your_token_here'\n\n")
  cat("将尝试匿名访问 (可能失败)...\n\n")
}

results_all <- list()

for (gene in core_genes) {
  cat("\n========== ", gene, " ==========\n")
  ensg <- ensg_map[[gene]]
  exp_id <- paste0("eqtl-a-", ensg)
  cat("eQTLGen ID:", exp_id, "\n")
  
  exp_dat <- tryCatch({
    extract_instruments(outcomes = exp_id, p1 = 5e-8, clump = TRUE)
  }, error = function(e) { 
    cat("❌ eQTL 提取错误:", conditionMessage(e), "\n")
    NULL 
  })
  
  if (is.null(exp_dat) || nrow(exp_dat) == 0) {
    cat("无显著工具变量 (p<5e-8), 尝试 p<5e-5...\n")
    exp_dat <- tryCatch({
      extract_instruments(outcomes = exp_id, p1 = 5e-5, clump = TRUE)
    }, error = function(e) { cat("❌ 失败:", conditionMessage(e), "\n"); NULL })
  }
  
  if (is.null(exp_dat) || nrow(exp_dat) == 0) {
    cat("⚠️ 无工具变量, 跳过该基因\n")
    next
  }
  cat("✓ 工具变量数:", nrow(exp_dat), "\n")
  
  out_dat <- tryCatch({
    extract_outcome_data(snps = exp_dat$SNP, outcomes = outcome_id)
  }, error = function(e) { 
    cat("❌ 结局提取错误:", conditionMessage(e), "\n")
    NULL 
  })
  
  if (is.null(out_dat) || nrow(out_dat) == 0) {
    cat("⚠️ 结局无匹配 SNP, 跳过\n")
    next
  }
  
  dat <- tryCatch({
    harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)
  }, error = function(e) { 
    cat("❌ harmonise 错误:", conditionMessage(e), "\n")
    NULL 
  })
  
  if (is.null(dat)) next
  dat <- dat[dat$mr_keep, ]
  if (nrow(dat) < 3) {
    cat("⚠️ harmonise 后 SNP<3, 跳过\n")
    next
  }
  cat("✓ 分析 SNP 数:", nrow(dat), "\n")
  
  res <- tryCatch({
    mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
  }, error = function(e) { 
    cat("❌ MR 错误:", conditionMessage(e), "\n")
    NULL 
  })
  
  if (is.null(res)) next
  res$gene <- gene
  results_all[[gene]] <- res
  
  write.csv(res, file.path(outdir, paste0("mr_result_", gene, ".csv")), row.names = FALSE)
  saveRDS(dat, file.path(outdir, paste0("mr_dat_", gene, ".rds")))
  
  # 敏感性分析
  tryCatch({
    het <- mr_heterogeneity(dat)
    write.csv(het, file.path(outdir, paste0("mr_heterogeneity_", gene, ".csv")), row.names = FALSE)
  }, error = function(e) {})
  
  tryCatch({
    pleio <- mr_pleiotropy_test(dat)
    write.csv(pleio, file.path(outdir, paste0("mr_pleiotropy_", gene, ".csv")), row.names = FALSE)
  }, error = function(e) {})
  
  # 森林图 - SVG (TwoSampleMR 默认 PDF, 改为 SVG 便于投稿)
  tryCatch({
    p <- mr_forest_plot(mr_singlesnp(dat))
    svg(file.path(outdir, paste0("mr_forest_", gene, ".svg")), width = 7, height = 6)
    print(p[[length(p)]])
    dev.off()
  }, error = function(e) {})
  
  # 散点图 - SVG
  tryCatch({
    p2 <- mr_scatter_plot(res, dat)
    svg(file.path(outdir, paste0("mr_scatter_", gene, ".svg")), width = 6, height = 5)
    print(p2[[length(p2)]])
    dev.off()
  }, error = function(e) {})
}

if (length(results_all) > 0) {
  all_res <- do.call(rbind, results_all)
  write.csv(all_res, file.path(outdir, "MR_summary_all_genes.csv"), row.names = FALSE)
  cat("\n===== MR 汇总 (IVW) =====\n")
  ivw <- all_res[all_res$method == "Inverse variance weighted", ]
  if (nrow(ivw) > 0) {
    print(ivw[, c("gene", "nsnp", "b", "se", "pval")])
    cat("\n显著 (P<0.05):\n")
    sig <- ivw[ivw$pval < 0.05, c("gene", "nsnp", "b", "se", "pval")]
    if (nrow(sig) > 0) print(sig) else cat("无\n")
  }
} else {
  cat("\n❌ 未获得任何 MR 结果 (可能因 OpenGWAS token 缺失)\n")
  cat("请按上述提示设置 OPENGWAS_JWT 后重跑\n")
}
cat("\n完成。结果目录:", outdir, "\n")
