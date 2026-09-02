# ============================================================================
# 方案A Step④扩展: Bayesian Colocalization (coloc, PPH4)
# 检验: eQTL 与 CAD GWAS 信号是否共享同一因果变异
# 数据: 完整 cis 摘要 (±500kb, 见 08_download_cis.R)
# ============================================================================
suppressMessages({
  library(coloc)
  library(data.table)
})

cisdir <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 4 免疫浸润与诊断模型/MR/cis"
outdir <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 4 免疫浸润与诊断模型/MR/coloc"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

genes <- c("AIM2","BMF","CASP1","CTSD","CTSG","CYBA","ITGA5","ITGAM","MPO","NLRP3")
s_cad <- 34541 / (34541 + 261984)  # CAD 病例比例

results <- list()
for (gene in genes) {
  exp_f <- file.path(cisdir, paste0("cis_exposure_", gene, ".csv"))
  out_f <- file.path(cisdir, paste0("cis_outcome_", gene, ".csv"))
  if (!file.exists(exp_f) || !file.exists(out_f)) next

  exp <- fread(exp_f); out <- fread(out_f)
  common <- intersect(exp$rsid, out$rsid)
  if (length(common) < 50) next

  e <- exp[match(common, exp$rsid), ]
  o <- out[match(common, out$rsid), ]

  D1 <- list(beta = e$beta, varbeta = e$se^2, snp = e$rsid, position = e$position,
             MAF = ifelse(e$eaf > 0.5, 1 - e$eaf, e$eaf), N = e$n, type = "quant", sdY = 1)
  D2 <- list(beta = o$beta, varbeta = o$se^2, snp = o$rsid, position = o$position,
             MAF = ifelse(o$eaf > 0.5, 1 - o$eaf, o$eaf), N = o$n, type = "cc", s = s_cad)

  keep <- complete.cases(D1$beta, D1$varbeta, D1$MAF, D2$beta, D2$varbeta, D2$MAF)
  vec_cols <- c("beta","varbeta","snp","position","MAF","N")
  for (cc in vec_cols) {
    if (length(D1[[cc]]) == length(keep)) D1[[cc]] <- D1[[cc]][keep]
    if (length(D2[[cc]]) == length(keep)) D2[[cc]] <- D2[[cc]][keep]
  }

  res <- tryCatch(coloc.abf(dataset1 = D1, dataset2 = D2), error = function(e) NULL)
  if (is.null(res)) next

  s <- res$summary; snp_res <- res$results
  top_snp <- snp_res[which.max(snp_res$SNP.PP.H4), ]
  results[[gene]] <- data.frame(
    gene = gene, nsnp_common = length(common),
    PP.H0 = s["PP.H0.abf"], PP.H1 = s["PP.H1.abf"], PP.H2 = s["PP.H2.abf"],
    PP.H3 = s["PP.H3.abf"], PP.H4 = s["PP.H4.abf"],
    top_SNP = top_snp$snp, top_SNP_PP.H4 = top_snp$SNP.PP.H4,
    stringsAsFactors = FALSE)
  write.csv(snp_res, file.path(outdir, paste0("coloc_snps_", gene, ".csv")), row.names = FALSE)
}

if (length(results) > 0) {
  all <- do.call(rbind, results)
  write.csv(all, file.path(outdir, "coloc_summary_all_genes.csv"), row.names = FALSE)
  cat("===== Coloc 汇总 =====\n")
  print(all[, c("gene","nsnp_common","PP.H1","PP.H3","PP.H4")])
}
cat("\n完成。目录:", outdir, "\n")
