# ============================================================================
# MR 结果整合可视化: 10基因 x 4方法 森林图 + 散点图汇总
# ============================================================================
suppressMessages({
  library(ggplot2)
  library(dplyr)
})

outdir <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 4 免疫浸润与诊断模型/MR"
df <- read.csv(file.path(outdir, "MR_summary_all_genes.csv"), stringsAsFactors = FALSE)

# 标准化方法名
df$method2 <- case_when(
  grepl("Inverse variance", df$method) ~ "IVW",
  grepl("Egger", df$method) ~ "MR-Egger",
  grepl("Weighted median", df$method) ~ "Weighted median",
  grepl("Weighted mode", df$method) ~ "Weighted mode",
  TRUE ~ df$method
)

# 计算 OR 和 CI
df$OR <- exp(df$b)
df$OR_low <- exp(df$b - 1.96 * df$se)
df$OR_high <- exp(df$b + 1.96 * df$se)
df$sig <- ifelse(df$pval < 0.05, "P<0.05", ifelse(df$pval < 0.1, "P<0.1", "ns"))

# 按 IVW 排序
ivw_order <- df %>% filter(method2 == "IVW") %>% arrange(pval) %>% pull(gene)
df$gene <- factor(df$gene, levels = rev(ivw_order))
df$method2 <- factor(df$method2, levels = c("IVW", "MR-Egger", "Weighted median", "Weighted mode"))

# ---- 森林图 ----
p <- ggplot(df, aes(x = OR, y = gene, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.3, position = position_dodge(width = 0.6)) +
  facet_wrap(~method2, ncol = 4) +
  scale_color_manual(values = c("P<0.05" = "#D7191C", "P<0.1" = "#FDAE61", "ns" = "grey60")) +
  scale_x_log10() +
  labs(title = "Mendelian Randomization: 17 Cell Death Genes -> CAD (Nikpay 2015)",
       subtitle = "eQTLGen exposure | ebi-a-GCST005194 outcome",
       x = "OR (95% CI) per SD of gene expression", y = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(outdir, "MR_forest_all_methods.png"), p, width = 12, height = 7, dpi = 300)
ggsave(file.path(outdir, "MR_forest_all_methods.svg"), p, width = 12, height = 7)
cat("✅ 森林图已保存\n")

# ---- IVW 明细表 ----
ivw <- df %>% filter(method2 == "IVW") %>% arrange(pval)
cat("\n===== IVW 结果排序 =====")
print(ivw[, c("gene", "nsnp", "OR", "OR_low", "OR_high", "pval")])

# 保存整理后的表
write.csv(ivw[, c("gene", "nsnp", "b", "se", "pval", "OR", "OR_low", "OR_high")],
          file.path(outdir, "MR_IVW_summary_table.csv"), row.names = FALSE)

# ---- 显著结果 (加权中位数/众数) ----
cat("\n===== 敏感性分析显著 (P<0.05) =====")
sig <- df %>% filter(pval < 0.05) %>% select(gene, method2, nsnp, b, se, pval, OR)
if (nrow(sig) > 0) print(sig) else cat("无\n")
