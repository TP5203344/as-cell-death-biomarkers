suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
})

# 加载已有 DE 结果
de <- read.csv("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 1 研究设计与转录组差异分析/火山图GSE100927/GSE100927_DE_results.csv",
               stringsAsFactors = FALSE)

# 提取显著 DEG 基因名
up_genes <- unique(de$external_gene_name[de$logFC > 0.5 & de$FDR < 0.05])
dn_genes <- unique(de$external_gene_name[de$logFC < -0.5 & de$FDR < 0.05])

up_genes <- up_genes[!is.na(up_genes) & up_genes != ""]
dn_genes <- dn_genes[!is.na(dn_genes) & dn_genes != ""]

cat("UP genes:", length(up_genes), "DN genes:", length(dn_genes), "\n")

# 双向合并，写入 gene.txt
all_genes <- unique(c(up_genes, dn_genes))
out <- data.frame(id = all_genes, stringsAsFactors = FALSE)

out_dir <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 1 研究设计与转录组差异分析/桑基气泡图GSE100927"
write.table(out, file.path(out_dir, "gene.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("gene.txt saved:", nrow(out), "genes\n")
