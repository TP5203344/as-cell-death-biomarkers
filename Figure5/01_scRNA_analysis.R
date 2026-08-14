# ============================================================================
# Figure 5 单细胞转录组定位核心基因 (GSE115469) - 修复注释逻辑
# 上次 UMAP/聚类已完成(21 clusters)，本次只修复注释与下游，从保存的RDS继续
# ============================================================================
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(pheatmap)
})

setwd("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 5 单细胞定位")
dir.create("Results", showWarnings = FALSE)

# ---- 从RDS恢复（若存在）----
if (file.exists("Results/seurat_processed.rds")) {
  cat("从 RDS 恢复 Seurat 对象...\n")
  so <- readRDS("Results/seurat_processed.rds")
} else {
  cat("重新运行全流程...\n")
  expr <- read.csv("C:/Users/19267/.qclaw/workspace-agent-b69b3ebd/as_bioinfo/_gse115469.csv.gz",
                   row.names = 1, check.names = FALSE)
  so <- CreateSeuratObject(counts = expr, project = "GSE115469", min.cells = 3, min.features = 200)
  so$sample <- sapply(strsplit(colnames(so), "_"), function(x) x[1])
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
  so <- NormalizeData(so)
  so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
  so <- ScaleData(so)
  so <- RunPCA(so, features = VariableFeatures(so), npcs = 30)
  so <- RunUMAP(so, dims = 1:20)
  so <- FindNeighbors(so, dims = 1:20)
  so <- FindClusters(so, resolution = 0.6)
}

cat("聚类数:", length(unique(so$seurat_clusters)), "\n")

# ---- UMAP 图 ----
png("Results/02_UMAP_clusters.png", width = 2200, height = 1800, res = 300)
print(DimPlot(so, reduction = "umap", label = TRUE, pt.size = 0.5) +
        labs(title = "UMAP - GSE115469 Clusters"))
dev.off()

png("Results/03_UMAP_sample.png", width = 2200, height = 1800, res = 300)
print(DimPlot(so, reduction = "umap", group.by = "sample", pt.size = 0.5) +
        labs(title = "UMAP by Sample"))
dev.off()

# ---- 细胞类型注释（修正：矩阵行=基因，列=cluster）----
markers <- list(
  Macrophage = c("CD68", "CD163", "MRC1", "MSR1"),
  Monocyte   = c("CD14", "LYZ", "S100A8", "S100A9"),
  Endothelial= c("VWF", "PECAM1", "CLDN5"),
  SMC        = c("ACTA2", "MYH11", "TAGLN"),
  T_cell     = c("CD3D", "CD3E", "CD2"),
  B_cell     = c("CD79A", "MS4A1", "CD19"),
  NK         = c("NKG7", "KLRD1", "GNLY")
)

cluster_means <- AverageExpression(so, features = unique(unlist(markers)), group.by = "seurat_clusters")$RNA
cat("AverageExpression 维度:", nrow(cluster_means), "基因 x", ncol(cluster_means), "cluster\n")
write.csv(cluster_means, "Results/04_cluster_marker_expression.csv")

# 每个 cluster：计算 7 类 marker 的平均表达，取最高类
clusters <- colnames(cluster_means)
annotation <- character(length(clusters))
names(annotation) <- clusters
for (cl in clusters) {
  scores <- sapply(markers, function(m) {
    hit <- intersect(m, rownames(cluster_means))
    if (length(hit) == 0) return(0)
    mean(cluster_means[hit, cl])
  })
  annotation[cl] <- ifelse(max(scores) == 0, "Unknown", names(which.max(scores)))
}
cat("\n细胞类型注释:\n")
print(table(annotation))

# cluster id 可能带 g 前缀 (如 g0)，映射回 seurat_clusters
ct_vec <- annotation[paste0("g", as.character(so$seurat_clusters))]
names(ct_vec) <- colnames(so)
so <- AddMetaData(so, metadata = ct_vec, col.name = "celltype")
so$celltype[is.na(so$celltype)] <- "Unknown"

write.csv(data.frame(cluster = names(annotation), celltype = unname(annotation)),
          "Results/05_cluster_annotation.csv", row.names = FALSE)

png("Results/06_UMAP_celltype.png", width = 2400, height = 1800, res = 300)
print(DimPlot(so, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 0.5) +
        labs(title = "UMAP - Cell Types"))
dev.off()

# ---- 核心基因定位 ----
core_genes <- c("AIM2", "CASP1", "NLRP3", "RIPK3", "PYCARD", "IL18",
                "ITGAM", "MPO", "NCF2", "CTSG", "CTSD", "PGD",
                "ZEB1", "PRKAA2", "BMF", "ITGA5", "CYBA")
available <- intersect(core_genes, rownames(so))
cat("\n核心基因在单细胞中命中:", length(available), "/", length(core_genes), "\n")
cat("命中:", paste(available, collapse = ", "), "\n")

if (length(available) > 0) {
  png("Results/07_DotPlot_core_genes.png", width = 2800, height = 1200, res = 300)
  print(DotPlot(so, features = available, group.by = "celltype") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(title = "Core Cell Death Genes - Cell Type Specificity"))
  dev.off()

  for (g in available) {
    png(paste0("Results/08_FeaturePlot_", g, ".png"), width = 1800, height = 1500, res = 300)
    print(FeaturePlot(so, features = g, pt.size = 0.4) + labs(title = g))
    dev.off()
  }

  avg <- AverageExpression(so, features = available, group.by = "celltype")$RNA
  png("Results/09_Heatmap_core_celltype.png", width = 2600, height = 1600, res = 300)
  pheatmap(as.matrix(avg), scale = "row", cluster_cols = TRUE,
           main = "Core Genes by Cell Type",
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
  dev.off()
}

# ---- M1/M2 巨噬极化评分 ----
m1_genes <- intersect(c("NOS2", "IL1B", "TNF", "CXCL9", "CXCL10", "CXCL11", "CD80", "CD86"), rownames(so))
m2_genes <- intersect(c("MRC1", "CD163", "MSR1", "IL10", "TGFB1", "ARG1", "CCL18", "CCL22"), rownames(so))
cat("M1 marker 命中:", paste(m1_genes, collapse=","), "\n")
cat("M2 marker 命中:", paste(m2_genes, collapse=","), "\n")

so$M1_score <- colMeans(GetAssayData(so, layer = "data")[m1_genes, , drop = FALSE])
so$M2_score <- colMeans(GetAssayData(so, layer = "data")[m2_genes, , drop = FALSE])

png("Results/10_M1_M2_scores.png", width = 2600, height = 1000, res = 300)
print(FeaturePlot(so, features = c("M1_score", "M2_score"), pt.size = 0.4) &
        scale_color_gradientn(colors = c("grey90", "orange", "red")))
dev.off()

# ---- 细胞类型组成占比饼图/条形图 ----
celltype_counts <- table(so$celltype)
write.csv(as.data.frame(celltype_counts), "Results/11_celltype_composition.csv")

png("Results/12_celltype_barplot.png", width = 2200, height = 1400, res = 300)
print(ggplot(data.frame(celltype = so$celltype), aes(x = celltype, fill = celltype)) +
        geom_bar() + labs(title = "Cell Type Composition", x = "", y = "Cells") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none"))
dev.off()

saveRDS(so, "Results/seurat_processed.rds")
cat("\n✅ Figure 5 单细胞分析全部完成\n")
cat("结果目录: Results/\n")
