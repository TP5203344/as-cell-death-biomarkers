# ============================================================================
# Figure 5 补充: 核心基因在细胞类型中的 VlnPlot + 细胞类型占比饼图
# 从已保存对象恢复 (seurat_processed.rds 在上次失败前未保存, 重新加载矩阵+注释)
# ============================================================================
suppressMessages({
  library(Seurat)
  library(ggplot2)
})

setwd("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 5 单细胞定位")
dir.create("Results", showWarnings = FALSE)

# 读取矩阵并重建
cat("读取表达矩阵...\n")
expr <- read.csv("C:/Users/19267/.qclaw/workspace-agent-b69b3ebd/as_bioinfo/_gse115469.csv.gz",
                 row.names = 1, check.names = FALSE)
so <- CreateSeuratObject(counts = expr, project = "GSE115469", min.cells = 3, min.features = 200)
so$sample <- sapply(strsplit(colnames(so), "_"), function(x) x[1])
so <- NormalizeData(so)

# 读取注释并赋值 (cluster -> celltype)
ann <- read.csv("Results/05_cluster_annotation.csv")
# 重新计算聚类 (resolution 0.6 结果应与之前一致, 但需要重新跑 PCA/UMAP/cluster)
cat("重新聚类...\n")
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
so <- ScaleData(so)
so <- RunPCA(so, features = VariableFeatures(so), npcs = 30)
so <- RunUMAP(so, dims = 1:20)
so <- FindNeighbors(so, dims = 1:20)
so <- FindClusters(so, resolution = 0.6)

# 用注释映射
ct_map <- setNames(ann$celltype, ann$cluster)  # g0, g1...
ids <- paste0("g", as.character(so$seurat_clusters))
ct_vec <- unname(ct_map[ids])
ct_vec[is.na(ct_vec)] <- "Unknown"
names(ct_vec) <- colnames(so)
so <- AddMetaData(so, metadata = ct_vec, col.name = "celltype")

# ---- 细胞类型占比 ----
compo <- as.data.frame(table(so$celltype))
colnames(compo) <- c("celltype", "count")
compo$pct <- round(100 * compo$count / sum(compo$count), 1)
write.csv(compo, "Results/11_celltype_composition.csv", row.names = FALSE)
cat("\n细胞类型占比:\n")
print(compo)

png("Results/12_celltype_barplot.png", width = 2200, height = 1400, res = 300)
print(ggplot(compo, aes(x = reorder(celltype, -count), y = count, fill = celltype)) +
        geom_col() + geom_text(aes(label = paste0(pct, "%")), vjust = -0.3) +
        labs(title = "Cell Type Composition (GSE115469)", x = "", y = "Cells") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none"))
dev.off()

# ---- 核心基因 VlnPlot (分组按细胞类型, 选取关键基因) ----
core_sel <- c("AIM2", "CASP1", "NLRP3", "RIPK3", "PYCARD", "IL18",
              "ITGAM", "MPO", "NCF2", "CTSG", "CTSD", "ZEB1")
available <- intersect(core_sel, rownames(so))
cat("VlnPlot 基因命中:", length(available), "\n")

# 分 3 组图 (每组4个基因)
for (i in 1:ceiling(length(available)/4)) {
  idx <- ((i-1)*4+1):min(i*4, length(available))
  png(paste0("Results/13_VlnPlot_core_genes_", i, ".png"), width = 3200, height = 1800, res = 300)
  print(VlnPlot(so, features = available[idx], group.by = "celltype", pt.size = 0,
                ncol = 2) & theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

# ---- UMAP 按细胞类型 (重新生成确保与注释一致) ----
png("Results/06_UMAP_celltype.png", width = 2400, height = 1800, res = 300)
print(DimPlot(so, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 0.5) +
        labs(title = "UMAP - Cell Types"))
dev.off()

# ---- 保存完整对象 ----
saveRDS(so, "Results/seurat_processed.rds")
cat("\n✅ 补充分析完成, 对象已保存\n")
