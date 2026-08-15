
# ============================================================
# Pseudotime trajectory analysis (Supplementary Figure S8)
# Monocle 3 trajectory of myeloid cells + cell-death scoring
# ============================================================
suppressMessages({
    library(Seurat)
    library(monocle3)
    library(dplyr)
    library(ggplot2)
})

so <- readRDS("Results/seurat_processed.rds")

# 1) Cell-death module scores (AddModuleScore)
pyro_genes <- c("AIM2", "CASP1", "IL18", "NLRP3", "PYCARD")
neto_genes <- c("CTSG", "CYBA", "ITGAM", "MPO", "NCF2")
so <- AddModuleScore(so, features = list(pyro_genes), name = "PyroptosisScore", ctrl = 50)
so <- AddModuleScore(so, features = list(neto_genes), name = "NETosisScore", ctrl = 50)

summary_df <- data.frame(
    celltype = so$celltype,
    pyro = so$PyroptosisScore1,
    neto = so$NETosisScore1
) %>% group_by(celltype) %>%
    summarise(n_cells = n(), Pyroptosis_mean = mean(pyro), NETosis_mean = mean(neto)) %>%
    arrange(desc(Pyroptosis_mean))
write.csv(as.data.frame(summary_df), "Results/death_scores_summary.csv", row.names = FALSE)

# 2) Myeloid trajectory
myeloid <- subset(so, celltype %in% c("Macrophage", "Monocyte"))
cds <- SeuratWrappers::as.cell_data_set(myeloid)
cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE)

root_idx <- which.max(myeloid$PyroptosisScore1)
cds <- order_cells(cds, root_cells = colnames(myeloid)[root_idx])
pseudotime <- pseudotime(cds)
write.csv(data.frame(cell = names(pseudotime), pseudotime = as.numeric(pseudotime)),
          "Results/pseudotime_myeloid.csv", row.names = FALSE)

pdf("Results/pseudotime_trajectory_v2.pdf", width = 9, height = 6.5)
print(plot_cells(cds, color_cells_by = "pseudotime",
                 label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE,
                 cell_size = 1.2) +
      ggtitle("Monocyte-macrophage trajectory inferred by Monocle 3") +
      scale_color_viridis_c(option = "plasma"))
dev.off()

# 3) Death scores along pseudotime
df_plot <- data.frame(
    pseudotime = as.numeric(pseudotime),
    pyroptosis = myeloid$PyroptosisScore1,
    netosis = myeloid$NETosisScore1,
    celltype = myeloid$celltype
)
write.csv(df_plot, "Results/pseudotime_death_scores.csv", row.names = FALSE)

pdf("Results/pseudotime_death_scores.pdf", width = 11, height = 4.5)
p1 <- ggplot(df_plot, aes(x = pseudotime, y = pyroptosis, color = celltype)) +
    geom_point(size = 0.4, alpha = 0.5) + geom_smooth(se = FALSE, linewidth = 1) +
    labs(title = "Pyroptosis score along myeloid pseudotime", y = "Pyroptosis score") + theme_bw()
p2 <- ggplot(df_plot, aes(x = pseudotime, y = netosis, color = celltype)) +
    geom_point(size = 0.4, alpha = 0.5) + geom_smooth(se = FALSE, linewidth = 1) +
    labs(title = "NETosis score along myeloid pseudotime", y = "NETosis score") + theme_bw()
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()
