# ============================================================
# Step 2: Load RDS → limma DE → GO enrichment → 环形热图+UpSet
# ============================================================
suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(readxl)
  library(writexl)
  library(RColorBrewer)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(svglite)
})

# ---------- 路径 ----------
run_args <- commandArgs(FALSE)
file_arg <- run_args[grepl("^--file=", run_args)]
if (length(file_arg) > 0) {
  out_dir <- normalizePath(dirname(sub("^--file=", "", file_arg[1])))
} else {
  out_dir <- normalizePath(".")
}
setwd(out_dir)

rds_file <- file.path(out_dir, "expr_norm.rds")
dat <- readRDS(rds_file)
expr_norm <- dat$expr_norm
gene_map  <- dat$gene_map
group_df  <- dat$group_df

cat("Loaded:", nrow(expr_norm), "probes x", ncol(expr_norm), "samples\n")
cat("Groups:", table(group_df$Group), "\n")

# ---------- limma DE ----------
group <- factor(group_df$Group, levels = c("Control", "Atherosclerotic"))
design <- model.matrix(~ group); colnames(design) <- c("Intercept", "Athero_vs_Control")
fit <- lmFit(expr_norm, design)
fit <- eBayes(fit, trend = TRUE)

de_all <- topTable(fit, coef = "Athero_vs_Control", number = Inf, sort.by = "P")
de_all$ProbeName <- rownames(de_all)
de_all$SYMBOL <- gene_map[de_all$ProbeName]
empty_gene <- is.na(de_all$SYMBOL) | de_all$SYMBOL == ""
de_all$SYMBOL[empty_gene] <- de_all$ProbeName[empty_gene]

fc_threshold <- 0.5; fdr_threshold <- 0.05
up_genes <- unique(de_all$SYMBOL[de_all$logFC >  fc_threshold & de_all$adj.P.Val < fdr_threshold])
dn_genes <- unique(de_all$SYMBOL[de_all$logFC < -fc_threshold & de_all$adj.P.Val < fdr_threshold])
cat("UP:", length(up_genes), "DN:", length(dn_genes), "\n")

# ---------- 选代表性样本 ----------
set.seed(42)
ctrl_idx <- which(group_df$Group == "Control")
ath_idx  <- which(group_df$Group == "Atherosclerotic")
selected_ctrl <- sample(colnames(expr_norm)[ctrl_idx], min(3, length(ctrl_idx)))
selected_ath  <- sample(colnames(expr_norm)[ath_idx], min(3, length(ath_idx)))
selected_samples <- c(selected_ctrl, selected_ath)
sample_group_vec <- setNames(c(rep("Control",length(selected_ctrl)), rep("Atherosclerotic",length(selected_ath))), selected_samples)
cat("Selected cells:", paste(selected_samples, collapse=", "), "\n")

# ---------- GO 富集 ----------
up_entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
dn_entrez <- bitr(dn_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up <- enrichGO(gene=up_entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP",
                   pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.1)
ego_dn <- enrichGO(gene=dn_entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP",
                   pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.1)

go_combined <- rbind(
  if (!is.null(ego_up) && nrow(ego_up@result) > 0) ego_up@result %>% mutate(direction="UP") else NULL,
  if (!is.null(ego_dn) && nrow(ego_dn@result) > 0) ego_dn@result %>% mutate(direction="DN") else NULL
)

if (nrow(go_combined) > 0) {
  go_combined <- go_combined %>% arrange(p.adjust)
  top_terms <- go_combined %>% slice_head(n=13)
  selected_terms <- top_terms$ID
  term_names <- sapply(top_terms$Description, function(x) {
    if (nchar(x) > 35) paste0(substr(x,1,32),"...") else x
  })
  term_directions <- top_terms$direction
  cat("GO terms selected:", length(term_names), "\n")
  for (i in seq_along(term_names)) cat(sprintf("  %s [%s]: %s\n", i, term_directions[i], term_names[i]))
} else {
  stop("GO enrichment returned no terms")
}

# ---------- 提取每个 term 的 top 基因 + 表达 ----------
top_n_per_group <- 30
heatmap_rows <- list()

for (i in seq_along(selected_terms)) {
  tid <- selected_terms[i]
  tname <- term_names[i]
  gene_list <- strsplit(top_terms$geneID[i], "/")[[1]]
  if (length(gene_list) == 0) next
  gene_sym <- bitr(gene_list, fromType="ENTREZID", toType="SYMBOL", OrgDb=org.Hs.eg.db)
  if (nrow(gene_sym) == 0) next
  
  # 仅保留 DEG
  deg_genes <- de_all[de_all$SYMBOL %in% gene_sym$SYMBOL & de_all$adj.P.Val < 0.05, ]
  if (nrow(deg_genes) == 0) next
  
  if (term_directions[i] == "DN") {
    deg_genes <- deg_genes %>% arrange(adj.P.Val, logFC)
  } else {
    deg_genes <- deg_genes %>% arrange(adj.P.Val, desc(logFC))
  }
  top_genes <- deg_genes %>% slice_head(n=top_n_per_group)
  if (nrow(top_genes) == 0) next
  
  for (g in top_genes$SYMBOL) {
    probes <- de_all$ProbeName[de_all$SYMBOL == g]
    if (length(probes) == 0) next
    p <- probes[1]
    heatmap_rows[[paste0(tname,"|",g)]] <- c(tname, g, as.numeric(expr_norm[p, selected_samples]))
  }
}

# 构建 data frame
if (length(heatmap_rows) == 0) stop("No heatmap rows generated")
heatmap_df <- do.call(rbind, heatmap_rows)
n_samps <- length(selected_samples)
col_names <- paste0("Exp_", 1:n_samps)

out_df <- data.frame(
  Group = heatmap_df[,1],
  Gene  = heatmap_df[,2],
  stringsAsFactors = FALSE
)
for (j in 1:n_samps) {
  out_df[[col_names[j]]] <- as.numeric(heatmap_df[, j+2])
}

# Z-score by row
expr_mat <- as.matrix(out_df[, col_names])
expr_z <- t(scale(t(expr_mat)))
expr_z[is.na(expr_z)] <- 0
expr_z[expr_z > 2] <- 2; expr_z[expr_z < -2] <- -2
out_df[, col_names] <- round(expr_z, 4)
out_df <- out_df[rowSums(is.na(out_df[,col_names])) == 0, ]
cat("Heatmap data:", nrow(out_df), "genes x", n_samps, "samples\n")

xlsx_file <- file.path(out_dir, "data_for_upset.xlsx")
write_xlsx(out_df, xlsx_file)
cat("XLSX saved:", xlsx_file, "\n")

# ============================================================
# 绘图（下面直接嵌入热图.R 的逻辑）
# ============================================================

data <- out_df
group_order <- unique(data$Group)
data$Group <- factor(data$Group, levels = group_order)
ng <- length(group_order)

data_matrix <- as.matrix(data[, col_names])
rownames(data_matrix) <- data$Gene

sample_group <- sample_group_vec[selected_samples]
colnames(data_matrix) <- col_names
names(sample_group) <- col_names

# ---- 颜色 ----
red_blue <- colorRamp2(
  seq(-2, 2, length.out = 100),
  colorRampPalette(brewer.pal(n = 9, name = "RdBu"))(100)
)

# 13 色（复用框架配色表，按需截取）
all_13_colors <- c("#D73027","#FC8D59","#FEE090","#91BFDB","#4575B4","#313695",
                   "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF")
group_colors <- setNames(all_13_colors[1:ng], group_order)

sample_colors <- c("Control" = "#66C2A5", "Atherosclerotic" = "#FC8D62")

# ---- 环形热图 ----
circos.clear()
circos.par(
  start.degree = 60,
  gap.after    = c(rep(2, ng - 1), 15),
  track.margin = c(0, 0.01),
  cell.padding = c(0, 0, 0, 0)
)

circos.heatmap(
  data_matrix,
  split          = data$Group,
  cluster        = FALSE,
  bg.border      = "black",
  bg.lwd         = 0.8,
  cell.border    = "white",
  cell.lwd       = 0.3,
  rownames.side  = "outside",
  rownames.cex   = 0.35,
  col            = red_blue,
  track.height   = 0.20
)

# Layer 2: 列名
circos.track(
  track.index = get.current.track.index(),
  bg.border   = NA,
  panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == ng) {
      cn <- colnames(data_matrix)
      n  <- length(cn)
      ch <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
      yc <- seq(CELL_META$cell.ylim[1] + ch/2,
                CELL_META$cell.ylim[2] - ch/2, length.out = n)
      for (i in 1:n) {
        circos.lines(
          c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")),
          c(yc[i], yc[i]), col = "black", lwd = 1.5
        )
      }
      circos.text(
        rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"),
        yc, rev(cn), cex = 1, adj = c(0, 0.5), facing = "inside"
      )
    }
  }
)

# Layer 3: 样本分组色带
circos.track(
  track.index = get.current.track.index(),
  bg.border   = NA,
  panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == ng) {
      n <- length(sample_group)
      ch <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
      yc <- seq(CELL_META$cell.ylim[1] + ch/2,
                CELL_META$cell.ylim[2] - ch/2, length.out = n)
      for (i in 1:n) {
        circos.rect(
          CELL_META$cell.xlim[2] + convert_x(4, "mm"),
          yc[i] - ch * 0.35,
          CELL_META$cell.xlim[2] + convert_x(5.5, "mm"),
          yc[i] + ch * 0.35,
          col = sample_colors[sample_group[i]], border = NA
        )
      }
    }
  }
)

# Layer 4: 扇区分组标签
circos.track(
  ylim         = c(0, 1),
  track.height = 0.08,
  bg.col       = adjustcolor(group_colors[levels(data$Group)], alpha.f = 0.3),
  bg.border    = NA,
  panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, 0.5, CELL_META$sector.index,
                facing = "bending.inside", niceFacing = TRUE,
                cex = 0.55, adj = c(0.5, 0.5))
  }
)

# ---- 图例 ----
heatmap_legend <- Legend(
  title = "Expression (Z-score)",
  col_fun = red_blue,
  at = seq(-2, 2, length.out = 6),
  title_position = "leftcenter-rot",
  title_gp = gpar(fontsize = 12),
  labels_gp = gpar(fontsize = 12)
)
draw(heatmap_legend,
     x = unit(0.95, "npc") - unit(3, "mm"),
     y = unit(0.95, "npc") - unit(3, "mm"),
     just = c("right", "top"))

sample_legend <- Legend(
  title = "Sample",
  labels = names(sample_colors),
  legend_gp = gpar(fill = sample_colors),
  title_gp = gpar(fontsize = 10),
  labels_gp = gpar(fontsize = 10)
)
draw(sample_legend,
     x = unit(0.95, "npc") - unit(3, "mm"),
     y = unit(0.80, "npc"),
     just = c("right", "top"))

# ---- UpSet ----
gene_sets_list <- lapply(group_order, function(g) {
  unique(data %>% filter(Group == g) %>% pull(Gene))
})
names(gene_sets_list) <- group_order

m <- make_comb_mat(gene_sets_list)
comb_sizes <- comb_size(m)
valid <- sort(comb_sizes[comb_sizes > 0], decreasing = TRUE)
top_n_upset <- min(15, length(valid))

if (top_n_upset > 0) {
  m_sub <- m[names(valid)[1:top_n_upset]]
  pushViewport(viewport(x = 0.5, y = 0.5, width = 0.44, height = 0.28))
  ht_upset <- UpSet(
    m_sub,
    set_order   = group_order,
    comb_col    = "#4393C3",
    pt_size     = unit(1.3, "mm"),
    lwd         = 0.6,
    left_annotation = upset_left_annotation(
      m_sub,
      gp = gpar(fill = adjustcolor(group_colors[names(gene_sets_list)], alpha.f = 0.3)),
      annotation_name_side = "bottom"
    ),
    right_annotation = NULL,
    top_annotation = upset_top_annotation(m_sub, gp = gpar(fill = "#4393C3")),
    row_names_gp = gpar(fontsize = 5, fontfamily = "serif"),
    column_title = "Intersection size",
    column_title_gp = gpar(fontsize = 7, fontfamily = "serif")
  )
  draw(ht_upset, newpage = FALSE, background = "transparent")
  grid.text("GSE100927\nPeripheral Artery\nAtherosclerosis",
            x = unit(0.5, "npc"), y = unit(-0.06, "npc"),
            gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "serif"))
  popViewport()
}

# ---- 保存 SVG ----
svg_file <- file.path(out_dir, "GSE100927_CircularHeatmap.svg")
svglite(svg_file, width = 12, height = 10, bg = "white")
# 需要重绘（当前设备是交互式，需要完整重画）
circos.clear()
circos.par(start.degree = 60, gap.after = c(rep(2, ng-1), 15),
           track.margin = c(0, 0.01), cell.padding = c(0, 0, 0, 0))
circos.heatmap(data_matrix, split = data$Group, cluster = FALSE,
               bg.border = "black", bg.lwd = 0.8, cell.border = "white",
               cell.lwd = 0.3, rownames.side = "outside", rownames.cex = 0.35,
               col = red_blue, track.height = 0.20)
circos.track(track.index = get.current.track.index(), bg.border = NA,
  panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == ng) {
      cn <- colnames(data_matrix); n <- length(cn)
      ch <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
      yc <- seq(CELL_META$cell.ylim[1] + ch/2, CELL_META$cell.ylim[2] - ch/2, length.out = n)
      for (i in 1:n) circos.lines(c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1,"mm")), c(yc[i], yc[i]), col="black", lwd=1.5)
      circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1.5,"mm"), yc, rev(cn), cex=1, adj=c(0,0.5), facing="inside")
    }
  })
circos.track(track.index = get.current.track.index(), bg.border = NA,
  panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == ng) {
      n <- length(sample_group); ch <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
      yc <- seq(CELL_META$cell.ylim[1] + ch/2, CELL_META$cell.ylim[2] - ch/2, length.out = n)
      for (i in 1:n) circos.rect(CELL_META$cell.xlim[2] + convert_x(4,"mm"), yc[i] - ch*0.35, CELL_META$cell.xlim[2] + convert_x(5.5,"mm"), yc[i] + ch*0.35, col = sample_colors[sample_group[i]], border = NA)
    }
  })
circos.track(ylim=c(0,1), track.height=0.08, bg.col=adjustcolor(group_colors[levels(data$Group)], alpha.f=0.3), bg.border=NA,
  panel.fun = function(x,y) circos.text(CELL_META$xcenter, 0.5, CELL_META$sector.index, facing="bending.inside", niceFacing=TRUE, cex=0.55, adj=c(0.5,0.5)))
draw(heatmap_legend, x=unit(0.95,"npc")-unit(3,"mm"), y=unit(0.95,"npc")-unit(3,"mm"), just=c("right","top"))
draw(sample_legend, x=unit(0.95,"npc")-unit(3,"mm"), y=unit(0.80,"npc"), just=c("right","top"))
if (top_n_upset > 0) {
  pushViewport(viewport(x=0.5, y=0.5, width=0.44, height=0.28))
  draw(ht_upset, newpage=FALSE, background="transparent")
  grid.text("GSE100927\nPeripheral Artery\nAtherosclerosis", x=unit(0.5,"npc"), y=unit(-0.06,"npc"), gp=gpar(fontsize=10, fontface="bold", fontfamily="serif"))
  popViewport()
}
dev.off()
cat("SVG saved:", svg_file, "\n")

# ---- 保存 PDF ----
pdf_file <- file.path(out_dir, "GSE100927_CircularHeatmap.pdf")
pdf(pdf_file, width = 12, height = 10)
# 重绘（circos 不同步，需完整重画）
circos.clear()
circos.par(start.degree = 60, gap.after = c(rep(2, ng-1), 15),
           track.margin = c(0, 0.01), cell.padding = c(0, 0, 0, 0))
circos.heatmap(data_matrix, split = data$Group, cluster = FALSE,
               bg.border = "black", bg.lwd = 0.8, cell.border = "white",
               cell.lwd = 0.3, rownames.side = "outside", rownames.cex = 0.35,
               col = red_blue, track.height = 0.20)
circos.track(track.index = get.current.track.index(), bg.border = NA,
  panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == ng) {
      cn <- colnames(data_matrix); n <- length(cn)
      ch <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
      yc <- seq(CELL_META$cell.ylim[1] + ch/2, CELL_META$cell.ylim[2] - ch/2, length.out = n)
      for (i in 1:n) circos.lines(c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1,"mm")), c(yc[i], yc[i]), col="black", lwd=1.5)
      circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1.5,"mm"), yc, rev(cn), cex=1, adj=c(0,0.5), facing="inside")
    }
  })
circos.track(track.index = get.current.track.index(), bg.border = NA,
  panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == ng) {
      n <- length(sample_group); ch <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
      yc <- seq(CELL_META$cell.ylim[1] + ch/2, CELL_META$cell.ylim[2] - ch/2, length.out = n)
      for (i in 1:n) circos.rect(CELL_META$cell.xlim[2] + convert_x(4,"mm"), yc[i] - ch*0.35, CELL_META$cell.xlim[2] + convert_x(5.5,"mm"), yc[i] + ch*0.35, col = sample_colors[sample_group[i]], border = NA)
    }
  })
circos.track(ylim=c(0,1), track.height=0.08, bg.col=adjustcolor(group_colors[levels(data$Group)], alpha.f=0.3), bg.border=NA,
  panel.fun = function(x,y) circos.text(CELL_META$xcenter, 0.5, CELL_META$sector.index, facing="bending.inside", niceFacing=TRUE, cex=0.55, adj=c(0.5,0.5)))
draw(heatmap_legend, x=unit(0.95,"npc")-unit(3,"mm"), y=unit(0.95,"npc")-unit(3,"mm"), just=c("right","top"))
draw(sample_legend, x=unit(0.95,"npc")-unit(3,"mm"), y=unit(0.80,"npc"), just=c("right","top"))
if (top_n_upset > 0) {
  pushViewport(viewport(x=0.5, y=0.5, width=0.44, height=0.28))
  draw(ht_upset, newpage=FALSE, background="transparent")
  grid.text("GSE100927\nPeripheral Artery\nAtherosclerosis", x=unit(0.5,"npc"), y=unit(-0.06,"npc"), gp=gpar(fontsize=10, fontface="bold", fontfamily="serif"))
  popViewport()
}
dev.off()
cat("PDF saved:", pdf_file, "\n")
cat("=== Done ===\n")
