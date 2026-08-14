# =============================================================================
# GSE100927 富集分析综合桑基气泡图 - GO(BP/CC/MF), DO, KEGG, Reactome
# 套用框架:桑吉气泡图.R
# 数据源:火山图GSE100927/GSE100927_DE_results.csv 的 DEG 基因列表
# =============================================================================

# =============================================================================
# 设置工作目录
# =============================================================================
base_dir <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 1 研究设计与转录组差异分析/桑基气泡图GSE100927"
setwd(base_dir)

# =============================================================================
# 数据处理与基因ID转换
# =============================================================================
rt <- read.table("gene.txt", header = TRUE, sep = "\t", check.names = FALSE)
colnames(rt)[1] <- "id"

suppressPackageStartupMessages({
  library(magick)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(tidyverse)
  library(ggalluvial)
  library(patchwork)
  library(RColorBrewer)
  library(Cairo)
  library(DOSE)
  library(GOSemSim)
  library(AnnotationDbi)
})

# =============================================================================
# 预加载 HDO 数据库(绕过 yulab-smu.top 防火墙)
# =============================================================================
hdo_file <- "C:/Users/19267/AppData/Local/GOSemSim/HDO.sqlite"
if (!file.exists(hdo_file)) {
  stop("HDO.sqlite not found at ", hdo_file)
}
db_hdo <- loadDb(hdo_file)
assign(".onto_HDO", db_hdo, envir = GOSemSim:::get_gosemsim_env())
cat("HDO pre-loaded successfully\n")

# =============================================================================
# 参数设置
# =============================================================================
pvalueFilter  <- 0.05
adjPvalFilter <- 0.05
showNum       <- 30
top_pathways  <- 20
colorSel <- "p.adjust"
if (adjPvalFilter > 0.05) colorSel <- "pvalue"

plot_width         <- 12
plot_height        <- 14
sankey_width_ratio <- 3
dot_width_ratio    <- 2
gene_y_pos         <- 1.1
pathway_y_pos      <- 1
sankey_text_size   <- 2.8
bubble_x_label     <- "GeneRatio"
bubble_p_label     <- "Pvalue"

pathway_colors    <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
bubble_low_color  <- "#2166ac"
bubble_high_color <- "#b2182b"
gene_palette_name <- "Set3"
bubble_size_range <- c(3, 8)
bubble_legend_size <- c(2, 4, 6, 8, 10)

windowsFonts(Arial = windowsFont("Arial"))
plot_font <- "Arial"

# ID转换
genes <- unique(as.vector(rt$id))
entrezIDs <- mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
gene_mapping <- data.frame(symbol = genes, entrez = entrezIDs, stringsAsFactors = FALSE)
gene_mapping <- gene_mapping[!is.na(gene_mapping$entrez), ]
rt_mapped <- merge(rt, gene_mapping, by.x = "id", by.y = "symbol")
gene <- unique(gene_mapping$entrez)
gene <- gene[grepl("^[0-9]+$", gene)]
cat("有效Entrez ID数量:", length(gene), "\n")

# =============================================================================
# 桑基气泡图函数
# =============================================================================
create_sankey_bubble <- function(enrich_df, rt_data, output_dir, analysis_name, top_n = top_pathways) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    cat(analysis_name, "无结果\n")
    return(NULL)
  }
  write.table(enrich_df, file.path(base_dir, paste0(analysis_name, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  enrich_top <- enrich_df %>%
    arrange(pvalue) %>%
    slice_head(n = min(top_n, nrow(enrich_df)))

  enrichData <- enrich_top %>%
    separate_rows(geneID, sep = "/") %>%
    mutate(
      GeneRatio_numeric = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"),
                                            function(x) as.numeric(x[1]) / as.numeric(x[2]))),
      Count = as.numeric(Count),
      Description = factor(Description, levels = unique(Description))
    )

  # ---- 过滤基因:只保留在 >=2 条通路中出现的 hub 基因,上限 50 ----
  gene_pathway_count <- enrichData %>%
    distinct(geneID, Description) %>%
    count(geneID, sort = TRUE) %>%
    filter(n >= 2)
  if (nrow(gene_pathway_count) > 50) {
    gene_pathway_count <- gene_pathway_count %>% slice_head(n = 50)
  }
  enrichData <- enrichData %>%
    filter(geneID %in% gene_pathway_count$geneID)
  cat(analysis_name, "过滤后基因数:", length(unique(enrichData$geneID)), "\n")

  dfForLodes <- enrichData %>%
    select(geneID, Description, Count) %>%
    rename(Gene = geneID, Pathway = Description, Freq = Count)

  sankeyData <- to_lodes_form(dfForLodes, key = "axis", axes = c("Gene", "Pathway")) %>%
    mutate(y_pos = ifelse(axis == "Gene", gene_y_pos, pathway_y_pos))

  insert_spacer_nodes <- function(data) {
    create_spacers <- function(template_df, spacer_names) {
      if (!length(spacer_names)) return(NULL)
      purrr::map_dfr(spacer_names, ~ {
        row <- template_df[1, , drop = FALSE]
        row$alluvium <- paste0("spacer_", .x)
        row$stratum <- .x
        row
      })
    }
    create_interspersed_nodes <- function(nodes, prefix) {
      if (length(nodes) <= 2) return(nodes)
      spacer_names <- paste0("SPACER_", prefix, "_", seq_len(length(nodes) - 1))
      result <- character(length(nodes) + length(spacer_names))
      result[seq(1, length(result), by = 2)] <- as.character(nodes)
      result[seq(2, length(result), by = 2)] <- spacer_names
      result
    }
    gene_data    <- filter(data, axis == "Gene")
    pathway_data <- filter(data, axis == "Pathway")
    new_gene_order    <- create_interspersed_nodes(unique(gene_data$stratum), "GENE")
    new_pathway_order <- create_interspersed_nodes(unique(pathway_data$stratum), "PATHWAY")
    gene_spacers    <- create_spacers(gene_data, setdiff(new_gene_order, unique(gene_data$stratum)))
    pathway_spacers <- create_spacers(pathway_data, setdiff(new_pathway_order, unique(pathway_data$stratum)))
    combined_data <- rbind(data, gene_spacers, pathway_spacers)
    combined_data$stratum <- factor(combined_data$stratum,
                                    levels = c(new_gene_order, new_pathway_order))
    combined_data
  }

  sankeyData <- insert_spacer_nodes(sankeyData)

  genes_in_data    <- unique(dfForLodes$Gene)
  pathways_in_data <- unique(dfForLodes$Pathway)
  spacer_strata    <- grep("SPACER_", unique(sankeyData$stratum), value = TRUE)

  gene_colors <- setNames(
    colorRampPalette(brewer.pal(min(12, brewer.pal.info[gene_palette_name, "maxcolors"]),
                                gene_palette_name))(length(genes_in_data)),
    genes_in_data)

  pathway_palette <- if (length(pathways_in_data) <= length(pathway_colors)) {
    setNames(pathway_colors[seq_along(pathways_in_data)], pathways_in_data)
  } else {
    setNames(colorRampPalette(c(pathway_colors, brewer.pal(9, "Set1")))(length(pathways_in_data)),
             pathways_in_data)
  }

  nodeColors <- c(gene_colors, pathway_palette,
                  setNames(rep("transparent", length(spacer_strata)), spacer_strata))

  sankeyData <- sankeyData %>%
    mutate(node_color = nodeColors[as.character(stratum)]) %>%
    group_by(alluvium) %>%
    mutate(to_node_name = ifelse(axis == "Pathway", as.character(stratum), NA)) %>%
    fill(to_node_name, .direction = "downup") %>%
    ungroup() %>%
    mutate(flow_color = nodeColors[to_node_name])

  base_theme <- theme(
    text = element_text(family = plot_font, face = "bold"),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

  sankeyPlot <- ggplot(sankeyData,
    aes(x = axis, stratum = stratum, alluvium = alluvium, y = y_pos, label = stratum)) +
    geom_flow(aes(fill = flow_color), alpha = 0.3, width = 0.05, knot.pos = 0.3,
              color = "transparent") +
    geom_stratum(aes(fill = node_color), color = NA, width = 0.05) +
    geom_text(stat = "stratum",
      aes(label = ifelse(grepl("SPACER_", as.character(after_stat(stratum))),
                         "", as.character(after_stat(stratum)))),
      size = sankey_text_size, hjust = 1, nudge_x = -0.04,
      family = plot_font, fontface = "bold") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0.3, 0, 0, 0)) +
    scale_fill_identity() +
    guides(fill = "none") +
    theme_void() + base_theme +
    labs(x = paste0("Gene-Term relationship (", analysis_name, ")")) +
    theme(axis.title.x = element_text(margin = margin(t = 8), size = 12))

  sankeyPlotData <- ggplot_build(sankeyPlot)
  rightNodes <- sankeyPlotData$data[[2]] %>%
    filter(x == max(x)) %>%
    arrange(ymin) %>%
    mutate(
      node_name = as.character(stratum),
      node_center_y = (ymin + ymax) / 2,
      node_height = ymax - ymin
    ) %>%
    filter(node_name %in% pathways_in_data & !grepl("SPACER_", node_name)) %>%
    select(node_name, node_center_y, node_height, ymin, ymax)

  bubbleData <- enrichData %>%
    distinct(Description, .keep_all = TRUE) %>%
    mutate(GeneRatio = GeneRatio_numeric) %>%
    left_join(rightNodes, by = c("Description" = "node_name"))

  bubblePlot <- ggplot(bubbleData,
    aes(x = GeneRatio, y = node_center_y, color = -log10(pvalue))) +
    geom_point(aes(size = Count), stroke = 0.5) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_gradient(low = bubble_low_color, high = bubble_high_color) +
    scale_radius(range = bubble_size_range, name = "Count", breaks = bubble_legend_size) +
    guides(
      color = guide_colorbar(order = 1, barwidth = 0.8, barheight = 4),
      size = guide_legend(order = 2, keywidth = 0.8, keyheight = 0.8, breaks = bubble_legend_size)
    ) +
    labs(size = "Count", color = paste0("-log10(", bubble_p_label, ")"),
         x = bubble_x_label) +
    theme_bw() + base_theme +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(), panel.grid = element_blank(),
      axis.text.x = element_text(margin = margin(t = 4), size = 12),
      axis.title.x = element_text(margin = margin(t = 8), size = 12),
      legend.title = element_text(size = 12), legend.text = element_text(size = 12),
      legend.position = "right", legend.justification = c(0, 0.15),
      legend.margin = margin(0, 0, 0, 10)
    )

  yRange <- sankeyPlotData$layout$panel_params[[1]]$y.range
  combinedPlot <- (sankeyPlot + coord_cartesian(clip = "off", ylim = yRange)) +
    (bubblePlot +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = max(rightNodes$ymax),
               fill = NA, color = "black", linewidth = 0.3) +
      coord_cartesian(ylim = yRange)) +
    plot_layout(widths = c(sankey_width_ratio, dot_width_ratio))

  out_pdf <- file.path(output_dir, paste0(analysis_name, "_sankey_bubble.pdf"))
  out_png <- file.path(output_dir, paste0(analysis_name, "_sankey_bubble.png"))
  out_svg <- file.path(output_dir, paste0(analysis_name, "_sankey_bubble.svg"))
  ggsave(out_pdf, combinedPlot, width = plot_width, height = plot_height + 5, device = cairo_pdf)
  ggsave(out_png, combinedPlot, width = plot_width, height = plot_height + 5, dpi = 300,
         bg = "white")
  ggsave(out_svg, combinedPlot, width = plot_width, height = plot_height + 5, device = svglite::svglite)
  cat(analysis_name, "桑基气泡图已保存:", out_pdf, "\n")
  return(combinedPlot)
}

# =============================================================================
# 富集分析执行
# =============================================================================

# KEGG
cat("\n======== KEGG富集分析 ========\n")
tryCatch({
  kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)
  KEGG <- as.data.frame(kk)
  KEGG$geneID <- as.character(sapply(KEGG$geneID, function(x) {
    ids <- strsplit(x, "/")[[1]]
    symbols <- gene_mapping$symbol[match(ids, gene_mapping$entrez)]
    symbols[is.na(symbols)] <- ids[is.na(symbols)]
    paste(symbols, collapse = "/")
  }))
  KEGG <- KEGG[KEGG$pvalue < pvalueFilter & KEGG$p.adjust < adjPvalFilter, ]
  cat("KEGG显著通路:", nrow(KEGG), "个\n")
  create_sankey_bubble(KEGG, rt_mapped, file.path(base_dir, "KEGG"), "KEGG")
}, error = function(e) cat("KEGG错误:", e$message, "\n"))

# GO BP
cat("\n======== GO BP富集分析 ========\n")
tryCatch({
  go_bp <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                    readable = TRUE)
  go_bp_df <- as.data.frame(go_bp)
  go_bp_df <- go_bp_df[go_bp_df$pvalue < pvalueFilter & go_bp_df$p.adjust < adjPvalFilter, ]
  cat("GO BP显著条目:", nrow(go_bp_df), "个\n")
  create_sankey_bubble(go_bp_df, rt_mapped, file.path(base_dir, "GO_BP"), "GO_BP")
}, error = function(e) cat("GO BP错误:", e$message, "\n"))

# GO CC
cat("\n======== GO CC富集分析 ========\n")
tryCatch({
  go_cc <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, ont = "CC",
                    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                    readable = TRUE)
  go_cc_df <- as.data.frame(go_cc)
  go_cc_df <- go_cc_df[go_cc_df$pvalue < pvalueFilter & go_cc_df$p.adjust < adjPvalFilter, ]
  cat("GO CC显著条目:", nrow(go_cc_df), "个\n")
  create_sankey_bubble(go_cc_df, rt_mapped, file.path(base_dir, "GO_CC"), "GO_CC")
}, error = function(e) cat("GO CC错误:", e$message, "\n"))

# GO MF
cat("\n======== GO MF富集分析 ========\n")
tryCatch({
  go_mf <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, ont = "MF",
                    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                    readable = TRUE)
  go_mf_df <- as.data.frame(go_mf)
  go_mf_df <- go_mf_df[go_mf_df$pvalue < pvalueFilter & go_mf_df$p.adjust < adjPvalFilter, ]
  cat("GO MF显著条目:", nrow(go_mf_df), "个\n")
  create_sankey_bubble(go_mf_df, rt_mapped, file.path(base_dir, "GO_MF"), "GO_MF")
}, error = function(e) cat("GO MF错误:", e$message, "\n"))

# DO
cat("\n======== DO富集分析 ========\n")
tryCatch({
  do_enrich <- enrichDO(gene = gene, pAdjustMethod = "BH",
                        pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  do_df <- as.data.frame(do_enrich)
  do_df <- do_df[do_df$pvalue < pvalueFilter & do_df$p.adjust < adjPvalFilter, ]
  cat("DO显著条目:", nrow(do_df), "个\n")
  create_sankey_bubble(do_df, rt_mapped, file.path(base_dir, "DO"), "DO")
}, error = function(e) cat("DO错误:", e$message, "\n"))

# Reactome - 使用在线 GMT 文件(无需 reactome.db 454MB 包)
cat("\n======== Reactome富集分析 ========\n")
tryCatch({
  reactome_gmt_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c2.cp.reactome.v2024.1.Hs.symbols.gmt"
  reactome_gmt_file <- file.path(base_dir, "c2.cp.reactome.v2024.1.Hs.symbols.gmt")

  # Convert gene symbols to Entrez for enricher
  gene_entrez <- bitr(rt$id, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)$ENTREZID
  cat("Mapped Entrez:", length(gene_entrez), "\n")

  # Convert GMT to Entrez-based via org.Hs.eg.db
  reactome_genesets_raw <- read.gmt(reactome_gmt_file)
  sym2entrez <- bitr(unique(reactome_genesets_raw$gene), fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  reactome_genesets <- reactome_genesets_raw %>%
    left_join(sym2entrez, by = c("gene" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID)) %>%
    select(term, gene = ENTREZID)
  cat("Reactome gene sets loaded:", length(unique(reactome_genesets$term)), "\n")

  reactome_enrich <- enricher(gene = gene_entrez, pvalueCutoff = 1, qvalueCutoff = 1,
                              pAdjustMethod = "BH", TERM2GENE = reactome_genesets)

  reactome_df <- as.data.frame(reactome_enrich)
  reactome_df <- reactome_df[reactome_df$pvalue < pvalueFilter &
                               reactome_df$p.adjust < adjPvalFilter, ]

  # Strip REACTOME_ prefix & convert underscores -> spaces for readability
  reactome_df$Description <- gsub("^REACTOME_", "", reactome_df$Description)
  reactome_df$Description <- gsub("_", " ", reactome_df$Description)

  # Convert geneID (Entrez) -> symbols for display consistency
  reactome_df$geneID <- as.character(sapply(reactome_df$geneID, function(x) {
    ids <- strsplit(x, "/")[[1]]
    symbols <- gene_mapping$symbol[match(ids, gene_mapping$entrez)]
    symbols[is.na(symbols)] <- ids[is.na(symbols)]
    paste(symbols, collapse = "/")
  }))

  cat("Reactome显著通路:", nrow(reactome_df), "个\n")
  create_sankey_bubble(reactome_df, rt_mapped, file.path(base_dir, "REACTOME"),
                       "Reactome")
}, error = function(e) cat("Reactome错误:", e$message, "\n"))

# =============================================================================
# 自动拼图 2×3
# =============================================================================
cat("\n======== 正在生成2x3组合图 ========\n")

# 使用 PNG 拼图(magick 原生支持,无需 pdftools)
png_files <- c(
  file.path(base_dir, "KEGG",   "KEGG_sankey_bubble.png"),
  file.path(base_dir, "GO_BP",  "GO_BP_sankey_bubble.png"),
  file.path(base_dir, "GO_CC",  "GO_CC_sankey_bubble.png"),
  file.path(base_dir, "GO_MF",  "GO_MF_sankey_bubble.png"),
  file.path(base_dir, "DO",     "DO_sankey_bubble.png"),
  file.path(base_dir, "REACTOME","Reactome_sankey_bubble.png")
)
labels <- c("A", "B", "C", "D", "E", "F")
existing_idx <- which(file.exists(png_files))
if (length(existing_idx) >= 2) {
  images <- list()
  for (i in existing_idx) {
    img <- image_read(png_files[i])
    img <- image_border(img, "white", "40x30")
    img <- image_annotate(img, labels[i], size = 42, color = "black",
                          weight = 700, font = "Arial", location = "+8+3")
    images[[as.character(i)]] <- img
  }
  first_img <- images[[1]]
  target_width  <- image_info(first_img)$width
  target_height <- image_info(first_img)$height
  images <- lapply(images, function(img) {
    image_resize(img, paste0(target_width, "x", target_height, "!"))
  })

  all_images <- list()
  for (i in 1:6) {
    if (as.character(i) %in% names(images)) {
      all_images[[i]] <- images[[as.character(i)]]
    } else {
      blank <- image_blank(target_width, target_height, color = "white")
      blank <- image_annotate(blank, paste0(labels[i], " (No data)"),
                              size = 36, color = "grey60", weight = 700,
                              font = "Arial", gravity = "center")
      all_images[[i]] <- blank
    }
  }
  row1 <- image_append(c(all_images[[1]], all_images[[2]], all_images[[3]]))
  row2 <- image_append(c(all_images[[4]], all_images[[5]], all_images[[6]]))
  combined <- image_append(c(row1, row2), stack = TRUE)

  out_pdf  <- file.path(base_dir, "Combined_sankey_bubble_2x3.pdf")
  out_tiff <- file.path(base_dir, "Combined_sankey_bubble_2x3.tiff")
  image_write(combined, out_pdf,  format = "pdf",  density = 300)
  image_write(combined, out_tiff, format = "tiff", density = 300, compression = "lzw")
  cat("组合图已保存:\n", out_pdf, "\n", out_tiff, "\n")
} else {
  cat("有效图不足,无法拼图\n")
}
cat("\n======== 全部分析完成 ========\n")
