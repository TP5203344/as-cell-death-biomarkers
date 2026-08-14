# ============================================================
# 仅重绘火山图（使用已有 DE 结果 CSV）
# ============================================================

packages <- c("ggplot2", "dplyr", "ggrepel", "svglite")
need_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(need_install) > 0) install.packages(need_install)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(svglite)

# ---------- 路径 ----------
run_args <- commandArgs(FALSE)
file_arg <- run_args[grepl("^--file=", run_args)]
if (length(file_arg) > 0) {
  out_dir <- normalizePath(dirname(sub("^--file=", "", file_arg[1])))
} else {
  out_dir <- normalizePath(".")
}
csv_file <- file.path(out_dir, "GSE100927_DE_results.csv")
svg_file <- file.path(out_dir, "GSE100927_Volcano.svg")
png_file <- file.path(out_dir, "GSE100927_Volcano.png")

# ---------- 参数 ----------
fc_threshold  <- 0.5
fdr_threshold <- 0.05
top_n         <- 200
Mark <- c("FABP4", "CD36", "MMP9", "ACTA2", "MYH11", "CNN1")

# ---------- 读取 ----------
data <- read.csv(csv_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

data_names <- colnames(data)
blank_name <- is.na(data_names) | trimws(data_names) == ""
if (any(blank_name)) {
  blank_index <- which(blank_name)
  for (k in seq_along(blank_index)) {
    if (blank_index[k] == 1) data_names[blank_index[k]] <- "ensembl_gene_id"
    else data_names[blank_index[k]] <- paste0("unnamed_column_", k)
  }
}
colnames(data) <- make.unique(data_names)

required_columns <- c("external_gene_name", "logFC", "FDR")
missing_columns  <- setdiff(required_columns, colnames(data))
if (length(missing_columns) > 0) stop("CSV缺少必要列：", paste(missing_columns, collapse = ", "))

data$external_gene_name <- trimws(as.character(data$external_gene_name))
data$logFC              <- suppressWarnings(as.numeric(data$logFC))
data$FDR                <- suppressWarnings(as.numeric(data$FDR))

data <- data %>%
  filter(!is.na(external_gene_name), external_gene_name != "",
         !is.na(logFC), !is.na(FDR)) %>%
  mutate(FDR_plot = pmax(FDR, 1e-300),
         negLog10FDR = -log10(FDR_plot))

# ---------- 分组 ----------
Atheroma_UP <- subset(data, logFC >  fc_threshold & FDR < fdr_threshold)
Atheroma_DN <- subset(data, logFC < -fc_threshold & FDR < fdr_threshold)

Top_UP <- Atheroma_UP %>% arrange(FDR, desc(logFC)) %>% slice_head(n = top_n)
Top_DN <- Atheroma_DN %>% arrange(FDR, logFC) %>% slice_head(n = top_n)

message("Atheroma_UP: ", nrow(Atheroma_UP))
message("Atheroma_DN: ", nrow(Atheroma_DN))

Atheroma_Sig_UP <- Atheroma_UP %>%
  arrange(FDR, desc(logFC)) %>%
  filter(!grepl("^[0-9]+$", external_gene_name)) %>%
  slice_head(n = top_n) %>% pull(external_gene_name)

Atheroma_Sig_DN <- Atheroma_DN %>%
  arrange(FDR, logFC) %>%
  filter(!grepl("^[0-9]+$", external_gene_name)) %>%
  slice_head(n = top_n) %>% pull(external_gene_name)

data <- data %>%
  mutate(group = case_when(
    external_gene_name %in% Mark                ~ "Mark",
    external_gene_name %in% Atheroma_Sig_UP     ~ "ATHEROMA_UP_SIGNATURE",
    external_gene_name %in% Atheroma_Sig_DN     ~ "ATHEROMA_DN_SIGNATURE",
    TRUE                                        ~ "Other"
  ))

# ---------- 矩形 ----------
y_rectangle_top <- -log10(min(c(Atheroma_UP$FDR_plot, Atheroma_DN$FDR_plot), na.rm = TRUE))

rect_up <- data.frame(
  xmin = min(Atheroma_UP$logFC, na.rm = TRUE),
  xmax = max(Atheroma_UP$logFC, na.rm = TRUE),
  ymin = -log10(fdr_threshold),
  ymax = y_rectangle_top)

rect_dn <- data.frame(
  xmin = min(Atheroma_DN$logFC, na.rm = TRUE),
  xmax = max(Atheroma_DN$logFC, na.rm = TRUE),
  ymin = -log10(fdr_threshold),
  ymax = y_rectangle_top)

rect_top_up <- data.frame(
  xmin = min(Top_UP$logFC, na.rm = TRUE),
  xmax = max(Top_UP$logFC, na.rm = TRUE),
  ymin = -log10(max(Top_UP$FDR_plot, na.rm = TRUE)),
  ymax = y_rectangle_top)

rect_top_dn <- data.frame(
  xmin = min(Top_DN$logFC, na.rm = TRUE),
  xmax = max(Top_DN$logFC, na.rm = TRUE),
  ymin = -log10(max(Top_DN$FDR_plot, na.rm = TRUE)),
  ymax = y_rectangle_top)

col_up      <- "#E64B35"
col_dn      <- "#4DBBD5"
col_sig_up  <- "#FAA500"
col_sig_dn  <- "#C4C4C4"
col_mark    <- "#15F6F8"

x_lim <- max(abs(data$logFC), na.rm = TRUE) * 1.05
y_lim <- y_rectangle_top * 1.08

# ---------- 绘图 ----------
p <- ggplot(data, aes(x = logFC, y = negLog10FDR)) +
  geom_rect(data = rect_dn, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#E6E7FC") +
  geom_rect(data = rect_up, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#FFE5E6") +
  geom_rect(data = rect_top_dn, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "transparent", color = "black",
            linewidth = 1, linetype = "dotted") +
  geom_rect(data = rect_top_up, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "transparent", color = "black",
            linewidth = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linewidth = 0.6, color = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.6, color = "grey60") +
  geom_hline(yintercept = -log10(fdr_threshold), linewidth = 1,
             linetype = "dotted", color = "#AD6D6F") +
  geom_point(data = subset(data, group == "Other"),
             shape = 21, color = "black", fill = "white",
             alpha = 0.1, size = 2, stroke = 1.2) +
  geom_point(data = subset(data, group == "ATHEROMA_UP_SIGNATURE"),
             shape = 21, fill = col_sig_up, color = "black",
             size = 2.8, stroke = 1.2) +
  geom_point(data = subset(data, group == "ATHEROMA_DN_SIGNATURE"),
             shape = 21, fill = col_sig_dn, color = "black",
             size = 2.8, stroke = 1.2) +
  geom_point(data = subset(data, group == "Mark"),
             shape = 21, fill = col_mark, color = "black",
             size = 2.8, stroke = 1.2) +
  geom_segment(data = subset(data, group == "ATHEROMA_DN_SIGNATURE"),
               aes(x = logFC, y = -1.4, xend = logFC, yend = -0.45),
               color = col_sig_dn, linewidth = 0.35) +
  geom_segment(data = subset(data, group == "ATHEROMA_UP_SIGNATURE"),
               aes(x = logFC, y = -2.1, xend = logFC, yend = -1.55),
               color = col_sig_up, linewidth = 0.35) +
  geom_label_repel(
    data = subset(data, group == "Mark"),
    aes(label = external_gene_name), fill = "white", alpha = 0.8,
    min.segment.length = unit(0.1, "lines"),
    box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines"),
    segment.color = "grey30", direction = "both", seed = 1, max.overlaps = Inf) +
  annotate("text", x = -1.5, y = y_rectangle_top * 0.92,
           label = "Atheroma DOWN\n(Control > Atherosclerotic)",
           lineheight = 0.8, size = 5, color = col_dn, fontface = "bold", vjust = 0) +
  annotate("text", x = 1.5, y = y_rectangle_top * 0.92,
           label = "Atheroma UP\n(Atherosclerotic > Control)",
           lineheight = 0.8, size = 5, color = col_up, fontface = "bold", vjust = 0) +
  annotate("text", x = min(Top_DN$logFC, na.rm = TRUE) + 0.05, y = y_rectangle_top * 0.82,
           label = paste0("Top ", top_n), lineheight = 1, size = 5,
           fontface = "bold", hjust = 0) +
  annotate("text", x = max(Top_UP$logFC, na.rm = TRUE) - 0.05, y = y_rectangle_top * 0.82,
           label = paste0("Top ", top_n), lineheight = 1, size = 5,
           fontface = "bold", hjust = 1) +
  annotate("text", x = -x_lim * 1.0, y = -0.45,
           label = "ATHEROMA_DN_SIGNATURE", lineheight = 1,
           size = 5, color = col_sig_dn, hjust = 0) +
  annotate("text", x = -x_lim * 1.0, y = -1.55,
           label = "ATHEROMA_UP_SIGNATURE", lineheight = 1,
           size = 5, color = col_sig_up, hjust = 0) +
  annotate("text", x = x_lim * 0.85, y = -0.45,
           label = paste0("n = ", length(Atheroma_Sig_DN)),
           lineheight = 0.8, size = 5, color = col_sig_dn, hjust = 1) +
  annotate("text", x = x_lim * 0.85, y = -1.55,
           label = paste0("n = ", length(Atheroma_Sig_UP)),
           lineheight = 0.8, size = 5, color = col_sig_up, hjust = 1) +
  annotate("text", x = x_lim * 0.92, y = -log10(fdr_threshold) - 0.05,
           label = as.expression(bquote(alpha == .(fdr_threshold))),
           lineheight = 1, size = 5, color = "#AD6D6F", hjust = 1, vjust = 1) +
  annotate("text", x = -x_lim * 0.90, y = -3.2,
           label = "GSE100927  |  Peripheral Artery  |  104 samples  |  Agilent G3",
           lineheight = 1, size = 3.8, color = "grey50", fontface = "italic", hjust = 0) +
  scale_x_continuous(breaks = pretty(c(-x_lim, x_lim), n = 6), expand = c(0, 0)) +
  scale_y_continuous(breaks = pretty(c(0, y_rectangle_top), n = 6),
                     expand = expansion(add = c(0, 0.35))) +
  coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(-3.5, y_lim), clip = "off") +
  labs(x = expression(log[2] ~ "(Atherosclerotic / Control)"),
       y = "Significance (-log10 FDR)") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16, face = "bold"),
        legend.title = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(10, 20, 10, 10))

ggsave(filename = svg_file, plot = p, height = 6.5, width = 8,
       device = svglite::svglite, bg = "white")
ggsave(filename = png_file, plot = p, height = 6.5, width = 8,
       dpi = 300, bg = "white")

message("SVG saved: ", svg_file)
message("PNG saved: ", png_file)
message("=== Done ===")
