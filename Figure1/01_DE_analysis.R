# ============================================================
# GSE43292 动脉粥样硬化火山图 — 完全适配版
# 基于"代码调试.R"框架，所有 KRAS 替换为 GSE43292 数据
# ============================================================

# ---------- 可调参数 ----------
fc_threshold  <- 0.5            # logFC 阈值
fdr_threshold <- 0.05           # FDR 阈值
top_n         <- 200            # 每侧 Top 基因数量
Mark <- c("FABP4", "CD36", "MMP9", "ACTA2", "MYH11", "CNN1")


# ---------- 1. 包 ----------
packages <- c("ggplot2", "dplyr", "ggrepel", "svglite")
need_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(need_install) > 0) {
  install.packages(need_install)
}
library(ggplot2)
library(dplyr)
library(ggrepel)
library(svglite)
library(grid)

# ---------- 2. 工作目录与文件名 ----------
run_args <- commandArgs(FALSE)
file_arg <- run_args[grepl("^--file=", run_args)]
if (length(file_arg) > 0) {
  work_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  work_dir <- normalizePath(".")
}
setwd(work_dir)

csv_file    <- "GSE43292_DE_results.csv"
output_file <- "GSE43292_Volcano_Volcano.svg"   # 与原命名一致
output_png  <- "GSE43292_Volcano_Volcano.png"

if (!file.exists(csv_file)) {
  stop("没有找到文件：", file.path(work_dir, csv_file))
}

# ---------- 3. 读取差异表达数据 ----------
data <- read.csv(
  csv_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 修复空列名
data_names <- colnames(data)
blank_name <- is.na(data_names) | trimws(data_names) == ""
if (any(blank_name)) {
  blank_index <- which(blank_name)
  for (i in seq_along(blank_index)) {
    current_index <- blank_index[i]
    if (current_index == 1) {
      data_names[current_index] <- "ensembl_gene_id"
    } else {
      data_names[current_index] <- paste0("unnamed_column_", i)
    }
  }
}
colnames(data) <- make.unique(data_names)

required_columns <- c("external_gene_name", "logFC", "FDR")
missing_columns  <- setdiff(required_columns, colnames(data))
if (length(missing_columns) > 0) {
  stop("CSV文件缺少必要列：", paste(missing_columns, collapse = ", "))
}

data$external_gene_name <- trimws(as.character(data$external_gene_name))
data$logFC              <- suppressWarnings(as.numeric(data$logFC))
data$FDR                <- suppressWarnings(as.numeric(data$FDR))

data <- data %>%
  filter(
    !is.na(external_gene_name),
    external_gene_name != "",
    !is.na(logFC),
    !is.na(FDR)
  ) %>%
  mutate(
    FDR_plot    = pmax(FDR, 1e-300),
    negLog10FDR = -log10(FDR_plot)
  )

# ---------- 4. 分组：右侧=斑块上调，左侧=斑块下调 ----------
# 右侧：logFC > 0.5  & FDR < 0.05  → Atheroma_UP   (Plaque > Intact)
# 左侧：logFC < -0.5 & FDR < 0.05  → Atheroma_DN   (Intact > Plaque)
Atheroma_UP   <- subset(data, logFC >  fc_threshold & FDR < fdr_threshold)
Atheroma_DN   <- subset(data, logFC < -fc_threshold & FDR < fdr_threshold)

# ---------- 5. Top 基因（按 FDR 排序） ----------
Top_UP <- Atheroma_UP %>%
  arrange(FDR, desc(logFC)) %>%
  slice_head(n = top_n)

Top_DN <- Atheroma_DN %>%
  arrange(FDR, logFC) %>%
  slice_head(n = top_n)

if (nrow(Top_UP) == 0) {
  stop("右侧没有满足阈值的基因，无法绘制 Top ", top_n, " 虚线框。")
}
if (nrow(Top_DN) == 0) {
  stop("左侧没有满足阈值的基因，无法绘制 Top ", top_n, " 虚线框。")
}

message("Atheroma_UP 基因数量：", nrow(Atheroma_UP))
message("Atheroma_DN 基因数量：", nrow(Atheroma_DN))
message("左侧 Top 基因数量：", nrow(Top_DN))
message("右侧 Top 基因数量：", nrow(Top_UP))

# ---------- 6. 构造基因集（替代 GMT 文件）----------
# 从 DE 结果构造 Atheroma 特征基因集（Top N by FDR）
Atheroma_Sig_UP <- Atheroma_UP %>%
  arrange(FDR, desc(logFC)) %>%
  filter(!grepl("^[0-9]+$", external_gene_name)) %>%
  slice_head(n = top_n) %>%
  pull(external_gene_name)

Atheroma_Sig_DN <- Atheroma_DN %>%
  arrange(FDR, logFC) %>%
  filter(!grepl("^[0-9]+$", external_gene_name)) %>%
  slice_head(n = top_n) %>%
  pull(external_gene_name)

message("ATHEROMA_UP_SIGNATURE  基因数量：", length(Atheroma_Sig_UP))
message("ATHEROMA_DN_SIGNATURE  基因数量：", length(Atheroma_Sig_DN))

# ---------- 7. 指定需要标记的基因 ----------
# AS/VCMC 经典标志物（定义见文件顶部可调参数）

data <- data %>%
  mutate(group = case_when(
  external_gene_name %in% Mark            ~ "Mark",
  external_gene_name %in% Atheroma_Sig_UP ~ "ATHEROMA_UP_SIGNATURE",
  external_gene_name %in% Atheroma_Sig_DN ~ "ATHEROMA_DN_SIGNATURE",
  TRUE                                    ~ "Other"
  )
)

# ---------- 8. 构建矩形数据 ----------
y_rectangle_top <- -log10(min(c(Atheroma_UP$FDR_plot, Atheroma_DN$FDR_plot), na.rm = TRUE))

rect_up <- data.frame(
  xmin = min(Atheroma_UP$logFC, na.rm = TRUE),
  xmax = max(Atheroma_UP$logFC, na.rm = TRUE),
  ymin = -log10(fdr_threshold),
  ymax = y_rectangle_top
)

rect_dn <- data.frame(
  xmin = min(Atheroma_DN$logFC, na.rm = TRUE),
  xmax = max(Atheroma_DN$logFC, na.rm = TRUE),
  ymin = -log10(fdr_threshold),
  ymax = y_rectangle_top
)

rect_top_up <- data.frame(
  xmin = min(Top_UP$logFC, na.rm = TRUE),
  xmax = max(Top_UP$logFC, na.rm = TRUE),
  ymin = -log10(max(Top_UP$FDR_plot, na.rm = TRUE)),
  ymax = y_rectangle_top
)

rect_top_dn <- data.frame(
  xmin = min(Top_DN$logFC, na.rm = TRUE),
  xmax = max(Top_DN$logFC, na.rm = TRUE),
  ymin = -log10(max(Top_DN$FDR_plot, na.rm = TRUE)),
  ymax = y_rectangle_top
)

# ---------- 9. 配色（NPG 风格）----------
col_up      <- "#E64B35"    # 红 — Atheroma UP
col_dn      <- "#4DBBD5"    # 蓝 — Atheroma DN
col_sig_up  <- "#FAA500"    # 橙 — UP Signature
col_sig_dn  <- "#C4C4C4"    # 灰 — DN Signature
col_mark    <- "#15F6F8"    # 青 — Mark 基因

# ---------- 10. 坐标 ----------
x_lim  <- max(abs(data$logFC), na.rm = TRUE) * 1.05
y_lim  <- y_rectangle_top * 1.08

# ---------- 11. 绘制火山图 ----------
p <- ggplot(data, aes(x = logFC, y = negLog10FDR)) +

  # 背景色块
  geom_rect(
    data = rect_dn,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#E6E7FC"
  ) +
  geom_rect(
    data = rect_up,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#FFE5E6"
  ) +

  # Top N 虚线框
  geom_rect(
    data = rect_top_dn,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "transparent",
    color = "black",
    linewidth = 1,
    linetype = "dotted"
  ) +
  geom_rect(
    data = rect_top_up,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "transparent",
    color = "black",
    linewidth = 1,
    linetype = "dotted"
  ) +

  # 参考线
  geom_hline(
    yintercept = 0,
    linewidth = 0.6,
    color = "grey60"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.6,
    color = "grey60"
  ) +
  geom_hline(
    yintercept = -log10(fdr_threshold),
    linewidth = 1,
    linetype = "dotted",
    color = "#AD6D6F"
  ) +

  # 散点 — Other
  geom_point(
    data = subset(data, group == "Other"),
    shape = 21,
    color = "black",
    fill  = "white",
    alpha = 0.1,
    size  = 2,
    stroke = 1.2
  ) +

  # 散点 — Atheroma UP Signature
  geom_point(
    data = subset(data, group == "ATHEROMA_UP_SIGNATURE"),
    shape = 21,
    fill  = col_sig_up,
    color = "black",
    size  = 2.8,
    stroke = 1.2
  ) +

  # 散点 — Atheroma DN Signature
  geom_point(
    data = subset(data, group == "ATHEROMA_DN_SIGNATURE"),
    shape = 21,
    fill  = col_sig_dn,
    color = "black",
    size  = 2.8,
    stroke = 1.2
  ) +

  # 散点 — Mark
  geom_point(
    data = subset(data, group == "Mark"),
    shape = 21,
    fill  = col_mark,
    color = "black",
    size  = 2.8,
    stroke = 1.2
  ) +

  # 底部条形码：覆盖文字标签位置
  geom_segment(
    data = subset(data, group == "ATHEROMA_DN_SIGNATURE"),
    aes(x = logFC, y = -1.0, xend = logFC, yend = -0.45),
    color = col_sig_dn,
    linewidth = 0.35
  ) +
  geom_segment(
    data = subset(data, group == "ATHEROMA_UP_SIGNATURE"),
    aes(x = logFC, y = -1.7, xend = logFC, yend = -1.15),
    color = col_sig_up,
    linewidth = 0.35
  ) +

  # Mark 基因标签
  geom_label_repel(
    data = subset(data, group == "Mark"),
    aes(label = external_gene_name),
    fill = "white",
    alpha = 0.8,
    min.segment.length = unit(0.1, "lines"),
    box.padding       = unit(0.35, "lines"),
    point.padding     = unit(0.5, "lines"),
    segment.color     = "grey30",
    direction         = "both",
    seed              = 1,
    max.overlaps      = Inf
  ) +

  # ---- 标注 ----

  # 左右标题
  annotate(
    "text",
    x    = -1.5,
    y    = y_rectangle_top * 0.92,
    label = "Atheroma DOWN\n(Intact > Plaque)",
    lineheight = 0.8,
    size   = 5,
    color  = col_dn,
    fontface = "bold",
    vjust    = 0
  ) +
  annotate(
    "text",
    x    = 1.5,
    y    = y_rectangle_top * 0.92,
    label = "Atheroma UP\n(Plaque > Intact)",
    lineheight = 0.8,
    size   = 5,
    color  = col_up,
    fontface = "bold",
    vjust    = 0
  ) +

  # Top N 标注
  annotate(
    "text",
    x    = min(Top_DN$logFC, na.rm = TRUE) + 0.05,
    y    = y_rectangle_top * 0.82,
    label = paste0("Top ", top_n),
    lineheight = 1,
    size   = 5,
    fontface = "bold",
    hjust    = 0
  ) +
  annotate(
    "text",
    x    = max(Top_UP$logFC, na.rm = TRUE) - 0.05,
    y    = y_rectangle_top * 0.82,
    label = paste0("Top ", top_n),
    lineheight = 1,
    size   = 5,
    fontface = "bold",
    hjust    = 1
  ) +

  # 底部图例（上移，避开条形码）
  annotate(
    "text",
    x    = -x_lim * 1.0,
    y    = -0.45,
    label = "ATHEROMA_DN_SIGNATURE",
    lineheight = 1,
    size   = 5,
    color  = col_sig_dn,
    hjust    = 0
  ) +
  annotate(
    "text",
    x    = -x_lim * 1.0,
    y    = -1.15,
    label = "ATHEROMA_UP_SIGNATURE",
    lineheight = 1,
    size   = 5,
    color  = col_sig_up,
    hjust    = 0
  ) +

  # 右侧基因数量
  annotate(
    "text",
    x    = x_lim * 0.85,
    y    = -0.45,
    label = paste0("n = ", length(Atheroma_Sig_DN)),
    lineheight = 0.8,
    size   = 5,
    color  = col_sig_dn,
    hjust    = 1
  ) +
  annotate(
    "text",
    x    = x_lim * 0.85,
    y    = -1.15,
    label = paste0("n = ", length(Atheroma_Sig_UP)),
    lineheight = 0.8,
    size   = 5,
    color  = col_sig_up,
    hjust    = 1
  ) +

  # FDR 阈值
  annotate(
    "text",
    x    = x_lim * 0.92,
    y    = -log10(fdr_threshold) - 0.05,
    label = as.expression(bquote(alpha == .(fdr_threshold))),
    lineheight = 1,
    size   = 5,
    color  = "#AD6D6F",
    hjust    = 1,
    vjust    = 1
  ) +

  # 数据源
  annotate(
    "text",
    x    = -x_lim * 0.90,
    y    = -2.15,
    label = "GSE43292  |  Carotid Atheroma  |  32 paired  |  RMA",
    lineheight = 1,
    size   = 3.8,
    color  = "grey50",
    fontface = "italic",
    hjust    = 0
  ) +

  scale_x_continuous(
    breaks = pretty(c(-x_lim, x_lim), n = 6),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = pretty(c(0, y_rectangle_top), n = 6),
    expand = expansion(add = c(0, 0.35))
  ) +
  coord_cartesian(
    xlim  = c(-x_lim, x_lim),
    ylim  = c(-2.3, y_lim),
    clip  = "off"
  ) +
  labs(
    x = expression(log[2] ~ "(Atheroma / Intact)"),
    y = "Significance (-log10 FDR)"
  ) +
  theme_classic() +
  theme(
    axis.text  = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(10, 20, 10, 10)
  )

if (interactive()) print(p)

# ---------- 12. 保存 ----------
ggsave(
  filename = output_file,
  plot     = p,
  height   = 6.5,
  width    = 8,
  device   = svglite::svglite,
  bg       = "white"
)
ggsave(
  filename = output_png,
  plot     = p,
  height   = 6.5,
  width    = 8,
  dpi      = 300,
  bg       = "white"
)

message("SVG图片已保存：", file.path(work_dir, output_file))
message("PNG图片已保存：", file.path(work_dir, output_png))
