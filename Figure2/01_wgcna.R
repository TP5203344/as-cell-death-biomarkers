# ==================================================================================
# ------------------------ 参数设置 ---------------------------
# 设置工作目录和输入文件
workDir <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 2 WGCNA鉴定斑块演进模块"
expFilePath <- "Sample Type Matrix.csv"

# 分析参数
min_module_size <- 50                    # 最小模块大小 (建议30-100)
merge_threshold <- 0.25                  # 模块合并阈值 (建议0.15-0.30)
soft_threshold_range <- 1:30             # 软阈值范围
cutHeight_sample <- 15000                # 样本聚类切割高度

# 输出参数
figure_width <- 12                       # 图片宽度 (增加宽度)
figure_height <- 8                       # 图片高度
figure_dpi <- 300                        # 图片分辨率
pdf_pointsize <- 12                      # PDF字体大小

# 统计显著性阈值
p_threshold_high <- 0.001                # 高显著性 (***)
p_threshold_medium <- 0.01               # 中等显著性 (**)
p_threshold_low <- 0.05                  # 低显著性 (*)

# 颜色设置
color_palette <- "RdBu"                  # 热图配色方案
use_custom_colors <- TRUE                # 是否使用自定义颜色

# ------------------------ 环境初始化与包加载 ---------------------------
cat("==== WGCNA高级版分析开始 ====\n")
start_time <- Sys.time()

# 保存参数到变量（避免清理环境时丢失）
temp_params <- list(
  workDir = workDir,
  expFilePath = expFilePath,
  min_module_size = min_module_size,
  merge_threshold = merge_threshold,
  soft_threshold_range = soft_threshold_range,
  cutHeight_sample = cutHeight_sample,
  figure_width = figure_width,
  figure_height = figure_height,
  figure_dpi = figure_dpi,
  pdf_pointsize = pdf_pointsize,
  p_threshold_high = p_threshold_high,
  p_threshold_medium = p_threshold_medium,
  p_threshold_low = p_threshold_low,
  color_palette = color_palette,
  use_custom_colors = use_custom_colors,
  start_time = start_time  # 保存开始时间
)

# 清理环境
rm(list = setdiff(ls(), "temp_params"))
gc()

# 变量恢复
list2env(temp_params, envir = .GlobalEnv)
rm(temp_params)

# 加载包
suppressPackageStartupMessages({
  library(WGCNA)
  library(limma)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  library(scales)
  library(viridis)
  library(corrplot)
  library(cowplot)
})

# 显式指定使用stats包的dist函数，避免冲突
dist <- stats::dist

# 创建输出目录结构
output_dirs <- c(
  "01_Quality_Control",          # 质量控制图表
  "02_Network_Analysis",         # 网络分析图表
  "03_Module_Detection",         # 模块检测图表
  "04_Module_Trait_Analysis",    # 模块-性状关联分析
  "05_Gene_Analysis",           # 基因水平分析
  "06_Data_Output",             # 数据输出文件
  "07_Module_Visualization"     # 模块可视化
)

# 确保在工作目录中创建目录
current_wd <- getwd()
setwd(workDir)  # 切换到工作目录

for(dir in output_dirs) {
  if(!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("创建目录:", dir, "\n")
  }
}

cat("当前工作目录:", getwd(), "\n")
cat("输出目录创建完成\n")

# 主题设置
theme_sci <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 1),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold")
    )
}

# 设置随机种子和进度跟踪
set.seed(12345)

# 进度跟踪
total_steps <- 12
current_step <- 0
update_progress <- function(message) {
  current_step <<- current_step + 1
  cat(sprintf("[%d/%d] %s\n", current_step, total_steps, message))
}

# ------------------------ 数据读取与预处理 ---------------------------
update_progress("数据读取与预处理")

# 检查文件是否存在
if (!file.exists(expFilePath)) {
  stop("表达矩阵文件不存在！请检查路径: ", file.path(getwd(), expFilePath))
}

# 自动检测分隔符并读取数据
tmp_head <- readLines(expFilePath, 1)
sep <- ifelse(grepl(",", tmp_head), ",", "\t")
rawData <- read.table(expFilePath, sep = sep, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(rawData) <- rawData[, 1]
exprData <- rawData[, -1]

# 数据预处理流程
dimNames <- list(rownames(exprData), colnames(exprData))
numericData <- matrix(as.numeric(as.matrix(exprData)), nrow = nrow(exprData), dimnames = dimNames)

# Log2转换
if (max(numericData, na.rm = TRUE) > 1000) {
  numericData <- log2(numericData + 1)
  cat("已进行Log2转换\n")
} else {
  cat("数据似乎已经过Log2转换，跳过\n")
}

# limma归一化
normalizedData <- normalizeBetweenArrays(numericData)

# 过滤低方差基因 
gene_sd <- apply(normalizedData, 1, sd, na.rm = TRUE)
filteredData <- normalizedData[gene_sd > 0.25, ]
cat("过滤后保留基因数:", nrow(filteredData), "\n")

# 自动分组识别
sample_names <- colnames(filteredData)
group_info <- ifelse(grepl("_con$", sample_names, ignore.case = TRUE), "Control",
                     ifelse(grepl("_tre$", sample_names, ignore.case = TRUE), "Treatment", "Unknown"))
if (any(group_info == "Unknown")) stop("存在未识别分组的样本名！")
numControl <- sum(group_info == "Control")
numTreatment <- sum(group_info == "Treatment")
cat("对照组样本数:", numControl, "，处理组样本数:", numTreatment, "\n")

# 数据转置为WGCNA格式
exprMatrix <- t(filteredData)

# ------------------------ 样本质量控制 ---------------------------
update_progress("样本质量控制")

# 样本和基因质量检查
enableWGCNAThreads()
sampleCheck <- goodSamplesGenes(exprMatrix, verbose = 3)
if (!sampleCheck$allOK) {
  exprMatrix <- exprMatrix[sampleCheck$goodSamples, sampleCheck$goodGenes]
  cat("移除低质量样本/基因后，维度:", dim(exprMatrix), "\n")
}

# 样本聚类分析
sampleDist <- dist(exprMatrix)
sampleDendro <- hclust(sampleDist, method = "average")

# 高质量样本聚类图（带红蓝分组色条 + 样本名标签）---- PDF + SVG
# 用 WGCNA 内置 plotDendroAndColors，同时画树 + 分组色条
plot_file <- file.path("01_Quality_Control", "01_Sample_Clustering_Analysis")
max_h <- max(sampleDendro$height)
sample_colors <- ifelse(group_info[match(rownames(exprMatrix), sample_names)] == "Control",
                        "#3498DB", "#E74C3C")
color_matrix <- cbind(Control = ifelse(sample_colors == "#3498DB", "#3498DB", "#E74C3C"))
rownames(color_matrix) <- rownames(exprMatrix)
# 缩短标签: 去掉 "GSM1" 前缀，只留后 6 位 + _con/_tre
dendroLabels <- gsub("^GSM1", "", rownames(exprMatrix))

pdf(paste0(plot_file, ".pdf"), width = figure_width, height = figure_height, pointsize = pdf_pointsize)
plotDendroAndColors(sampleDendro, colors = color_matrix,
                    groupLabels = c("Group"),
                    dendroLabels = dendroLabels, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Sample Clustering for Outlier Detection",
                    cex.colorLabels = 1.0, cex.dendroLabels = 0.55,
                    marAll = c(1, 8, 4, 2))
dev.off()

svg(paste0(plot_file, ".svg"), width = figure_width, height = figure_height)
plotDendroAndColors(sampleDendro, colors = color_matrix,
                    groupLabels = c("Group"),
                    dendroLabels = dendroLabels, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Sample Clustering for Outlier Detection",
                    cex.colorLabels = 1.0, cex.dendroLabels = 0.55,
                    marAll = c(1, 8, 4, 2))
dev.off()

# 异常样本检测
cutHeight <- cutHeight_sample
clusterCut <- cutreeStatic(sampleDendro, cutHeight = cutHeight, minSize = 3)
if (length(unique(clusterCut)) > 1) {
  keepSamples <- (clusterCut == 1)
  exprMatrix <- exprMatrix[keepSamples, ]
  cat("移除异常样本后，样本数:", nrow(exprMatrix), "\n")
}

# ------------------------ 临床性状数据准备 ---------------------------
update_progress("临床性状数据准备")

# 创建临床数据
clinicalData <- data.frame(
  Control = ifelse(group_info[match(rownames(exprMatrix), sample_names)] == "Control", 1, 0),
  Disease = ifelse(group_info[match(rownames(exprMatrix), sample_names)] == "Treatment", 1, 0))
rownames(clinicalData) <- rownames(exprMatrix)

# 样本-性状热图 ---- PDF + SVG
plot_file <- file.path("01_Quality_Control", "02_Sample_Trait_Heatmap")
sampleColors <- numbers2colors(clinicalData, signed = FALSE)
pdf(paste0(plot_file, ".pdf"), width = figure_width + 2, height = figure_height)
plotDendroAndColors(sampleDendro, sampleColors,
                    groupLabels = names(clinicalData),
                    main = "Sample Dendrogram and Clinical Traits")
dev.off()
svg(paste0(plot_file, ".svg"), width = figure_width + 2, height = figure_height)
plotDendroAndColors(sampleDendro, sampleColors,
                    groupLabels = names(clinicalData),
                    main = "Sample Dendrogram and Clinical Traits")
dev.off()

# ------------------------ 软阈值选择 ---------------------------
update_progress("软阈值选择与网络拓扑分析")

# 软阈值计算
powerVector <- soft_threshold_range
sftResult <- pickSoftThreshold(exprMatrix, powerVector = powerVector, verbose = 5)
if (is.null(sftResult$powerEstimate)) {
  optimalPower <- 6
  warning("无法自动确定最佳软阈值，使用默认值6")
} else {
  optimalPower <- sftResult$powerEstimate
}
cat("选择的软阈值:", optimalPower, "\n")

sft_data <- data.frame(
  Power = sftResult$fitIndices[, 1],
  SFT_R2 = -sign(sftResult$fitIndices[, 3]) * sftResult$fitIndices[, 2],
  Mean_Connectivity = sftResult$fitIndices[, 5])

p1 <- ggplot(sft_data, aes(x = Power, y = SFT_R2)) +
  geom_point(size = 3, color = "#2C3E50") +
  geom_text(aes(label = Power), hjust = -0.3, vjust = -0.3, size = 3.5) +
  geom_hline(yintercept = 0.9, color = "#E74C3C", linetype = "dashed", size = 1) +
  labs(x = "Soft Threshold (power)", y = "Scale Free Topology Model Fit (signed R²)",
       title = "Scale Independence") + theme_sci() + ylim(c(0, 1))

p2 <- ggplot(sft_data, aes(x = Power, y = Mean_Connectivity)) +
  geom_point(size = 3, color = "#8E44AD") +
  geom_text(aes(label = Power), hjust = -0.3, vjust = -0.3, size = 3.5) +
  labs(x = "Soft Threshold (power)", y = "Mean Connectivity",
       title = "Mean Connectivity") + theme_sci()

# 合并图表 ---- PDF + SVG
plot_file <- file.path("02_Network_Analysis", "01_Soft_Threshold_Analysis")
pdf(paste0(plot_file, ".pdf"), width = figure_width + 2, height = figure_height - 2)
grid.arrange(p1, p2, ncol = 2)
dev.off()
svg(paste0(plot_file, ".svg"), width = figure_width + 2, height = figure_height - 2)
grid.arrange(p1, p2, ncol = 2)
dev.off()

# ------------------------ 网络构建与模块检测 ---------------------------
update_progress("网络构建与TOM计算")

adjacencyMatrix <- adjacency(exprMatrix, power = optimalPower)
TOMMatrix <- TOMsimilarity(adjacencyMatrix)
dissTOM <- 1 - TOMMatrix
geneDendro <- hclust(as.dist(dissTOM), method = "average")

# 基因聚类树图 ---- PDF + SVG
plot_file <- file.path("02_Network_Analysis", "02_Gene_Clustering_Dendrogram")
pdf(paste0(plot_file, ".pdf"), width = figure_width + 5, height = figure_height)
par(cex = 0.8, mar = c(2, 6, 4, 2))
plot(geneDendro, xlab = "", sub = "", main = "Gene Clustering Based on TOM Dissimilarity",
     labels = FALSE, hang = 0.04, cex.main = 1.5, cex.lab = 1.3)
dev.off()
svg(paste0(plot_file, ".svg"), width = figure_width + 5, height = figure_height)
par(cex = 0.8, mar = c(2, 6, 4, 2))
plot(geneDendro, xlab = "", sub = "", main = "Gene Clustering Based on TOM Dissimilarity",
     labels = FALSE, hang = 0.04, cex.main = 1.5, cex.lab = 1.3)
dev.off()

# ------------------------ 动态模块检测 ---------------------------
update_progress("动态模块检测与合并")

minModuleSize <- min_module_size
dynamicMods <- cutreeDynamic(dendro = geneDendro, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
moduleColors <- labels2colors(dynamicMods)
cat("检测到", length(unique(moduleColors)), "个模块\n")

MEs0 <- moduleEigengenes(exprMatrix, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
mergeThreshold <- merge_threshold

# 模块合并分析图 ---- PDF + SVG
plot_file <- file.path("03_Module_Detection", "01_Module_Merging_Analysis")
# 优化: 关闭默认叶标签(避免低高度标签掉到x轴下方), 底部统一加模块名(按模块颜色着色), 阈值注释放右上角
me_order <- METree$order
me_labels <- gsub("^ME", "", METree$labels)
me_colors <- tolower(me_labels)  # 模块名本身就是颜色名
me_colors[me_colors == "turquoise"] <- "#1B9E77"  # turquoise 可读性更好
# 提取树中每个叶子的颜色 (模块名=颜色名)
leaf_colors <- me_colors[me_order]
ylim_top <- max(METree$height) * 1.25

pdf(paste0(plot_file, ".pdf"), width = figure_width, height = figure_height - 2)
par(cex = 0.8, mar = c(6, 6, 4, 2))
plot(METree, main = "Clustering of Module Eigengenes", xlab = "", sub = "",
     labels = FALSE, hang = 0.05, cex.main = 1.5, cex.axis = 1.1,
     ylim = c(0, ylim_top))
abline(h = mergeThreshold, col = "#E74C3C", lwd = 2, lty = 2)
# 底部统一标签（按模块颜色着色）
x_positions <- seq_along(me_order)
mtext(me_labels[me_order], side = 1, line = 1.2, at = x_positions,
      col = leaf_colors, cex = 1.1)
# 阈值注释贴虚线（右上角空白区）
usr <- par("usr")
text(x = usr[2] - (usr[2] - usr[1]) * 0.02, y = mergeThreshold * 1.22,
     labels = paste("Merge threshold =", mergeThreshold), col = "#E74C3C", cex = 1.2, adj = 1)
dev.off()
svg(paste0(plot_file, ".svg"), width = figure_width, height = figure_height - 2)
par(cex = 0.8, mar = c(6, 6, 4, 2))
plot(METree, main = "Clustering of Module Eigengenes", xlab = "", sub = "",
     labels = FALSE, hang = 0.05, cex.main = 1.5, cex.axis = 1.1,
     ylim = c(0, ylim_top))
abline(h = mergeThreshold, col = "#E74C3C", lwd = 2, lty = 2)
x_positions <- seq_along(me_order)
mtext(me_labels[me_order], side = 1, line = 1.2, at = x_positions,
      col = leaf_colors, cex = 1.1)
usr <- par("usr")
text(x = usr[2] - (usr[2] - usr[1]) * 0.02, y = mergeThreshold * 1.22,
     labels = paste("Merge threshold =", mergeThreshold), col = "#E74C3C", cex = 1.2, adj = 1)
dev.off()

merge <- mergeCloseModules(exprMatrix, moduleColors, cutHeight = mergeThreshold, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# 合并前后对比图 ---- PDF + SVG
plot_file <- file.path("03_Module_Detection", "02_Module_Colors_Comparison")
pdf(paste0(plot_file, ".pdf"), width = figure_width + 5, height = figure_height)
plotDendroAndColors(geneDendro, cbind(moduleColors, mergedColors),
                    c("Original", "Merged"), dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Module Assignment: Before and After Merging")
dev.off()
svg(paste0(plot_file, ".svg"), width = figure_width + 5, height = figure_height)
plotDendroAndColors(geneDendro, cbind(moduleColors, mergedColors),
                    c("Original", "Merged"), dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Module Assignment: Before and After Merging")
dev.off()

moduleColors <- mergedColors
MEs <- mergedMEs
cat("合并后模块数:", length(unique(moduleColors)), "\n")

# 重新绘制基因聚类树图（底部附合并后模块颜色条）---- PDF + SVG
# 覆盖第5步的纯树版本，形成“树 + 模块色条”标准 WGCNA 图
plot_file <- file.path("02_Network_Analysis", "02_Gene_Clustering_Dendrogram")
pdf(paste0(plot_file, ".pdf"), width = figure_width + 5, height = figure_height)
plotDendroAndColors(geneDendro, colors = cbind(Module = moduleColors),
                    groupLabels = c("Module Colors"),
                    dendroLabels = FALSE, hang = 0.04,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Clustering Based on TOM Dissimilarity",
                    cex.colorLabels = 1.0, cex.dendroLabels = 0.6,
                    marAll = c(1, 8, 4, 2))
dev.off()
svg(paste0(plot_file, ".svg"), width = figure_width + 5, height = figure_height)
plotDendroAndColors(geneDendro, colors = cbind(Module = moduleColors),
                    groupLabels = c("Module Colors"),
                    dendroLabels = FALSE, hang = 0.04,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Clustering Based on TOM Dissimilarity",
                    cex.colorLabels = 1.0, cex.dendroLabels = 0.6,
                    marAll = c(1, 8, 4, 2))
dev.off()

# ------------------------ 模块-性状关联分析 ---------------------------
update_progress("模块-性状关联分析")

moduleTraitCor <- cor(MEs, clinicalData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(exprMatrix))

# 方法1: 经典WGCNA labeledHeatmap风格 ---- PDF + SVG
plot_file <- file.path("04_Module_Trait_Analysis", "01_Module_Trait_Classic_WGCNA")
pdf(paste0(plot_file, ".pdf"), width = 8, height = 10)
par(mar = c(6, 8.5, 3, 3))
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(clinicalData),
               yLabels = names(MEs), ySymbols = names(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 0.8, cex.lab.x = 1.2, cex.lab.y = 1.0,
               zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()
svg(paste0(plot_file, ".svg"), width = 8, height = 10)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(clinicalData),
               yLabels = names(MEs), ySymbols = names(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 0.8, cex.lab.x = 1.2, cex.lab.y = 1.0,
               zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()

# 方法2: 增强版经典风格（带模块颜色条）---- PDF + SVG
cor_matrix <- as.matrix(moduleTraitCor)
p_matrix <- as.matrix(moduleTraitPvalue)
text_labels <- matrix(paste(sprintf("%.2f", cor_matrix),
                            "\n(", ifelse(p_matrix < 1e-100, "< 1e-100", sprintf("%.0e", p_matrix)), ")",
                            sep = ""), nrow = nrow(cor_matrix))
module_names <- rownames(cor_matrix)
module_colors <- gsub("ME", "", module_names)
color_annotation <- data.frame(Module = factor(module_colors, levels = module_colors))
rownames(color_annotation) <- module_names
module_color_map <- structure(module_colors, names = module_colors)

plot_file <- file.path("04_Module_Trait_Analysis", "02_Module_Trait_Enhanced_Classic")
pdf(paste0(plot_file, ".pdf"), width = 10, height = 12)
pheatmap(cor_matrix,
         color = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                                    "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                                    "#D6604D", "#B2182B", "#67001F"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize = 12, fontsize_row = 12, fontsize_col = 14,
         main = "Module-trait relationships", border_color = "white",
         cellwidth = 80, cellheight = 35,
         display_numbers = text_labels, number_color = "black", fontsize_number = 10,
         annotation_row = color_annotation, annotation_colors = list(Module = module_color_map),
         annotation_names_row = FALSE,
         legend = TRUE, legend_breaks = c(-1, -0.5, 0, 0.5, 1),
         legend_labels = c("-1", "-0.5", "0", "0.5", "1"))
dev.off()
svg(paste0(plot_file, ".svg"), width = 10, height = 12)
pheatmap(cor_matrix,
         color = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                                    "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                                    "#D6604D", "#B2182B", "#67001F"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize = 12, fontsize_row = 12, fontsize_col = 14,
         main = "Module-trait relationships", border_color = "white",
         cellwidth = 80, cellheight = 35,
         display_numbers = text_labels, number_color = "black", fontsize_number = 10,
         annotation_row = color_annotation, annotation_colors = list(Module = module_color_map),
         annotation_names_row = FALSE,
         legend = TRUE, legend_breaks = c(-1, -0.5, 0, 0.5, 1),
         legend_labels = c("-1", "-0.5", "0", "0.5", "1"))
dev.off()

# 方法3: ggplot2版本 ---- PDF + SVG
library(ggplot2)
library(reshape2)
cor_melted <- melt(moduleTraitCor)
p_melted <- melt(moduleTraitPvalue)
combined_data <- merge(cor_melted, p_melted, by = c("Var1", "Var2"))
colnames(combined_data) <- c("Module", "Trait", "Correlation", "Pvalue")
combined_data$ModuleColor <- gsub("ME", "", combined_data$Module)
combined_data$Label <- paste(sprintf("%.2f", combined_data$Correlation),
                             "\n(", ifelse(combined_data$Pvalue < 1e-100, "< 1e-100",
                                           sprintf("%.0e", combined_data$Pvalue)), ")", sep = "")
module_order <- names(MEs)
combined_data$Module <- factor(combined_data$Module, levels = rev(module_order))

p_classic <- ggplot(combined_data, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white", size = 0.8) +
  scale_fill_gradient2(low = "#053061", mid = "white", high = "#67001F",
                       midpoint = 0, limits = c(-1, 1), name = "",
                       breaks = c(-1, -0.5, 0, 0.5, 1),
                       labels = c("-1", "-0.5", "0", "0.5", "1")) +
  geom_text(aes(label = Label), size = 3.2, color = "black", fontface = "bold") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14, color = "black", face = "bold"),
        axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.title = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black", margin = margin(b = 20)),
        panel.background = element_rect(fill = "white", color = NA), panel.grid = element_blank(), panel.border = element_blank(),
        legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 12, face = "bold"),
        legend.key.width = unit(1.2, "cm"), legend.key.height = unit(4, "cm"), legend.margin = margin(l = 20),
        plot.margin = margin(20, 30, 20, 20)) +
  labs(title = "Module-trait relationships")

plot_file <- file.path("04_Module_Trait_Analysis", "03_Module_Trait_ggplot2_Classic")
ggsave(paste0(plot_file, ".pdf"), p_classic, width = 8, height = 10, device = "pdf")
ggsave(paste0(plot_file, ".svg"), p_classic, width = 8, height = 10, device = "svg")

# ------------------------ 基因重要性分析 ---------------------------
update_progress("基因重要性与模块归属度分析")

geneModuleMembership <- as.data.frame(cor(exprMatrix, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(exprMatrix)))
geneTraitSignificance <- as.data.frame(cor(exprMatrix, clinicalData, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(exprMatrix)))

unique_modules <- unique(moduleColors)
target_trait <- "Disease"

for (module in unique_modules) {
  if (module == "grey") next
  moduleGenes <- (moduleColors == module)
  mm_column <- paste0("ME", module)
  if (mm_column %in% colnames(geneModuleMembership)) {
    MM <- abs(geneModuleMembership[moduleGenes, mm_column])
    GS <- abs(geneTraitSignificance[moduleGenes, target_trait])
    cor_test <- cor.test(MM, GS, method = "pearson")
    cor_value <- round(cor_test$estimate, 3)
    p_value <- format(cor_test$p.value, scientific = TRUE, digits = 2)
    plot_data <- data.frame(MM = MM, GS = GS)
    p_scatter <- ggplot(plot_data, aes(x = MM, y = GS)) +
      geom_point(color = module, alpha = 0.7, size = 2) +
      geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
      labs(x = paste("Module Membership in", module, "module"),
           y = paste("Gene Significance for", target_trait),
           title = paste0("Module: ", module, "\nCorrelation = ", cor_value, ", P-value = ", p_value)) +
      theme_sci() +
      annotate("text", x = Inf, y = Inf, label = paste("n =", sum(moduleGenes)), hjust = 1.1, vjust = 1.1, size = 4)
    plot_file <- file.path("05_Gene_Analysis", paste0("MM_vs_GS_", module, "_module"))
    ggsave(paste0(plot_file, ".pdf"), p_scatter, width = 6, height = 6, device = "pdf")
    ggsave(paste0(plot_file, ".svg"), p_scatter, width = 6, height = 6, device = "svg")
  }
}

# ------------------------ 模块可视化 ---------------------------
update_progress("模块可视化与网络图")

module_sizes <- table(moduleColors)
size_data <- data.frame(Module = names(module_sizes), GeneCount = as.numeric(module_sizes))
size_data$Module <- factor(size_data$Module, levels = size_data$Module[order(size_data$GeneCount, decreasing = TRUE)])

p_sizes <- ggplot(size_data, aes(x = Module, y = GeneCount, fill = Module)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  scale_fill_identity() +
  labs(title = "Gene Counts per Module", x = "Module", y = "Number of Genes") +
  theme_sci() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  geom_text(aes(label = GeneCount), vjust = -0.3, size = 3)

plot_file <- file.path("07_Module_Visualization", "01_Module_Gene_Counts")
ggsave(paste0(plot_file, ".pdf"), p_sizes, width = figure_width, height = figure_height - 2, device = "pdf")
ggsave(paste0(plot_file, ".svg"), p_sizes, width = figure_width, height = figure_height - 2, device = "svg")

cat("开始生成模块表达热图...\n")
if (!exists("geneModuleMembership") || !exists("moduleColors") || !exists("group_info")) {
  cat("警告：缺少必要的变量，跳过模块热图生成\n")
} else {
  valid_modules <- size_data$Module[size_data$Module != "grey"]
  top_modules <- head(valid_modules, 6)
  cat("将为以下模块生成热图:", paste(top_modules, collapse = ", "), "\n")
  for (module in top_modules) {
    tryCatch({
      module_genes <- colnames(exprMatrix)[moduleColors == module]
      if (length(module_genes) < 5) { cat("模块", module, "基因数少于5，跳过\n"); next }
      cat("处理模块:", module, "，基因数:", length(module_genes), "\n")
      mm_col <- paste0("ME", module)
      if (mm_col %in% colnames(geneModuleMembership)) {
        module_gene_indices <- which(moduleColors == module)
        gene_mm <- abs(geneModuleMembership[module_gene_indices, mm_col])
        names(gene_mm) <- colnames(exprMatrix)[module_gene_indices]
        top_genes <- names(sort(gene_mm, decreasing = TRUE))[1:min(50, length(gene_mm))]
        module_expr <- t(exprMatrix[, top_genes, drop = FALSE])
        expr_samples <- colnames(module_expr)
        matched_groups <- group_info[match(expr_samples, sample_names)]
        display_groups <- ifelse(matched_groups == "Treatment", "Disease", matched_groups)
        sample_annotation <- data.frame(Group = factor(display_groups, levels = c("Control", "Disease")))
        rownames(sample_annotation) <- expr_samples
        if (any(is.na(matched_groups))) {
          cat("警告：模块", module, "存在未匹配的样本，使用默认分组\n")
          sample_annotation$Group <- factor(rep(c("Control", "Disease"), length.out = ncol(module_expr)))
        }
        ann_colors <- list(Group = c(Control = "#3498DB", Disease = "#E74C3C"))
        plot_file <- file.path("07_Module_Visualization", paste0("Module_", module, "_Expression_Heatmap"))
        pdf(paste0(plot_file, ".pdf"), width = figure_width + 2, height = figure_width)
        pheatmap(module_expr, annotation_col = sample_annotation, annotation_colors = ann_colors,
                 scale = "row", clustering_distance_rows = "correlation", clustering_distance_cols = "euclidean",
                 show_rownames = FALSE, show_colnames = TRUE, fontsize = 10, fontsize_col = 8,
                 color = colorRampPalette(c("#2166AC", "white", "#D73027"))(100),
                 main = paste("Expression Heatmap for", module, "Module"), border_color = NA)
        dev.off()
        svg(paste0(plot_file, ".svg"), width = figure_width + 2, height = figure_width)
        pheatmap(module_expr, annotation_col = sample_annotation, annotation_colors = ann_colors,
                 scale = "row", clustering_distance_rows = "correlation", clustering_distance_cols = "euclidean",
                 show_rownames = FALSE, show_colnames = TRUE, fontsize = 10, fontsize_col = 8,
                 color = colorRampPalette(c("#2166AC", "white", "#D73027"))(100),
                 main = paste("Expression Heatmap for", module, "Module"), border_color = NA)
        dev.off()
        cat("✓ 模块", module, "热图已生成\n")
      } else {
        cat("警告：未找到模块", module, "的特征基因列，跳过\n")
      }
    }, error = function(e) {
      cat("错误：生成模块", module, "热图时出错:", e$message, "\n")
      if (dev.cur() > 1) dev.off()
    })
  }
}

# ------------------------ 数据输出 ---------------------------
update_progress("数据输出与文件整理")

gene_info <- data.frame(Gene = colnames(exprMatrix), Module = moduleColors, stringsAsFactors = FALSE)
for (mod in colnames(MEs)) {
  gene_info[, paste0("MM_", mod)] <- geneModuleMembership[, mod]
  gene_info[, paste0("MM_pvalue_", mod)] <- MMPvalue[, mod]
}
for (trait in colnames(clinicalData)) {
  gene_info[, paste0("GS_", trait)] <- geneTraitSignificance[, trait]
  gene_info[, paste0("GS_pvalue_", trait)] <- GSPvalue[, trait]
}

write.table(gene_info, file = file.path("06_Data_Output", "01_Gene_Module_Information.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(MEs, file = file.path("06_Data_Output", "02_Module_Eigengenes.txt"),
            sep = "\t", row.names = TRUE, quote = FALSE)
write.table(moduleTraitCor, file = file.path("06_Data_Output", "03_Module_Trait_Correlations.txt"),
            sep = "\t", row.names = TRUE, quote = FALSE)
write.table(moduleTraitPvalue, file = file.path("06_Data_Output", "04_Module_Trait_Pvalues.txt"),
            sep = "\t", row.names = TRUE, quote = FALSE)

for (module in unique(moduleColors)) {
  module_genes <- gene_info$Gene[gene_info$Module == module]
  write.table(module_genes, file = file.path("06_Data_Output", paste0("05_Module_", module, "_Genes.txt")),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# ------------------------ 分析摘要报告 ---------------------------
update_progress("生成分析摘要报告")

summary_stats <- list(
  total_genes = ncol(exprMatrix), total_samples = nrow(exprMatrix),
  control_samples = sum(clinicalData$Control), treatment_samples = sum(clinicalData$Disease),
  total_modules = length(unique(moduleColors)), soft_threshold = optimalPower,
  largest_module = names(sort(table(moduleColors), decreasing = TRUE))[1],
  largest_module_size = max(table(moduleColors)))

summary_text <- paste(
  "=== WGCNA分析摘要 ===",
  paste("分析时间:", Sys.time()),
  paste("总基因数:", summary_stats$total_genes),
  paste("总样本数:", summary_stats$total_samples),
  paste("对照组样本数:", summary_stats$control_samples),
  paste("处理组样本数:", summary_stats$treatment_samples),
  paste("检测到的模块数:", summary_stats$total_modules),
  paste("使用的软阈值:", summary_stats$soft_threshold),
  paste("最大模块:", summary_stats$largest_module),
  paste("最大模块基因数:", summary_stats$largest_module_size),
  "", "=== 输出文件说明 ===",
  "01_Quality_Control/: 质量控制图表",
  "02_Network_Analysis/: 网络分析图表",
  "03_Module_Detection/: 模块检测图表",
  "04_Module_Trait_Analysis/: 模块-性状关联分析",
  "05_Gene_Analysis/: 基因水平分析",
  "06_Data_Output/: 数据输出文件",
  "07_Module_Visualization/: 模块可视化", "",
  sep = "\n")
writeLines(summary_text, file.path("06_Data_Output", "00_Analysis_Summary.txt"))

# ------------------------ 完成 ---------------------------
end_time <- Sys.time()
total_time <- round(difftime(end_time, start_time, units = "mins"), 2)
cat("\n==== WGCNA分析完成 ====\n")
cat("总用时:", total_time, "分钟\n")
cat("所有结果已保存到相应文件夹中\n")
cat("请查看 06_Data_Output/00_Analysis_Summary.txt 获取详细摘要\n")

rm(list = setdiff(ls(), c("exprMatrix", "MEs", "moduleColors", "clinicalData")))
gc()