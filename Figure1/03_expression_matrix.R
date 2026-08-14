# Step 1: 读取 raw → 表达矩阵 → 保存 RDS
library(data.table)
library(limma)

raw_dir  <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 1 研究设计与转录组差异分析/GSE100927_RAW"
rds_file <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 1 研究设计与转录组差异分析/热图GSE100927/expr_norm.rds"

gz_files <- list.files(raw_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
cat("Found", length(gz_files), "raw files\n")

read_agilent_one <- function(fname) {
  con <- gzfile(fname, "rt")
  header_line <- 0; line_num <- 0
  while (TRUE) {
    l <- readLines(con, n = 1, warn = FALSE)
    if (length(l) == 0) break
    line_num <- line_num + 1
    if (grepl("^FEATURES", l)) { header_line <- line_num; break }
  }
  close(con)
  if (header_line == 0) return(NULL)
  dat <- tryCatch(
    fread(fname, skip = header_line - 1, header = TRUE, sep = "\t",
          quote = "", fill = TRUE, showProgress = FALSE, nThread = 1,
          select = c("ProbeName","GeneName","gProcessedSignal")),
    error = function(e) NULL
  )
  if (is.null(dat)) return(NULL)
  dat$gProcessedSignal <- suppressWarnings(as.numeric(dat$gProcessedSignal))
  dat <- dat[!is.na(dat$gProcessedSignal), ]
  dat <- dat[!grepl("^(GE_BrightCorner|DarkCorner|ETG|E1A|r60|rRNA|Ladder|Control|.*_control_|.*_pos_|.*_neg_)",
                    dat$ProbeName, ignore.case = TRUE), ]
  dat
}

expr_list <- list()
sample_names <- character(0)

for (i in seq_along(gz_files)) {
  f <- gz_files[i]
  gsm <- sub("^(GSM\\d+)_.*", "\\1", basename(f))
  if (i %% 5 == 0) cat(sprintf("  [%d/%d] %s\n", i, length(gz_files), gsm))
  dat <- read_agilent_one(f)
  if (is.null(dat) || nrow(dat) == 0) { cat("  SKIP:", gsm, "\n"); next }
  expr_list[[gsm]] <- dat
  sample_names <- c(sample_names, gsm)
}

cat("Merging matrix...\n")
all_probes <- unique(unlist(lapply(expr_list, function(x) x$ProbeName)))
expr_mat <- matrix(NA, nrow = length(all_probes), ncol = length(sample_names))
rownames(expr_mat) <- all_probes; colnames(expr_mat) <- sample_names
gene_map <- rep("", length(all_probes)); names(gene_map) <- all_probes

for (j in seq_along(sample_names)) {
  dat <- expr_list[[sample_names[j]]]
  m <- match(dat$ProbeName, all_probes)
  expr_mat[m, j] <- dat$gProcessedSignal
  empty <- gene_map[m] == ""
  gene_map[m][empty] <- dat$GeneName[empty]
}
cat("Matrix:", nrow(expr_mat), "x", ncol(expr_mat), "\n")

# Filter + normalize
min_samples <- ceiling(ncol(expr_mat) * 0.25)
keep <- rowSums(expr_mat > 10, na.rm = TRUE) >= min_samples
expr_filt <- expr_mat[keep, ]
expr_norm <- normalizeBetweenArrays(log2(expr_filt + 1), method = "quantile")
cat("Filtered:", nrow(expr_norm), "probes\n")

# 分组
series_gz <- "C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 1 研究设计与转录组差异分析/GSE100927_series_matrix.txt.gz"
con <- gzfile(series_gz, "rt"); series_lines <- readLines(con, warn = FALSE); close(con)
parse_tab <- function(x) { x <- sub("^[^\t]+\t", "", x); strsplit(x, "\t")[[1]] |> gsub(pattern='"', replacement='') }
titles <- parse_tab(series_lines[grep("^!Sample_title", series_lines)])
gsms   <- parse_tab(series_lines[grep("^!Sample_geo_accession", series_lines)])
group_df <- data.frame(GSM=gsms, Title=titles, Group=ifelse(grepl("Control",titles,ignore.case=TRUE),"Control","Atherosclerotic"), stringsAsFactors=FALSE)
group_df <- group_df[group_df$GSM %in% sample_names, ]
group_df <- group_df[match(sample_names, group_df$GSM), ]

saveRDS(list(expr_norm = expr_norm, gene_map = gene_map, group_df = group_df), rds_file)
cat("Saved:", rds_file, "\nDone.\n")
