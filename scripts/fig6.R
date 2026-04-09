# ================================================================
# Fig6 — CTEPH external validation heatmap using Fig5 consensus genes
#
# Main logic:
#   1) use relative paths only
#   2) signature genes are read from:
#        outputs/Fig5_outputs/Figure5_source_data.xlsx -> sheet "Consensus_genes"
#   3) plotting style kept consistent with your previous code
#   4) remove duplicated comparison Fig 4f
#   5) reverse Fig 4e direction so it becomes:
#        B-prePEA_septum vs B-prePEA_RV
#
# Output:
#   outputs/Fig6_outputs/
#     1) Fig6_43gene_detailed_clean.xlsx
#     2) Fig6_43gene_heatmap_sigcount.pdf
#     3) Fig6_43gene_heatmap_sigcount.png
#     4) Fig6_heatmap_matrices.xlsx
# ================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

# ------------------------------------------------
# 0. Auto-detect script directory
# ------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  idx <- grep(file_arg, args)
  
  if (length(idx) > 0) {
    path <- sub(file_arg, "", args[idx[1]])
    return(dirname(normalizePath(path)))
  }
  
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(dirname(normalizePath(ctx$path)))
    }
  }
  
  normalizePath(getwd())
}

script_dir  <- get_script_dir()
project_dir <- normalizePath(file.path(script_dir, ".."))
raw_dir     <- file.path(project_dir, "data", "raw")
output_root <- file.path(project_dir, "outputs")

# 输出目录
outdir <- file.path(output_root, "Fig6_outputs")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ------------------------------------------------
# 1. Read consensus genes from Fig5 source data
# ------------------------------------------------
fig5_dir   <- file.path(output_root, "Fig5_outputs")
fig5_file  <- file.path(fig5_dir, "Figure5_source_data.xlsx")
fig5_sheet <- "Consensus_genes"

if (!file.exists(fig5_file)) {
  stop("Cannot find Fig5 source data file: ", fig5_file)
}

sig_tab <- read.xlsx(fig5_file, sheet = fig5_sheet)

# 优先用 Gene 列，否则用第一列
gcol <- if ("Gene" %in% colnames(sig_tab)) "Gene" else colnames(sig_tab)[1]

signature_genes <- toupper(na.omit(trimws(as.character(sig_tab[[gcol]]))))
signature_genes <- unique(signature_genes)

cat("Number of signature genes from Fig5 Consensus_genes:", length(signature_genes), "\n")
print(head(signature_genes, 10))

# ------------------------------------------------
# 2. Input files and valid sheets (relative paths)
# ------------------------------------------------
files <- list(
  Fig2 = file.path(raw_dir, "44161_2025_672_MOESM4_ESM.xlsx"),
  Fig3 = file.path(raw_dir, "44161_2025_672_MOESM5_ESM.xlsx"),
  Fig4 = file.path(raw_dir, "44161_2025_672_MOESM6_ESM.xlsx"),
  Fig5 = file.path(raw_dir, "44161_2025_672_MOESM7_ESM.xlsx")
)

# 已删除重复的 Fig 4f
sheets_map <- list(
  Fig2 = c("Fig 2c", "Fig 2d", "Fig 2e"),
  Fig3 = c("Fig 3c", "Fig 3d", "Fig 3e"),
  Fig4 = c("Fig 4c", "Fig 4e"),
  Fig5 = c("Fig 5c", "Fig 5d", "Fig 5e", "Fig 5f")
)

# 检查文件是否存在
for (nm in names(files)) {
  if (!file.exists(files[[nm]])) {
    stop("Cannot find file: ", files[[nm]])
  }
}

# ------------------------------------------------
# 3. Contrast labels
#    Fig 4e is displayed in reversed direction:
#    original: B-prePEA_RV / B-prePEA_septum
#    display : B-prePEA_septum vs B-prePEA_RV
# ------------------------------------------------
contrast_labels <- c(
  "Fig 2c" = "B-prePEAm_RV vs B-prePEAs_RV",
  "Fig 2d" = "B-prePEAm_RV vs B-prePEAi_RV",
  "Fig 2e" = "B-prePEAi_RV vs B-prePEAs_RV",
  "Fig 3c" = "B-prePEAm_RV vs Control_RV",
  "Fig 3d" = "B-prePEAi_RV vs Control_RV",
  "Fig 3e" = "B-prePEAs_RV vs Control_RV",
  "Fig 4c" = "B-prePEA_septum vs Control_septum",
  "Fig 4e" = "B-prePEA_septum vs B-prePEA_RV",
  "Fig 5c" = "B-postPEA_septum vs B-prePEA_RV",
  "Fig 5d" = "B-postPEAm_septum vs B-prePEAm_RV",
  "Fig 5e" = "B-postPEAi_septum vs B-prePEAi_RV",
  "Fig 5f" = "B-postPEAs_septum vs B-prePEAs_RV"
)

contrast_order <- unname(contrast_labels)

# 需要翻转 log2FC 方向的 contrast
reverse_contrasts <- c("Fig 4e")

# ------------------------------------------------
# 4. Helper: extract target genes from one sheet
#    force priority to use Ensembl.gene if present
# ------------------------------------------------
extract_sig_genes_full <- function(df, contrast_name, signature_genes, reverse_contrasts = NULL) {
  names(df) <- trimws(names(df))
  
  # 优先选择 Ensembl.gene；否则找 Gene / SYMBOL
  if ("Ensembl.gene" %in% names(df)) {
    gene_col <- "Ensembl.gene"
  } else {
    gene_col <- grep("Gene|SYMBOL", names(df), value = TRUE, ignore.case = TRUE)[1]
  }
  
  lfc_col  <- grep("log2|logFC", names(df), value = TRUE, ignore.case = TRUE)[1]
  padj_col <- grep("padj|FDR|adj", names(df), value = TRUE, ignore.case = TRUE)[1]
  bm_col   <- grep("baseMean", names(df), value = TRUE, ignore.case = TRUE)[1]
  
  if (is.na(gene_col) || length(gene_col) == 0) {
    stop("No gene column found in contrast: ", contrast_name)
  }
  if (is.na(lfc_col) || length(lfc_col) == 0) {
    stop("No log2FC column found in contrast: ", contrast_name)
  }
  if (is.na(padj_col) || length(padj_col) == 0) {
    stop("No adjusted p-value column found in contrast: ", contrast_name)
  }
  
  if (is.na(bm_col) || length(bm_col) == 0) {
    df$baseMean <- NA
    bm_col <- "baseMean"
  }
  
  df2 <- df %>%
    dplyr::select(
      Gene = !!sym(gene_col),
      baseMean = !!sym(bm_col),
      log2FoldChange = !!sym(lfc_col),
      padj = !!sym(padj_col)
    ) %>%
    dplyr::mutate(
      Gene = toupper(trimws(as.character(Gene))),
      baseMean = suppressWarnings(as.numeric(baseMean)),
      log2FoldChange = suppressWarnings(as.numeric(log2FoldChange)),
      padj = suppressWarnings(as.numeric(padj)),
      Contrast = contrast_name
    ) %>%
    dplyr::filter(!is.na(Gene), Gene != "", Gene %in% signature_genes)
  
  # 若某个基因重复，保留 padj 最小的一条
  df2 <- df2 %>%
    dplyr::arrange(Gene, padj, dplyr::desc(abs(log2FoldChange))) %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  # 若需要按展示方向翻转 log2FC，则乘以 -1
  if (!is.null(reverse_contrasts) && contrast_name %in% reverse_contrasts) {
    df2 <- df2 %>%
      dplyr::mutate(log2FoldChange = -log2FoldChange)
  }
  
  return(df2)
}

# ------------------------------------------------
# 5. Traverse all sheets
# ------------------------------------------------
all_sig_full <- list()

for (fig in names(files)) {
  for (sh in sheets_map[[fig]]) {
    cat("Reading:", fig, "-", sh, "\n")
    df <- read.xlsx(files[[fig]], sheet = sh)
    subdf <- extract_sig_genes_full(
      df = df,
      contrast_name = sh,
      signature_genes = signature_genes,
      reverse_contrasts = reverse_contrasts
    )
    all_sig_full[[paste(fig, sh, sep = ":")]] <- subdf
  }
}

all_sig_full <- bind_rows(all_sig_full)

# 增加直观 contrast 标签
all_sig_full$Contrast_label <- contrast_labels[all_sig_full$Contrast]

# 按原始顺序整理
all_sig_full$Contrast_label <- factor(all_sig_full$Contrast_label, levels = contrast_order)

# ------------------------------------------------
# 6. Fill missing gene-contrast combinations
#    (to ensure all consensus genes appear in heatmap)
# ------------------------------------------------
full_grid <- expand.grid(
  Gene = signature_genes,
  Contrast_label = contrast_order,
  stringsAsFactors = FALSE
)

df_full <- full_grid %>%
  dplyr::left_join(
    all_sig_full %>%
      dplyr::select(Gene, Contrast_label, baseMean, log2FoldChange, padj),
    by = c("Gene", "Contrast_label")
  )

# ------------------------------------------------
# 7. Export detailed Excel
# ------------------------------------------------
outfile_xlsx <- file.path(outdir, "Fig6_43gene_detailed_clean.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "AllData_raw_found")
writeData(wb, "AllData_raw_found", all_sig_full)

addWorksheet(wb, "AllData_full_grid")
writeData(wb, "AllData_full_grid", df_full)

addWorksheet(wb, "Signature_genes")
writeData(wb, "Signature_genes", data.frame(Gene = signature_genes))

saveWorkbook(wb, outfile_xlsx, overwrite = TRUE)

cat("✅ Detailed table saved to:", outfile_xlsx, "\n")

# ------------------------------------------------
# 8. Build log2FC matrix
# ------------------------------------------------
mat <- df_full %>%
  dplyr::select(Gene, Contrast_label, log2FoldChange) %>%
  tidyr::pivot_wider(names_from = Contrast_label, values_from = log2FoldChange) %>%
  as.data.frame()

rownames(mat) <- mat$Gene
mat <- mat[, -1, drop = FALSE]
mat <- as.matrix(mat)

# 保证列顺序固定
mat <- mat[, contrast_order, drop = FALSE]

# 缺失值保留，后面 scale 后再处理
mat_scaled <- t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)] <- 0

# ------------------------------------------------
# 9. Build significance star matrix
# ------------------------------------------------
sig_mat <- df_full %>%
  dplyr::select(Gene, Contrast_label, padj) %>%
  tidyr::pivot_wider(names_from = Contrast_label, values_from = padj) %>%
  as.data.frame()

rownames(sig_mat) <- sig_mat$Gene
sig_mat <- sig_mat[, -1, drop = FALSE]
sig_mat <- as.matrix(sig_mat)

sig_mat <- sig_mat[, contrast_order, drop = FALSE]
sig_mat_num <- apply(sig_mat, 2, as.numeric)

sig_mat_chr <- ifelse(
  is.na(sig_mat_num), "",
  ifelse(sig_mat_num < 0.001, "***",
         ifelse(sig_mat_num < 0.01, "**",
                ifelse(sig_mat_num < 0.05, "*", "")))
)

sig_mat_chr <- matrix(
  sig_mat_chr,
  nrow = nrow(sig_mat_num),
  ncol = ncol(sig_mat_num),
  dimnames = dimnames(sig_mat)
)

sig_mat_chr <- sig_mat_chr[rownames(mat_scaled), colnames(mat_scaled), drop = FALSE]

# ------------------------------------------------
# 10. Count significant contrasts for each gene
#     significance: padj < 0.05 & |log2FC| >= 1
# ------------------------------------------------
n_cols <- ncol(mat_scaled)

sig_count_df <- df_full %>%
  dplyr::mutate(is_sig = !is.na(padj) & padj < 0.05 & !is.na(log2FoldChange) & abs(log2FoldChange) >= 1) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(sig_count = sum(is_sig, na.rm = TRUE), .groups = "drop")

rn <- rownames(mat_scaled)
sig_vec <- sig_count_df$sig_count[match(rn, sig_count_df$Gene)]
sig_vec[is.na(sig_vec)] <- 0L

labels_row <- paste0(rn, " (", sig_vec, "/", n_cols, ")")

annotation_row <- data.frame(`Sig.count` = sig_vec)
rownames(annotation_row) <- rn

ann_colors <- list(
  `Sig.count` = colorRampPalette(c("#f0f0f0", "#b81f25"))(50)
)

# ------------------------------------------------
# 11. Colors and row clustering
# ------------------------------------------------
my_colors <- colorRampPalette(c("blue", "#ffff7f", "orange", "red"))(50)

# ------------------------------------------------
# 12. Plot heatmap
# ------------------------------------------------
outfile_pdf <- file.path(outdir, "Fig6_43gene_heatmap_sigcount.pdf")
outfile_png <- file.path(outdir, "Fig6_43gene_heatmap_sigcount.png")

pdf(outfile_pdf, width = 6, height = 6.5)
pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_row = labels_row,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  color = my_colors,
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 7,
  angle_col = 90,
  main = paste0("CTEPH ", length(signature_genes), "-gene signature across contrasts"),
  display_numbers = sig_mat_chr,
  number_color = "black",
  fontsize_number = 8
)
dev.off()

png(outfile_png, width = 6, height = 6.5, units = "in", res = 300)
pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_row = labels_row,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  color = my_colors,
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 7,
  angle_col = 90,
  main = paste0("CTEPH ", length(signature_genes), "-gene signature across contrasts"),
  display_numbers = sig_mat_chr,
  number_color = "black",
  fontsize_number = 8
)
dev.off()

cat("✅ Heatmap PDF saved to:", outfile_pdf, "\n")
cat("✅ Heatmap PNG saved to:", outfile_png, "\n")

# ------------------------------------------------
# 13. Also export the final matrices for checking
# ------------------------------------------------
outfile_matrix <- file.path(outdir, "Fig6_heatmap_matrices.xlsx")

wb2 <- createWorkbook()

addWorksheet(wb2, "mat_log2FC")
writeData(wb2, "mat_log2FC", data.frame(Gene = rownames(mat), mat, check.names = FALSE))

addWorksheet(wb2, "mat_scaled")
writeData(wb2, "mat_scaled", data.frame(Gene = rownames(mat_scaled), mat_scaled, check.names = FALSE))

addWorksheet(wb2, "sig_stars")
writeData(wb2, "sig_stars", data.frame(Gene = rownames(sig_mat_chr), sig_mat_chr, check.names = FALSE))

addWorksheet(wb2, "sig_count")
writeData(wb2, "sig_count", data.frame(
  Gene = rn,
  Sig_count = sig_vec,
  Row_label = labels_row
))

saveWorkbook(wb2, outfile_matrix, overwrite = TRUE)

cat("✅ Heatmap matrices saved to:", outfile_matrix, "\n")
