
# ⚠️ NOTE:
# GSE186989 file is not included in the repository due to size limitation.
# Please download it from GEO and place it under:
# data/GSE186989_CTRL_vs._SU_HX_vs._Chrysin.xlsx

file.path(data_dir, "GSE186989_CTRL_vs._SU_HX_vs._Chrysin.xlsx")



# =========================================================
# Figure 4 — Cross-model and cross-species conservation
# FINAL ABC VERSION
#
# Gene sets are read from:
#   outputs/Fig1_outputs/Figure1_source_data.xlsx
#
# External datasets are read from:
#   data/raw/
#
# Panels:
#   A. Replication rate of Shared / LV-only / RV-only
#   B. Replicated gene counts across datasets
#   C. Rat vs Human summary for Shared genes
#
# Output:
#   outputs/Fig4_outputs/
#     - Figure4_ABC_final.pdf
#     - Figure4_ABC_final.png
#     - Figure4_source_data.xlsx
#     - Fig4A_replication_rate.pdf
#     - Fig4B_replication_counts.pdf
#     - Fig4C_cross_species_summary.pdf
# =========================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(readxl)
  library(readr)
  library(DESeq2)
  library(ggplot2)
  library(patchwork)
  library(AnnotationDbi)
  library(org.Rn.eg.db)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

set.seed(1)

# =========================================================
# 0. Helper: locate script directory
# =========================================================
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
out_dir     <- file.path(output_root, "Fig4_outputs")

if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =========================================================
# 1. Global parameters
# =========================================================
alpha   <- 0.05
lfc_thr <- 1
base_size <- 11

theme_fig <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(face = "bold", color = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

TOUP <- function(x) toupper(trimws(as.character(x)))

clean_sym <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == "" | x == "NA"] <- NA_character_
  x
}

add_sig_col <- function(df, alpha = 0.05, lfc_thr = 1) {
  df %>%
    dplyr::mutate(
      Significant = !is.na(padj) & padj < alpha & abs(log2FC) >= lfc_thr,
      Direction = dplyr::case_when(
        Significant & log2FC > 0  ~ "Up",
        Significant & log2FC < 0  ~ "Down",
        TRUE ~ "NS"
      )
    )
}

safe_round_counts <- function(mat) {
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  round(mat, 0)
}

# =========================================================
# 2. Paths
# =========================================================
fig1_file   <- file.path(output_root, "Fig1_outputs", "Figure1_source_data.xlsx")
file_198618 <- file.path(raw_dir, "GSE198618_Normalized_Counts_RV_ALL.csv")
file_242014 <- file.path(raw_dir, "GSE242014_processed-data-pab.xlsx")
file_240923 <- file.path(raw_dir, "GSE240923_processed-data-mct.xlsx")
file_186989 <- file.path(raw_dir, "GSE186989_CTRL_vs._SU_HX_vs._Chrysin.xlsx")
file_133402 <- file.path(raw_dir, "GSE133402_Rat_RVandLV_HypoxiaNormoxia_Processed.xlsx")
file_240921 <- file.path(raw_dir, "GSE240921_processed-data-human.xlsx")

stopifnot(
  file.exists(fig1_file),
  file.exists(file_198618),
  file.exists(file_242014),
  file.exists(file_240923),
  file.exists(file_186989),
  file.exists(file_133402),
  file.exists(file_240921)
)

# =========================================================
# 3. Read gene sets from Fig1 output
# =========================================================
read_gene_set_from_fig1 <- function(path, sheet_name) {
  df <- openxlsx::read.xlsx(path, sheet = sheet_name)
  gcol <- intersect(c("Genename", "SYMBOL", "Gene", "GeneSymbol", "Gene.name"), names(df))
  if (length(gcol) == 0) {
    stop("No gene-name column found in sheet: ", sheet_name)
  }
  genes <- TOUP(unique(na.omit(df[[gcol[1]]])))
  genes[genes != ""]
}

shared_genes  <- read_gene_set_from_fig1(fig1_file, "Shared_genes")
lv_only_genes <- read_gene_set_from_fig1(fig1_file, "LV_only_genes")
rv_only_genes <- read_gene_set_from_fig1(fig1_file, "RV_only_genes")

cat("Gene set sizes from Fig1:\n")
cat("Shared =", length(shared_genes), "\n")
cat("LV-only =", length(lv_only_genes), "\n")
cat("RV-only =", length(rv_only_genes), "\n")

# =========================================================
# 4. Utility functions
# =========================================================
map_rat_ensembl_to_symbol <- function(ids) {
  ids2 <- sub("\\..*$", "", ids)
  sy <- AnnotationDbi::mapIds(
    org.Rn.eg.db,
    keys = ids2,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  unname(sy[ids2])
}

map_human_ensembl_to_symbol <- function(ids) {
  ids2 <- sub("\\..*$", "", ids)
  sy <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ids2,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  unname(sy[ids2])
}

collapse_by_symbol <- function(mat, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  mat2 <- mat[keep, , drop = FALSE]
  symbols2 <- symbols[keep]
  rownames(mat2) <- symbols2
  
  if (anyDuplicated(rownames(mat2))) {
    mat2 <- rowsum(mat2, group = rownames(mat2)) /
      as.vector(table(rownames(mat2)))
  }
  mat2
}

run_deseq_simple <- function(counts, groups, contrast_vec, is_normalized = FALSE) {
  coldata <- data.frame(condition = factor(groups))
  rownames(coldata) <- colnames(counts)
  
  cmat <- counts
  if (is_normalized) {
    cmat <- round(cmat, 0)
  }
  cmat <- safe_round_counts(cmat)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cmat,
    colData = coldata,
    design = ~ condition
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  
  res <- DESeq2::results(dds, contrast = contrast_vec)
  df <- as.data.frame(res)
  df$Gene <- rownames(df)
  df <- df[!is.na(df$padj), , drop = FALSE]
  
  df %>%
    dplyr::transmute(
      Gene = TOUP(Gene),
      log2FC = log2FoldChange,
      padj = padj
    ) %>%
    add_sig_col(alpha = alpha, lfc_thr = lfc_thr)
}

calc_replication_detail <- function(deg_df, gene_set) {
  gene_set <- unique(TOUP(gene_set))
  
  deg_sig <- deg_df %>%
    dplyr::filter(Significant, Direction %in% c("Up", "Down"))
  
  merged <- dplyr::inner_join(
    tibble::tibble(Gene = gene_set),
    deg_sig %>% dplyr::select(Gene, log2FC, padj, Direction, dataset),
    by = "Gene"
  )
  
  list(
    n_set = length(gene_set),
    n_rep = dplyr::n_distinct(merged$Gene),
    rate = ifelse(length(gene_set) == 0, NA_real_,
                  dplyr::n_distinct(merged$Gene) / length(gene_set)),
    genes = unique(merged$Gene),
    merged = merged
  )
}

# =========================================================
# 5. External dataset readers
# =========================================================

# ---------- GSE198618 ----------
read_GSE198618 <- function(file) {
  expr_raw <- read.csv(file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  ens_raw  <- expr_raw[[1]]
  mat      <- as.matrix(expr_raw[, -1, drop = FALSE])
  rownames(mat) <- sub("\\..*$", "", ens_raw)
  mode(mat) <- "numeric"
  
  sy <- map_human_ensembl_to_symbol(rownames(mat))
  mat_sym <- collapse_by_symbol(mat, sy)
  
  samples <- colnames(mat_sym)
  group <- rep(NA_character_, length(samples))
  group[grepl("(?i)(^control|^ctrl$|normal|healthy|rv[-_ ]?normal)", samples, perl = TRUE)] <- "Ctrl"
  group[grepl("(?i)(compensated|^crv$|rv[-_ ]?comp)", samples, perl = TRUE)] <- "cRV"
  group[grepl("(?i)(decomp|^drv$|rv[-_ ]?fail|fail|failure)", samples, perl = TRUE)] <- "dRV"
  
  if (any(is.na(group))) {
    stop("GSE198618: unmatched sample names:\n", paste(samples[is.na(group)], collapse = ", "))
  }
  
  deg <- run_deseq_simple(mat_sym, group, c("condition", "dRV", "Ctrl"), is_normalized = TRUE)
  deg$dataset <- "GSE198618"
  deg
}

# ---------- GSE242014 ----------
read_GSE242014 <- function(file) {
  expr_pab <- openxlsx::read.xlsx(file, sheet = 2)
  ens <- expr_pab[[1]]
  mat <- as.matrix(expr_pab[, -1, drop = FALSE])
  mode(mat) <- "numeric"
  rownames(mat) <- ens
  
  group <- ifelse(grepl("^CTRL", colnames(mat), ignore.case = TRUE), "Ctrl",
                  ifelse(grepl("^Comp", colnames(mat), ignore.case = TRUE), "cRV",
                         ifelse(grepl("^Decomp", colnames(mat), ignore.case = TRUE), "dRV", NA_character_)))
  if (any(is.na(group))) stop("GSE242014: unmatched sample names")
  
  sy <- map_rat_ensembl_to_symbol(rownames(mat))
  mat_sym <- collapse_by_symbol(mat, sy)
  
  deg <- run_deseq_simple(mat_sym, group, c("condition", "dRV", "Ctrl"), is_normalized = FALSE)
  deg$dataset <- "GSE242014"
  deg
}

# ---------- GSE240923 batch1 ----------
read_GSE240923_batch1 <- function(file) {
  count1 <- openxlsx::read.xlsx(file, sheet = 2, rowNames = TRUE)
  count1 <- as.matrix(count1)
  mode(count1) <- "numeric"
  
  sy <- map_rat_ensembl_to_symbol(rownames(count1))
  mat_sym <- collapse_by_symbol(count1, sy)
  
  samples <- colnames(mat_sym)
  group <- dplyr::case_when(
    grepl("^CTRL", samples) ~ "Ctrl",
    grepl("^COMP", samples) ~ "cRV",
    grepl("^DECOMP-FP", samples) ~ "dRV",
    grepl("^DECOMP-MCT69", samples) ~ "dRV",
    TRUE ~ NA_character_
  )
  if (any(is.na(group))) stop("GSE240923 batch1: unmatched sample names")
  
  deg <- run_deseq_simple(mat_sym, group, c("condition", "dRV", "Ctrl"), is_normalized = FALSE)
  deg$dataset <- "GSE240923_batch1"
  deg
}

# ---------- GSE240923 batch2 ----------
read_GSE240923_batch2 <- function(file) {
  count2 <- openxlsx::read.xlsx(file, sheet = 3, rowNames = TRUE)
  count2 <- as.matrix(count2)
  mode(count2) <- "numeric"
  
  sy <- map_rat_ensembl_to_symbol(rownames(count2))
  mat_sym <- collapse_by_symbol(count2, sy)
  
  samples <- colnames(mat_sym)
  group <- dplyr::case_when(
    grepl("^CTRL", samples, ignore.case = TRUE) ~ "Ctrl",
    grepl("^COMP", samples, ignore.case = TRUE) ~ "cRV",
    grepl("^DECOMP", samples, ignore.case = TRUE) ~ "dRV",
    TRUE ~ NA_character_
  )
  sex <- dplyr::case_when(
    grepl("-f(_|$)", samples, ignore.case = TRUE) ~ "F",
    grepl("-m(_|$)", samples, ignore.case = TRUE) ~ "M",
    TRUE ~ NA_character_
  )
  
  if (any(is.na(group)) || any(is.na(sex))) stop("GSE240923 batch2: unmatched sample names")
  
  coldata <- data.frame(
    Sex = factor(sex, levels = c("F", "M")),
    Group = factor(group, levels = c("Ctrl", "cRV", "dRV"))
  )
  rownames(coldata) <- colnames(mat_sym)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = safe_round_counts(mat_sym),
    colData = coldata,
    design = ~ Sex + Group
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  
  res <- DESeq2::results(dds, contrast = c("Group", "dRV", "Ctrl"))
  deg <- as.data.frame(res)
  deg$Gene <- rownames(deg)
  deg <- deg[!is.na(deg$padj), , drop = FALSE] %>%
    dplyr::transmute(
      Gene = TOUP(Gene),
      log2FC = log2FoldChange,
      padj = padj
    ) %>%
    add_sig_col(alpha = alpha, lfc_thr = lfc_thr)
  
  deg$dataset <- "GSE240923_batch2"
  deg
}

# ---------- GSE186989 ----------
read_GSE186989 <- function(file) {
  df <- readxl::read_excel(file, sheet = 1)
  deg <- df %>%
    dplyr::transmute(
      Gene = TOUP(`Feature ID`),
      log2FC = as.numeric(`EDGE test: CTRL vs SU/HX, tagwise dispersions - Fold change`),
      padj   = as.numeric(`EDGE test: CTRL vs SU/HX, tagwise dispersions - FDR p-value correction`)
    ) %>%
    dplyr::filter(!is.na(Gene), Gene != "", !is.na(log2FC), !is.na(padj)) %>%
    add_sig_col(alpha = alpha, lfc_thr = lfc_thr)
  
  deg$dataset <- "GSE186989"
  deg
}

# ---------- GSE133402 ----------
read_GSE133402 <- function(file) {
  raw <- readxl::read_excel(file, sheet = 1, col_names = FALSE)
  
  hdr <- as.character(unlist(raw[1, ]))
  count_idx <- which(tolower(hdr) == "counts")
  sample_idx <- count_idx - 2
  valid <- sample_idx >= 1
  
  count_idx <- count_idx[valid]
  sample_idx <- sample_idx[valid]
  
  sample_names <- trimws(hdr[sample_idx])
  
  gene_col_idx <- sample_idx[1]
  row_start <- 2
  
  gene_ids <- as.character(raw[[gene_col_idx]][row_start:nrow(raw)])
  count_mat <- sapply(count_idx, function(j) as.numeric(raw[[j]][row_start:nrow(raw)]))
  count_mat <- as.matrix(count_mat)
  colnames(count_mat) <- sample_names
  rownames(count_mat) <- gene_ids
  
  keep <- !is.na(rownames(count_mat)) & rownames(count_mat) != ""
  count_mat <- count_mat[keep, , drop = FALSE]
  
  rn <- rownames(count_mat)
  is_ens <- mean(grepl("^ENS", rn)) > 0.5
  if (is_ens) {
    sy <- map_rat_ensembl_to_symbol(rn)
  } else {
    sy <- rn
  }
  mat_sym <- collapse_by_symbol(count_mat, TOUP(sy))
  
  group <- ifelse(grepl("RMCH", colnames(mat_sym), ignore.case = TRUE), "Hypoxia",
                  ifelse(grepl("RMNT", colnames(mat_sym), ignore.case = TRUE), "Ctrl", NA_character_))
  chamber <- ifelse(grepl("LV", colnames(mat_sym), ignore.case = TRUE), "LV",
                    ifelse(grepl("RV", colnames(mat_sym), ignore.case = TRUE), "RV", NA_character_))
  
  if (any(is.na(group)) || any(is.na(chamber))) {
    stop("GSE133402: sample names could not be parsed")
  }
  
  keep_rv <- chamber == "RV"
  deg <- run_deseq_simple(
    mat_sym[, keep_rv, drop = FALSE],
    group[keep_rv],
    c("condition", "Hypoxia", "Ctrl"),
    is_normalized = FALSE
  )
  deg$dataset <- "GSE133402_RV"
  deg
}

# ---------- GSE240921 ----------
read_GSE240921 <- function(file) {
  df <- readxl::read_excel(file, sheet = 2)
  
  expr <- as.matrix(df[, -1, drop = FALSE])
  rownames(expr) <- as.character(df[[1]])
  mode(expr) <- "numeric"
  
  samp <- colnames(expr)
  use_rv <- grepl("^RV-", samp)
  if (sum(use_rv) >= 6) {
    expr <- expr[, use_rv, drop = FALSE]
    samp <- colnames(expr)
  }
  
  group <- ifelse(grepl("Control|RV-Normal", samp, ignore.case = TRUE), "Ctrl",
                  ifelse(grepl("Compensated|RV-Compen", samp, ignore.case = TRUE), "cRV",
                         ifelse(grepl("Decompensated|RV-Failing", samp, ignore.case = TRUE), "dRV", NA_character_)))
  
  if (any(is.na(group))) {
    stop("GSE240921: unmatched sample names")
  }
  
  rn <- rownames(expr)
  if (mean(grepl("^ENS", rn)) > 0.5) {
    sy <- map_human_ensembl_to_symbol(rn)
  } else {
    sy <- rn
  }
  expr_sym <- collapse_by_symbol(expr, TOUP(sy))
  
  deg <- run_deseq_simple(expr_sym, group, c("condition", "dRV", "Ctrl"), is_normalized = FALSE)
  deg$dataset <- "GSE240921"
  deg
}

# =========================================================
# 6. Run all readers
# =========================================================
deg_198618    <- read_GSE198618(file_198618)
deg_242014    <- read_GSE242014(file_242014)
deg_240923_b1 <- read_GSE240923_batch1(file_240923)
deg_240923_b2 <- read_GSE240923_batch2(file_240923)
deg_186989    <- read_GSE186989(file_186989)
deg_133402    <- read_GSE133402(file_133402)
deg_240921    <- read_GSE240921(file_240921)

all_deg <- dplyr::bind_rows(
  deg_198618,
  deg_242014,
  deg_240923_b1,
  deg_240923_b2,
  deg_186989,
  deg_133402,
  deg_240921
)

cat("\nDatasets loaded:\n")
print(table(all_deg$dataset))

# =========================================================
# 7. Replication calculations
# =========================================================
external_datasets <- unique(all_deg$dataset)

rep_list <- list()

for (ds in external_datasets) {
  deg_df <- all_deg %>% dplyr::filter(dataset == ds)
  
  r1 <- calc_replication_detail(deg_df, shared_genes)
  r2 <- calc_replication_detail(deg_df, lv_only_genes)
  r3 <- calc_replication_detail(deg_df, rv_only_genes)
  
  rep_list[[ds]] <- data.frame(
    dataset = ds,
    Shared_rate = r1$rate,
    LV_only_rate = r2$rate,
    RV_only_rate = r3$rate,
    Shared_n = r1$n_rep,
    LV_only_n = r2$n_rep,
    RV_only_n = r3$n_rep,
    stringsAsFactors = FALSE
  )
}

rep_df <- dplyr::bind_rows(rep_list)
rep_df$Species <- ifelse(grepl("198618|240921", rep_df$dataset), "Human", "Rat")

# =========================================================
# 8. Panel A — replication rate
# =========================================================
plotA_df <- rep_df %>%
  dplyr::select(dataset, Shared_rate, LV_only_rate, RV_only_rate) %>%
  tidyr::pivot_longer(-dataset, names_to = "Set", values_to = "Rate") %>%
  dplyr::mutate(
    Set = factor(
      Set,
      levels = c("Shared_rate", "LV_only_rate", "RV_only_rate"),
      labels = c("Shared", "LV-only", "RV-only")
    ),
    dataset = factor(dataset, levels = external_datasets)
  )

pA <- ggplot2::ggplot(plotA_df, ggplot2::aes(x = dataset, y = Rate, fill = Set)) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.72),
    width = 0.68,
    color = "black",
    linewidth = 0.2
  ) +
  ggplot2::scale_fill_manual(values = c("Shared" = "#D73027", "LV-only" = "#1F78B4", "RV-only" = "#33A02C")) +
  ggplot2::labs(
    title = "A. Replication rate of ventricular gene sets",
    x = NULL,
    y = "Same-direction replication rate"
  ) +
  theme_fig(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# =========================================================
# 9. Panel B — replicated counts
# =========================================================
plotB_df <- rep_df %>%
  dplyr::select(dataset, Shared_n, LV_only_n, RV_only_n) %>%
  tidyr::pivot_longer(-dataset, names_to = "Set", values_to = "Count") %>%
  dplyr::mutate(
    Set = factor(
      Set,
      levels = c("Shared_n", "LV_only_n", "RV_only_n"),
      labels = c("Shared", "LV-only", "RV-only")
    ),
    dataset = factor(dataset, levels = external_datasets)
  )

pB <- ggplot2::ggplot(plotB_df, ggplot2::aes(x = dataset, y = Count, fill = Set)) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.72),
    width = 0.68,
    color = "black",
    linewidth = 0.2
  ) +
  ggplot2::scale_fill_manual(values = c("Shared" = "#D73027", "LV-only" = "#1F78B4", "RV-only" = "#33A02C")) +
  ggplot2::labs(
    title = "B. Replicated gene counts across datasets",
    x = NULL,
    y = "Replicated gene number"
  ) +
  theme_fig(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# =========================================================
# 10. Panel C — Rat vs Human summary
# =========================================================
plotC_df <- rep_df %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    Mean_rate = mean(Shared_rate, na.rm = TRUE),
    Mean_n = mean(Shared_n, na.rm = TRUE),
    .groups = "drop"
  )

pC <- ggplot2::ggplot(plotC_df, ggplot2::aes(x = Species, y = Mean_rate, fill = Species)) +
  ggplot2::geom_col(width = 0.65, color = "black", linewidth = 0.2) +
  ggplot2::geom_text(ggplot2::aes(label = round(Mean_rate, 2)), vjust = -0.35, size = 4) +
  ggplot2::scale_fill_manual(values = c("Rat" = "#4DAF4A", "Human" = "#984EA3")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
  ggplot2::labs(
    title = "C. Cross-species support of the Shared set",
    x = NULL,
    y = "Mean replication rate"
  ) +
  theme_fig(base_size = 10) +
  ggplot2::theme(legend.position = "none")

# =========================================================
# 11. Compose Figure 4 ABC
# =========================================================
fig4_ABC <- (pA + pB + pC) +
  patchwork::plot_layout(widths = c(1.15, 1.15, 0.8)) +
  patchwork::plot_annotation(
    title = "Figure 4. Cross-model and cross-species conservation of ventricular remodeling gene sets",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

# =========================================================
# 12. Save figure
# =========================================================
ggplot2::ggsave(file.path(out_dir, "Figure4_ABC_final.pdf"), fig4_ABC, width = 13, height = 4.2)
ggplot2::ggsave(file.path(out_dir, "Figure4_ABC_final.png"), fig4_ABC, width = 13, height = 4.2, dpi = 300)

ggplot2::ggsave(file.path(out_dir, "Fig4A_replication_rate.pdf"), pA, width = 3.8, height = 3.3)
ggplot2::ggsave(file.path(out_dir, "Fig4B_replication_counts.pdf"), pB, width = 3.8, height = 3.3)
ggplot2::ggsave(file.path(out_dir, "Fig4C_cross_species_summary.pdf"), pC, width = 2.6, height = 3.0)

# =========================================================
# 13. Export source data
# =========================================================
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Shared_genes")
openxlsx::writeData(wb, "Shared_genes", data.frame(Gene = shared_genes))

openxlsx::addWorksheet(wb, "LV_only_genes")
openxlsx::writeData(wb, "LV_only_genes", data.frame(Gene = lv_only_genes))

openxlsx::addWorksheet(wb, "RV_only_genes")
openxlsx::writeData(wb, "RV_only_genes", data.frame(Gene = rv_only_genes))

openxlsx::addWorksheet(wb, "All_external_DEG")
openxlsx::writeData(wb, "All_external_DEG", all_deg)

openxlsx::addWorksheet(wb, "Replication_summary")
openxlsx::writeData(wb, "Replication_summary", rep_df)

openxlsx::addWorksheet(wb, "PanelA_data")
openxlsx::writeData(wb, "PanelA_data", plotA_df)

openxlsx::addWorksheet(wb, "PanelB_data")
openxlsx::writeData(wb, "PanelB_data", plotB_df)

openxlsx::addWorksheet(wb, "PanelC_data")
openxlsx::writeData(wb, "PanelC_data", plotC_df)

openxlsx::addWorksheet(wb, "Meta")
openxlsx::writeData(
  wb, "Meta",
  data.frame(
    Parameter = c("alpha", "lfc_thr", "Fig1_input"),
    Value = c(alpha, lfc_thr, fig1_file)
  )
)

openxlsx::addWorksheet(wb, "sessionInfo")
openxlsx::writeData(wb, "sessionInfo", capture.output(sessionInfo()))

openxlsx::saveWorkbook(wb, file.path(out_dir, "Figure4_source_data.xlsx"), overwrite = TRUE)

message("✅ Figure 4 ABC finished. Outputs saved to: ", out_dir)