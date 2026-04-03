# ============================================
# Figure 1 — Predominantly shared ventricular remodeling in PH
# Discovery dataset: GSE266139 (MCT)
# Threshold: padj < 0.05 and |log2FC| >= 1
# GitHub-ready version
# Panels:
#   A. DEG counts across 4 contrasts
#   B. Correlation of LV and RV responses
#   C. DEG set composition
#   D. Heatmap of top shared genes
#
# Output:
#   outputs/Fig1_outputs/
#     - Figure1_final.pdf
#     - Figure1_final.png
#     - Figure1_source_data.xlsx
#     - Fig1A_DEG_counts.pdf
#     - Fig1B_correlation.pdf
#     - Fig1C_set_composition.pdf
#     - Fig1D_heatmap.pdf
#
# Also copy:
#   data/processed/Figure1_source_data.xlsx
# ============================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(ggplotify)
  library(grid)
  library(scales)
  library(matrixStats)
})

set.seed(1)

# ---------- Helper: locate script directory ----------
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  idx <- grep(file_arg, cmd_args)
  
  if (length(idx) > 0) {
    path <- sub(file_arg, "", cmd_args[idx[1]])
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

script_dir   <- get_script_dir()
project_dir  <- normalizePath(file.path(script_dir, ".."))
raw_dir      <- file.path(project_dir, "data", "raw")
processed_dir <- file.path(project_dir, "data", "processed")
output_root  <- file.path(project_dir, "outputs")

if (!dir.exists(processed_dir)) dir.create(processed_dir, recursive = TRUE)
if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)

# ---------- Paths ----------
# 建议你把 raw 里的文件统一命名成这个
input_file <- file.path(raw_dir, "GSE266139_12967_2025_6792_MOESM2_ESM.xlsx")

# 如果你当前 raw 文件夹里还是旧名字，就自动兼容一下
if (!file.exists(input_file)) {
  fallback_file <- file.path(raw_dir, "12967_2025_6792_MOESM2_ESM.xlsx")
  if (file.exists(fallback_file)) {
    input_file <- fallback_file
  }
}

outdir <- file.path(output_root, "Fig1_outputs")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# ---------- Parameters ----------
alpha   <- 0.05
lfc_thr <- 1
n_heat  <- 18

# ---------- Colors ----------
my_group_cols <- c(
  "LP" = "#E41A1C",
  "RP" = "#377EB8",
  "LM" = "#4DAF4A",
  "RM" = "#984EA3"
)

my_set_cols <- c(
  "NS"           = "grey80",
  "Shared"       = "#D73027",
  "LV-only"      = "#1F78B4",
  "RV-only"      = "#33A02C",
  "Opposite-dir" = "#FF7F00"
)

# ---------- Theme ----------
theme_fig <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

# ---------- Helpers ----------
clean_sym <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == "" | x == "NA"] <- NA_character_
  x
}

deg_count <- function(res, alpha = 0.05, lfc_thr = 1) {
  up   <- sum(res$log2FoldChange >  lfc_thr & res$padj < alpha, na.rm = TRUE)
  down <- sum(res$log2FoldChange < -lfc_thr & res$padj < alpha, na.rm = TRUE)
  data.frame(Up = up, Down = down, Total = up + down)
}

# ---------- Read data ----------
expr_full <- openxlsx::read.xlsx(input_file, sheet = 1)
stopifnot(all(c("Gene.stable.ID", "Gene.name") %in% colnames(expr_full)))

gene_annot <- expr_full[, c("Gene.stable.ID", "Gene.name")]
colnames(gene_annot) <- c("Geneid", "Genename")
gene_annot$Geneid   <- clean_sym(gene_annot$Geneid)
gene_annot$Genename <- clean_sym(gene_annot$Genename)

# ---------- Strictly select count columns only ----------
sample_cols <- grep("^(LP|LM|RP|RM)_\\d+$", colnames(expr_full), value = TRUE)
if (length(sample_cols) == 0) {
  stop("No sample columns found matching LP_1 / LM_1 / RP_1 / RM_1 style names.")
}

cat("Selected sample columns:\n")
print(sample_cols)
cat("Number of samples:", length(sample_cols), "\n")

expr <- expr_full[, sample_cols, drop = FALSE]

# force numeric
expr[] <- lapply(expr, function(x) as.numeric(as.character(x)))

# check NA
na_per_col <- colSums(is.na(expr))
cat("NA count per sample column:\n")
print(na_per_col)

na_rows <- rowSums(is.na(expr)) > 0
cat("Genes with any NA across selected samples:", sum(na_rows), "\n")

if (sum(na_rows) > 0) {
  expr <- expr[!na_rows, , drop = FALSE]
  gene_annot <- gene_annot[!na_rows, , drop = FALSE]
}

# build matrix
expr_mat <- as.matrix(expr)
storage.mode(expr_mat) <- "numeric"

if (any(expr_mat < 0, na.rm = TRUE)) {
  stop("Negative values detected in count matrix.")
}
if (any(is.na(expr_mat))) {
  stop("NA values still remain in count matrix after filtering.")
}

expr_mat <- round(expr_mat)
storage.mode(expr_mat) <- "integer"
rownames(expr_mat) <- gene_annot$Geneid

cat("Final matrix dimension:", dim(expr_mat)[1], "genes x", dim(expr_mat)[2], "samples\n")
cat("Any NA left? ", any(is.na(expr_mat)), "\n")

# ---------- Sample annotation ----------
samples <- colnames(expr_mat)
group <- ifelse(grepl("^LP_", samples), "LP",
                ifelse(grepl("^RP_", samples), "RP",
                       ifelse(grepl("^LM_", samples), "LM",
                              ifelse(grepl("^RM_", samples), "RM", NA_character_))))

stopifnot(!any(is.na(group)))

coldata <- data.frame(
  sample = samples,
  condition = factor(group, levels = c("LP", "RP", "LM", "RM"))
)
rownames(coldata) <- samples

# ---------- Low-count filter ----------
# filter low-count genes to improve DESeq2 estimation
keep <- rowSums(expr_mat) >= 10
expr_mat <- expr_mat[keep, , drop = FALSE]
gene_annot <- gene_annot[match(rownames(expr_mat), gene_annot$Geneid), , drop = FALSE]

cat("After low-count filter:", dim(expr_mat)[1], "genes x", dim(expr_mat)[2], "samples\n")

# ---------- DESeq helper ----------
run_deseq <- function(counts, coldata, grp1, grp2) {
  keep_samples <- coldata$condition %in% c(grp1, grp2)
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts[, keep_samples, drop = FALSE],
    colData   = droplevels(coldata[keep_samples, , drop = FALSE]),
    design    = ~ condition
  )
  
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("condition", grp1, grp2)))
  res$Geneid <- rownames(res)
  res <- res[!is.na(res$padj), , drop = FALSE]
  res
}

add_annotation <- function(res, gene_annot) {
  out <- dplyr::left_join(res, gene_annot, by = "Geneid")
  out$Genename <- clean_sym(out$Genename)
  out
}

# ---------- Run DE ----------
res_LMvsLP <- add_annotation(run_deseq(expr_mat, coldata, "LM", "LP"), gene_annot)
res_RMvsRP <- add_annotation(run_deseq(expr_mat, coldata, "RM", "RP"), gene_annot)
res_RPvsLP <- add_annotation(run_deseq(expr_mat, coldata, "RP", "LP"), gene_annot)
res_RMvsLM <- add_annotation(run_deseq(expr_mat, coldata, "RM", "LM"), gene_annot)

# ---------- DEG counts ----------
cnt_LMvsLP <- deg_count(res_LMvsLP, alpha, lfc_thr)
cnt_RMvsRP <- deg_count(res_RMvsRP, alpha, lfc_thr)
cnt_RPvsLP <- deg_count(res_RPvsLP, alpha, lfc_thr)
cnt_RMvsLM <- deg_count(res_RMvsLM, alpha, lfc_thr)

cat("=== DEG counts by contrast (padj < 0.05, |log2FC| >= 1) ===\n")
print(cbind(Contrast = "LMvsLP", cnt_LMvsLP))
print(cbind(Contrast = "RMvsRP", cnt_RMvsRP))
print(cbind(Contrast = "RPvsLP", cnt_RPvsLP))
print(cbind(Contrast = "RMvsLM", cnt_RMvsLM))

# ---------- Panel A: DEG counts ----------
deg_summary <- data.frame(
  Contrast = c("LM vs LP", "RM vs RP", "RP vs LP", "RM vs LM"),
  Up   = c(cnt_LMvsLP$Up, cnt_RMvsRP$Up, cnt_RPvsLP$Up, cnt_RMvsLM$Up),
  Down = c(cnt_LMvsLP$Down, cnt_RMvsRP$Down, cnt_RPvsLP$Down, cnt_RMvsLM$Down)
)

deg_long <- deg_summary %>%
  mutate(Down = -Down) %>%
  pivot_longer(cols = c("Up", "Down"), names_to = "Direction", values_to = "Count")

deg_long$Contrast <- factor(
  deg_long$Contrast,
  levels = c("LM vs LP", "RM vs RP", "RP vs LP", "RM vs LM")
)

pA <- ggplot(deg_long, aes(x = Contrast, y = Count, fill = Direction)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    x = NULL,
    y = "Number of DEGs"
  ) +
  theme_fig() +
  theme(legend.position = "top")

# ---------- Merge LV/RV for shared structure ----------
lv_tab <- res_LMvsLP %>%
  select(Geneid, Genename, log2FoldChange, padj) %>%
  rename(lfc_lv = log2FoldChange, padj_lv = padj, Genename_lv = Genename)

rv_tab <- res_RMvsRP %>%
  select(Geneid, Genename, log2FoldChange, padj) %>%
  rename(lfc_rv = log2FoldChange, padj_rv = padj, Genename_rv = Genename)

m <- full_join(lv_tab, rv_tab, by = "Geneid") %>%
  mutate(
    Genename = dplyr::coalesce(Genename_lv, Genename_rv),
    sig_lv = !is.na(padj_lv) & padj_lv < alpha & abs(lfc_lv) >= lfc_thr,
    sig_rv = !is.na(padj_rv) & padj_rv < alpha & abs(lfc_rv) >= lfc_thr,
    same_dir = !is.na(lfc_lv) & !is.na(lfc_rv) & sign(lfc_lv) == sign(lfc_rv),
    set = case_when(
      sig_lv & sig_rv & same_dir  ~ "Shared",
      sig_lv & !sig_rv            ~ "LV-only",
      sig_rv & !sig_lv            ~ "RV-only",
      sig_lv & sig_rv & !same_dir ~ "Opposite-dir",
      TRUE                        ~ "NS"
    )
  )

m$Genename <- clean_sym(m$Genename)
m$set <- factor(m$set, levels = c("NS", "Shared", "LV-only", "RV-only", "Opposite-dir"))

# ---------- Panel B: Correlation scatter ----------
# keep original display logic but compute correlation from observed values only
cor_df <- m %>%
  mutate(
    lfc_lv_plot = ifelse(is.na(lfc_lv), 0, lfc_lv),
    lfc_rv_plot = ifelse(is.na(lfc_rv), 0, lfc_rv)
  )

cor_val <- cor(m$lfc_lv, m$lfc_rv, use = "pairwise.complete.obs", method = "pearson")

label_df <- cor_df %>%
  filter(set != "NS") %>%
  mutate(rank_score = abs(lfc_lv_plot) + abs(lfc_rv_plot)) %>%
  arrange(desc(rank_score)) %>%
  slice_head(n = 12)

pB <- ggplot(cor_df, aes(x = lfc_lv_plot, y = lfc_rv_plot, color = set)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_hline(yintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_point(size = 1.5, alpha = 0.72) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  geom_text_repel(
    data = label_df,
    aes(label = Genename),
    size = 3,
    color = "black",
    box.padding = 0.25,
    point.padding = 0.15,
    max.overlaps = Inf
  ) +
  annotate(
    "text",
    x = min(cor_df$lfc_lv_plot, na.rm = TRUE),
    y = max(cor_df$lfc_rv_plot, na.rm = TRUE),
    label = paste0("r = ", round(cor_val, 2)),
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  scale_color_manual(values = my_set_cols, drop = FALSE) +
  labs(
    x = "log2FC (LM vs LP)",
    y = "log2FC (RM vs RP)"
  ) +
  theme_fig() +
  theme(legend.position = "right")

# ---------- Panel C: Set composition ----------
set_counts <- m %>%
  filter(set != "NS") %>%
  count(set) %>%
  mutate(set = factor(set, levels = c("LV-only", "RV-only", "Shared", "Opposite-dir")))

pC <- ggplot(set_counts, aes(x = set, y = n, fill = set)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  geom_text(aes(label = n), vjust = -0.35, size = 4) +
  scale_fill_manual(values = my_set_cols[names(my_set_cols) %in% levels(set_counts$set)]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    x = NULL,
    y = "Gene number"
  ) +
  theme_fig() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

# ---------- Panel D: Heatmap of top shared genes ----------
shared_same <- m %>%
  filter(set == "Shared") %>%
  mutate(rank_score = abs(lfc_lv) + abs(lfc_rv)) %>%
  arrange(desc(rank_score))

top_shared_ids <- shared_same %>%
  slice_head(n = n_heat) %>%
  pull(Geneid) %>%
  unique()

# vst for heatmap
dds_all <- DESeqDataSetFromMatrix(
  countData = expr_mat,
  colData   = coldata,
  design    = ~ condition
)
dds_all <- estimateSizeFactors(dds_all)
vst_mat <- assay(vst(dds_all, blind = TRUE))

heat_mat <- vst_mat[rownames(vst_mat) %in% top_shared_ids, , drop = FALSE]
heat_mat <- heat_mat[match(top_shared_ids, rownames(heat_mat)), , drop = FALSE]

heat_gene_names <- gene_annot$Genename[match(rownames(heat_mat), gene_annot$Geneid)]
heat_gene_names[is.na(heat_gene_names)] <- rownames(heat_mat)[is.na(heat_gene_names)]
heat_gene_names <- make.unique(heat_gene_names)
rownames(heat_mat) <- heat_gene_names

heat_mat_z <- t(scale(t(heat_mat)))
heat_mat_z[!is.finite(heat_mat_z)] <- 0

anno_col <- data.frame(Group = coldata[colnames(heat_mat_z), "condition", drop = TRUE])
rownames(anno_col) <- colnames(heat_mat_z)

ht <- pheatmap(
  heat_mat_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = anno_col,
  annotation_colors = list(Group = my_group_cols),
  color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100),
  border_color = NA,
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  silent = TRUE
)

pD <- ggplotify::as.ggplot(ht$gtable)

# ---------- Compose Figure 1 ----------
fig1 <- (pA + pB) / (pC + pD) +
  plot_layout(widths = c(1, 1.15), heights = c(1, 1.2)) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

# ---------- Save figure ----------
ggsave(file.path(outdir, "Figure1_final.pdf"), fig1, width = 13, height = 11)
ggsave(file.path(outdir, "Figure1_final.png"), fig1, width = 13, height = 11, dpi = 300)

# ---------- Export source data ----------
out_xlsx <- file.path(outdir, "Figure1_source_data.xlsx")
wb <- createWorkbook()

addWorksheet(wb, "LMvsLP")
writeData(wb, "LMvsLP", res_LMvsLP[, c("Geneid", "Genename", "baseMean", "log2FoldChange", "padj")])

addWorksheet(wb, "RMvsRP")
writeData(wb, "RMvsRP", res_RMvsRP[, c("Geneid", "Genename", "baseMean", "log2FoldChange", "padj")])

addWorksheet(wb, "RPvsLP")
writeData(wb, "RPvsLP", res_RPvsLP[, c("Geneid", "Genename", "baseMean", "log2FoldChange", "padj")])

addWorksheet(wb, "RMvsLM")
writeData(wb, "RMvsLM", res_RMvsLM[, c("Geneid", "Genename", "baseMean", "log2FoldChange", "padj")])

addWorksheet(wb, "LV_RV_merge")
writeData(wb, "LV_RV_merge", m)

addWorksheet(wb, "DEG_summary")
writeData(wb, "DEG_summary", deg_summary)

addWorksheet(wb, "Set_counts")
writeData(wb, "Set_counts", set_counts)

addWorksheet(wb, "RV_only_genes")
writeData(
  wb, "RV_only_genes",
  m %>%
    filter(set == "RV-only") %>%
    select(Geneid, Genename, lfc_lv, padj_lv, lfc_rv, padj_rv, set) %>%
    arrange(desc(abs(lfc_rv)), Genename)
)

addWorksheet(wb, "LV_only_genes")
writeData(
  wb, "LV_only_genes",
  m %>%
    filter(set == "LV-only") %>%
    select(Geneid, Genename, lfc_lv, padj_lv, lfc_rv, padj_rv, set) %>%
    arrange(desc(abs(lfc_lv)), Genename)
)

addWorksheet(wb, "Shared_genes")
writeData(
  wb, "Shared_genes",
  m %>%
    filter(set == "Shared") %>%
    mutate(rank_score = abs(lfc_lv) + abs(lfc_rv)) %>%
    select(Geneid, Genename, lfc_lv, padj_lv, lfc_rv, padj_rv, set, rank_score) %>%
    arrange(desc(rank_score), Genename)
)

addWorksheet(wb, "Top_shared_genes")
writeData(
  wb, "Top_shared_genes",
  shared_same %>%
    slice_head(n = n_heat) %>%
    select(Geneid, Genename, lfc_lv, padj_lv, lfc_rv, padj_rv, rank_score)
)

addWorksheet(wb, "Meta")
writeData(
  wb, "Meta",
  data.frame(
    Parameter = c("alpha", "lfc_thr", "n_heat", "low_count_filter"),
    Value = c(alpha, lfc_thr, n_heat, "rowSums(counts) >= 10")
  )
)

addWorksheet(wb, "sessionInfo")
writeData(wb, "sessionInfo", capture.output(sessionInfo()))

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

# ---------- Copy source data to data/processed ----------
file.copy(
  from = out_xlsx,
  to   = file.path(processed_dir, "Figure1_source_data.xlsx"),
  overwrite = TRUE
)

message("✅ Figure 1 finished. Outputs saved to: ", outdir)

# =========================
# Export individual panels
# =========================
ggsave(
  file.path(outdir, "Fig1A_DEG_counts.pdf"),
  pA,
  width = 3,
  height = 3,
  units = "in"
)

ggsave(
  file.path(outdir, "Fig1B_correlation.pdf"),
  pB,
  width = 5.5,
  height = 3.2,
  units = "in"
)

ggsave(
  file.path(outdir, "Fig1C_set_composition.pdf"),
  pC,
  width = 2.5,
  height = 2.5,
  units = "in"
)

ggsave(
  file.path(outdir, "Fig1D_heatmap.pdf"),
  pD,
  width = 7,
  height = 3.5,
  units = "in"
)

