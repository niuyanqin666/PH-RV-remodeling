## ==========================================================
## Figure 3 — LV-only / RV-only / Shared signature scores
## Gene sets are read from:
##   outputs/Fig1_outputs/Figure1_source_data.xlsx
##     - LV_only_genes
##     - RV_only_genes
##     - Shared_genes
##
## Input:
##   data/raw/GSE240923_processed-data-mct.xlsx
##   outputs/Fig1_outputs/Figure1_source_data.xlsx
##
## Output:
##   outputs/Fig3_outputs/
##     - Fig3A_Batch1_SignatureScores.pdf
##     - Fig3A_Batch1_SignatureScores.png
##     - Fig3B_Batch2_SignatureScores.pdf
##     - Fig3B_Batch2_SignatureScores.png
##     - Fig3_source_data.xlsx
## ==========================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(AnnotationDbi)
  library(org.Rn.eg.db)
  library(rstatix)
})

set.seed(123)

TOUP <- function(x) toupper(trimws(as.character(x)))

## ----------------------------------------------------------
## 0) 自动定位脚本目录
## ----------------------------------------------------------
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

if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)

## ----------------------------------------------------------
## 1) 路径
## ----------------------------------------------------------
infile_expr <- file.path(raw_dir, "GSE240923_processed-data-mct.xlsx")
fig1_xlsx   <- file.path(output_root, "Fig1_outputs", "Figure1_source_data.xlsx")
out_dir     <- file.path(output_root, "Fig3_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(infile_expr)) stop("Expression file not found: ", infile_expr)
if (!file.exists(fig1_xlsx))   stop("Figure1 source data file not found: ", fig1_xlsx)

## ----------------------------------------------------------
## 2) 公共主题与尺寸
## ----------------------------------------------------------
base_size <- 12
theme_fig3 <- theme_classic(base_size = base_size) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

fig_w <- 8.8
fig_h <- 3.6

## ----------------------------------------------------------
## 3) 读表达数据 + 行名转 SYMBOL
## ----------------------------------------------------------
meta   <- read.xlsx(infile_expr, sheet = 1)
count1 <- read.xlsx(infile_expr, sheet = 2, rowNames = TRUE)  # Batch1
count2 <- read.xlsx(infile_expr, sheet = 3, rowNames = TRUE)  # Batch2
count1 <- as.matrix(count1); mode(count1) <- "numeric"
count2 <- as.matrix(count2); mode(count2) <- "numeric"

ens2symbol_and_collapse <- function(mat_counts) {
  ens_ids <- rownames(mat_counts)
  ens_ids_clean <- sub("\\.\\d+$", "", ens_ids)
  
  symbol_map <- AnnotationDbi::mapIds(
    org.Rn.eg.db,
    keys = ens_ids_clean,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  map_df <- data.frame(
    ENSEMBL = ens_ids,
    SYMBOL = unname(symbol_map),
    stringsAsFactors = FALSE
  )
  
  df <- as.data.frame(mat_counts, check.names = FALSE)
  df$Gene <- ifelse(is.na(map_df$SYMBOL) | map_df$SYMBOL == "", map_df$ENSEMBL, map_df$SYMBOL)
  
  df_sym <- df %>%
    group_by(Gene) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    as.data.frame()
  
  rownames(df_sym) <- df_sym$Gene
  df_sym <- df_sym[, setdiff(colnames(df_sym), "Gene"), drop = FALSE]
  
  rn_up <- TOUP(rownames(df_sym))
  mat_sym_upper <- df_sym
  rownames(mat_sym_upper) <- rn_up
  
  if (any(duplicated(rownames(mat_sym_upper)))) {
    mat_sym_upper <- rowsum(mat_sym_upper, group = rownames(mat_sym_upper)) /
      as.vector(table(rownames(mat_sym_upper)))
  }
  
  list(
    mat_sym = df_sym,
    mat_sym_upper = mat_sym_upper,
    map_df = map_df
  )
}

res_b1 <- ens2symbol_and_collapse(count1)
res_b2 <- ens2symbol_and_collapse(count2)

count1_sym <- as.matrix(res_b1$mat_sym_upper); mode(count1_sym) <- "numeric"
count2_sym <- as.matrix(res_b2$mat_sym_upper); mode(count2_sym) <- "numeric"

## ----------------------------------------------------------
## 4) 构建元信息
## ----------------------------------------------------------
# Batch1
samples_b1 <- colnames(count1_sym)
meta_b1 <- data.frame(
  Sample = samples_b1,
  Group = dplyr::case_when(
    grepl("^CTRL",         samples_b1) ~ "Ctrl",
    grepl("^COMP",         samples_b1) ~ "cRV",
    grepl("^DECOMP-FP",    samples_b1) ~ "early_dRV",
    grepl("^DECOMP-MCT69", samples_b1) ~ "late_dRV",
    TRUE ~ NA_character_
  ),
  stringsAsFactors = FALSE
)
meta_b1$Group <- factor(meta_b1$Group, levels = c("Ctrl", "cRV", "early_dRV", "late_dRV"))
rownames(meta_b1) <- meta_b1$Sample

# Batch2
samples_b2 <- colnames(count2_sym)
meta_b2 <- data.frame(
  Sample = samples_b2,
  Group = dplyr::case_when(
    grepl("^CTRL",   samples_b2, ignore.case = TRUE) ~ "Ctrl",
    grepl("^COMP",   samples_b2, ignore.case = TRUE) ~ "cRV",
    grepl("^DECOMP", samples_b2, ignore.case = TRUE) ~ "dRV",
    TRUE ~ NA_character_
  ),
  Sex = dplyr::case_when(
    grepl("-f(_|$)", samples_b2, ignore.case = TRUE) ~ "F",
    grepl("-m(_|$)", samples_b2, ignore.case = TRUE) ~ "M",
    TRUE ~ NA_character_
  ),
  stringsAsFactors = FALSE
)
meta_b2$Group <- factor(meta_b2$Group, levels = c("Ctrl", "cRV", "dRV"))
meta_b2$Sex   <- factor(meta_b2$Sex, levels = c("F", "M"))
rownames(meta_b2) <- meta_b2$Sample

## ----------------------------------------------------------
## 5) 从 Figure1_source_data.xlsx 读取三类基因
## ----------------------------------------------------------
read_gene_set_from_fig1 <- function(path, sheet_name) {
  stopifnot(file.exists(path))
  df <- openxlsx::read.xlsx(path, sheet = sheet_name)
  
  gcol <- intersect(c("Genename", "SYMBOL", "Gene", "GeneSymbol", "Gene.name"), names(df))
  if (length(gcol) == 0) {
    stop("No gene-name column found in sheet: ", sheet_name)
  }
  
  genes <- TOUP(unique(na.omit(df[[gcol[1]]])))
  genes[genes != ""]
}

sets <- list(
  LV_only = read_gene_set_from_fig1(fig1_xlsx, "LV_only_genes"),
  RV_only = read_gene_set_from_fig1(fig1_xlsx, "RV_only_genes"),
  Shared  = read_gene_set_from_fig1(fig1_xlsx, "Shared_genes")
)

message(sprintf(
  "[Figure1_source_data] LV-only=%d, RV-only=%d, Shared=%d",
  length(sets$LV_only), length(sets$RV_only), length(sets$Shared)
))

## ----------------------------------------------------------
## 6) 计算 signature 分数
##    行向 z-score -> 每个样本对集合取均值
## ----------------------------------------------------------
score_sets <- function(expr_sym, meta_df, sets_list, group_levels, has_sex = FALSE) {
  stopifnot(all(colnames(expr_sym) == rownames(meta_df)))
  
  calc_one <- function(glist, label) {
    keep <- rownames(expr_sym) %in% TOUP(glist)
    if (!any(keep)) return(NULL)
    
    mat <- expr_sym[keep, , drop = FALSE]
    mat <- t(scale(t(as.matrix(mat))))
    mat[!is.finite(mat)] <- 0
    
    data.frame(
      Sample = colnames(mat),
      Score = colMeans(mat),
      Set = label,
      stringsAsFactors = FALSE
    )
  }
  
  df <- dplyr::bind_rows(
    calc_one(sets_list$LV_only, "LV-only"),
    calc_one(sets_list$RV_only, "RV-only"),
    calc_one(sets_list$Shared,  "Shared")
  )
  
  meta2 <- meta_df
  if (!"Sample" %in% names(meta2)) meta2$Sample <- rownames(meta2)
  
  df <- dplyr::left_join(df, meta2, by = "Sample")
  df$Group <- factor(df$Group, levels = group_levels)
  if (has_sex && "Sex" %in% names(df)) {
    df$Sex <- droplevels(factor(df$Sex, levels = c("F", "M")))
  }
  df$Set <- factor(df$Set, levels = c("LV-only", "RV-only", "Shared"))
  df
}

## ----------------------------------------------------------
## 7) 单图绘制：三个面板并列
## ----------------------------------------------------------
plot_scores_faceted <- function(df_scores, title, group_levels, outfile,
                                ctrl_label = "Ctrl", same_y_scale = TRUE,
                                width_in = 8.8, height_in = 3.6,
                                anova_x = 0.6) {
  
  comps <- lapply(setdiff(group_levels, ctrl_label), function(x) c(ctrl_label, x))
  
  y_stats <- df_scores %>%
    dplyr::group_by(Set) %>%
    dplyr::summarise(
      y_min = min(Score, na.rm = TRUE),
      y_max = max(Score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(y_rng = y_max - y_min)
  
  pair_df <- dplyr::bind_rows(lapply(levels(df_scores$Set), function(set_i) {
    dat <- df_scores[df_scores$Set == set_i, , drop = FALSE]
    if (nrow(dat) == 0) return(NULL)
    
    tt <- rstatix::t_test(dat, Score ~ Group, comparisons = comps) %>%
      rstatix::add_significance("p") %>%
      dplyr::mutate(Set = set_i)
    
    tt <- dplyr::left_join(tt, y_stats, by = "Set") %>%
      dplyr::group_by(Set) %>%
      dplyr::mutate(
        .idx = dplyr::row_number(),
        y.position = y_max + (.idx) * 0.08 * y_rng
      ) %>%
      dplyr::ungroup()
    
    dplyr::transmute(tt, Set, group1, group2, label = p.signif, y.position)
  }))
  
  anova_text_df <- dplyr::bind_rows(lapply(levels(df_scores$Set), function(set_i) {
    dat <- df_scores[df_scores$Set == set_i, , drop = FALSE]
    if (nrow(dat) == 0) return(NULL)
    
    pval <- rstatix::anova_test(dat, Score ~ Group)$p
    ys <- y_stats[y_stats$Set == set_i, ]
    n_pairs <- length(comps)
    y_pos <- ys$y_max + (n_pairs + 1) * 0.08 * ys$y_rng + 0.02 * ys$y_rng
    
    data.frame(
      Set = set_i,
      x = anova_x,
      y.position = y_pos,
      label = sprintf("p = %.2e", pval),
      stringsAsFactors = FALSE
    )
  }))
  
  y_top_global <- max(c(pair_df$y.position, anova_text_df$y.position), na.rm = TRUE)
  y_min_global <- min(df_scores$Score, na.rm = TRUE)
  y_rng_global <- diff(range(df_scores$Score, na.rm = TRUE))
  y_top_global <- y_top_global + 0.04 * y_rng_global
  
  gp <- ggplot(df_scores, aes(Group, Score, fill = Group, color = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.70, width = 0.65) +
    geom_jitter(width = 0.12, size = 1.8, alpha = 0.90) +
    facet_wrap(~Set, nrow = 1, scales = if (same_y_scale) "fixed" else "free_y") +
    theme_fig3 +
    theme(
      plot.margin = margin(t = 16, r = 12, b = 6, l = 10)
    ) +
    labs(title = title, x = NULL, y = "Z-score mean") +
    ggpubr::stat_pvalue_manual(
      pair_df,
      label = "label",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.6,
      hide.ns = TRUE
    ) +
    geom_text(
      data = anova_text_df,
      mapping = aes(x = x, y = y.position, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0,
      size = 4
    ) +
    coord_cartesian(
      ylim = c(y_min_global - 0.05 * y_rng_global, y_top_global),
      clip = "off"
    ) +
    guides(fill = "none")
  
  ggsave(outfile, gp, width = width_in, height = height_in)
  gp
}

## ----------------------------------------------------------
## 8) 生成 Batch1 / Batch2 图
## ----------------------------------------------------------
fig3_b1 <- score_sets(
  count1_sym, meta_b1, sets,
  group_levels = c("Ctrl", "cRV", "early_dRV", "late_dRV"),
  has_sex = FALSE
)

p3A <- plot_scores_faceted(
  fig3_b1,
  "Batch1 Signature Scores",
  group_levels = c("Ctrl", "cRV", "early_dRV", "late_dRV"),
  outfile = file.path(out_dir, "Fig3A_Batch1_SignatureScores.pdf"),
  same_y_scale = TRUE,
  width_in = fig_w+2,
  height_in = fig_h
)

fig3_b2 <- score_sets(
  count2_sym, meta_b2, sets,
  group_levels = c("Ctrl", "cRV", "dRV"),
  has_sex = TRUE
)

p3B <- plot_scores_faceted(
  fig3_b2,
  "Batch2 Signature Scores",
  group_levels = c("Ctrl", "cRV", "dRV"),
  outfile = file.path(out_dir, "Fig3B_Batch2_SignatureScores.pdf"),
  same_y_scale = TRUE,
  width_in = fig_w+2,
  height_in = fig_h
)

ggsave(file.path(out_dir, "Fig3A_Batch1_SignatureScores.png"), p3A, width = fig_w, height = fig_h, dpi = 300)
ggsave(file.path(out_dir, "Fig3B_Batch2_SignatureScores.png"), p3B, width = fig_w, height = fig_h, dpi = 300)

## ----------------------------------------------------------
## 9) 导出结果和基因列表
## ----------------------------------------------------------
wb <- createWorkbook()

addWorksheet(wb, "Batch1_scores")
writeData(wb, "Batch1_scores", fig3_b1)

addWorksheet(wb, "Batch2_scores")
writeData(wb, "Batch2_scores", fig3_b2)

addWorksheet(wb, "LV_only_genes")
writeData(wb, "LV_only_genes", data.frame(Gene = sets$LV_only))

addWorksheet(wb, "RV_only_genes")
writeData(wb, "RV_only_genes", data.frame(Gene = sets$RV_only))

addWorksheet(wb, "Shared_genes")
writeData(wb, "Shared_genes", data.frame(Gene = sets$Shared))

addWorksheet(wb, "Overlap_summary")
writeData(
  wb, "Overlap_summary",
  data.frame(
    Set = c("LV-only", "RV-only", "Shared"),
    N_in_Figure1 = c(length(sets$LV_only), length(sets$RV_only), length(sets$Shared)),
    N_in_Batch1 = c(
      sum(TOUP(sets$LV_only) %in% rownames(count1_sym)),
      sum(TOUP(sets$RV_only) %in% rownames(count1_sym)),
      sum(TOUP(sets$Shared)  %in% rownames(count1_sym))
    ),
    N_in_Batch2 = c(
      sum(TOUP(sets$LV_only) %in% rownames(count2_sym)),
      sum(TOUP(sets$RV_only) %in% rownames(count2_sym)),
      sum(TOUP(sets$Shared)  %in% rownames(count2_sym))
    )
  )
)

addWorksheet(wb, "Meta")
writeData(
  wb, "Meta",
  data.frame(
    Parameter = c(
      "Expression_input",
      "Figure1_input",
      "Batch1_groups",
      "Batch2_groups",
      "fig_w",
      "fig_h"
    ),
    Value = c(
      infile_expr,
      fig1_xlsx,
      "Ctrl,cRV,early_dRV,late_dRV",
      "Ctrl,cRV,dRV",
      fig_w,
      fig_h
    )
  )
)

addWorksheet(wb, "sessionInfo")
writeData(wb, "sessionInfo", capture.output(sessionInfo()))

saveWorkbook(wb, file.path(out_dir, "Fig3_source_data.xlsx"), overwrite = TRUE)

message("✅ Figure 3 finished. Outputs saved to: ", out_dir)