# =========================================================
# Figure 2 — Functional enrichment for Shared and LV-only DEGs
# Parameters and plotting style fully matched to original Fig2 code
#
# Main figure panels:
#   A. Shared DEGs – GO BP dotplot
#   B. LV-only DEGs – GO BP dotplot
#   C. Shared DEGs – GO BP cnet
#   D. LV-only DEGs – GO BP cnet
#
# Notes:
#   - RV-only enrichment is RUN but NOT plotted in main figure
#   - RV-only GO results are exported to Figure2_source_data.xlsx
#
# Input:
#   outputs/Fig1_outputs/Figure1_source_data.xlsx
#
# Output:
#   outputs/Fig2_outputs/
#     - Figure2_final.pdf
#     - Figure2_final.png
#     - Figure2_source_data.xlsx
#     - Fig2A_Shared_DEGs_GO_BP_dot.pdf
#     - Fig2B_LVonly_DEGs_GO_BP_dot.pdf
#     - Fig2C_Shared_DEGs_GO_BP_cnet.pdf
#     - Fig2D_LVonly_DEGs_GO_BP_cnet.pdf
# =========================================================

options(stringsAsFactors = FALSE)
set.seed(1)

# ---------- 0) Packages ----------
suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(clusterProfiler)
  library(org.Rn.eg.db)
  library(enrichplot)
  library(GOSemSim)
})

# ---------- 1) Helper: locate script directory ----------
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

script_dir   <- get_script_dir()
project_dir  <- normalizePath(file.path(script_dir, ".."))
output_root  <- file.path(project_dir, "outputs")

if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)

# ---------- 2) Paths ----------
fig1_file <- file.path(output_root, "Fig1_outputs", "Figure1_source_data.xlsx")
out_dir   <- file.path(output_root, "Fig2_outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

if (!file.exists(fig1_file)) {
  stop("Required input not found: ", fig1_file, "\nPlease run Fig1 first.")
}

# ---------- 3) Original style parameters ----------
top_n_go     <- 12
show_n_cnet  <- 10
y_wrap       <- 30
pt_range     <- c(1.7, 5.0)
pal4         <- c("#440153", "#30678e", "#34b779", "#fce724")
base_size    <- 10
fig2_w       <- 12.5
fig2_h       <- 9.0

# enrichment thresholds
pval_cut <- 0.05
qval_cut <- 0.05
min_gs   <- 5
max_gs   <- 500

options(timeout = 600)

# ---------- 4) Read gene sets from Fig1 ----------
clean_vec <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[is.na(x) | x == "" | x == "NA"] <- NA_character_
  unique(stats::na.omit(x))
}

read_gene_sheet <- function(file, sheet_name) {
  df <- readxl::read_excel(file, sheet = sheet_name, col_names = TRUE)
  gene_col <- intersect(c("Genename", "Gene", "SYMBOL", "GeneSymbol", "Gene.name"), names(df))
  if (length(gene_col) == 0) {
    stop("No gene-name column found in sheet: ", sheet_name)
  }
  clean_vec(df[[gene_col[1]]])
}

genes_shared <- read_gene_sheet(fig1_file, "Shared_genes")
genes_LVonly <- read_gene_sheet(fig1_file, "LV_only_genes")
genes_RVonly <- read_gene_sheet(fig1_file, "RV_only_genes")

message(sprintf("Input gene counts — Shared=%d | LV-only=%d | RV-only=%d",
                length(genes_shared), length(genes_LVonly), length(genes_RVonly)))

# ---------- 5) ID mapping to ENTREZ (rat) ----------
guess_id_type <- function(x) {
  x <- na.omit(x)
  if (!length(x)) return(NA_character_)
  x0 <- head(x, min(length(x), 1000))
  x0n <- sub("\\..*$", "", x0)
  frac_ensembl <- mean(grepl("^ENS[A-Za-z]*", x0n))
  frac_entrez  <- mean(grepl("^[0-9]+$", x0n))
  if (is.na(frac_ensembl)) frac_ensembl <- 0
  if (is.na(frac_entrez))  frac_entrez  <- 0
  if (frac_ensembl >= 0.6) return("ENSEMBL")
  if (frac_entrez  >= 0.6) return("ENTREZID")
  "SYMBOL"
}

sym2entrez <- function(gvec) {
  gvec <- clean_vec(gvec)
  if (!length(gvec)) return(list(input_type = NA_character_, map = tibble(), entrez = character()))
  t <- guess_id_type(gvec)
  gclean <- if (t == "ENSEMBL") sub("\\..*$", "", gvec) else gvec
  message("Detected ID type: ", t, "  preview: ", paste(head(gvec, 3), collapse = ", "))
  if (t == "ENTREZID") {
    entrez <- unique(gclean)
    mapdf <- tibble::tibble(ENTREZID = entrez)
  } else {
    mapdf <- suppressMessages(
      clusterProfiler::bitr(gclean, fromType = t, toType = "ENTREZID", OrgDb = org.Rn.eg.db)
    )
    if (!nrow(mapdf)) stop("Mapping failed from ", t, ". Check species/ID type.")
    entrez <- unique(mapdf$ENTREZID)
  }
  list(input_type = t, map = tibble::as_tibble(mapdf), entrez = entrez)
}

m_shared <- sym2entrez(genes_shared)
m_LV     <- sym2entrez(genes_LVonly)
m_RV     <- sym2entrez(genes_RVonly)

message(sprintf("Mapped ENTREZ — Shared=%d | LV-only=%d | RV-only=%d",
                length(m_shared$entrez), length(m_LV$entrez), length(m_RV$entrez)))

# ---------- 6) GO enrichment ----------
run_go <- function(entrez_vec) {
  if (length(entrez_vec) < min_gs) return(NULL)
  clusterProfiler::enrichGO(
    gene         = entrez_vec,
    OrgDb        = org.Rn.eg.db,
    keyType      = "ENTREZID",
    ont          = "BP",
    pvalueCutoff = pval_cut,
    qvalueCutoff = qval_cut,
    pAdjustMethod = "BH",
    minGSSize    = min_gs,
    maxGSSize    = max_gs,
    readable     = TRUE
  )
}

ego_shared <- run_go(m_shared$entrez)
ego_LV     <- run_go(m_LV$entrez)
ego_RV     <- run_go(m_RV$entrez)

# ---------- 7) simplify + cnet（only category labels） ----------
simplify_safe <- function(ego) {
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  suppressMessages(
    clusterProfiler::simplify(
      ego,
      cutoff = 0.5,
      by = "p.adjust",
      select_fun = min,
      measure = "Wang"
    )
  )
}

ego_shared_s <- simplify_safe(ego_shared)
ego_LV_s     <- simplify_safe(ego_LV)
ego_RV_s     <- simplify_safe(ego_RV)

make_cnet_terms_only <- function(obj, title, base_size = 10, show_n = 12) {
  if (is.null(obj)) {
    return(
      ggplot() +
        annotate("text", x = .5, y = .5,
                 label = paste0(title, "\nNo significant terms"),
                 size = base_size / 3) +
        theme_void(base_size = base_size)
    )
  }
  
  if ("node_label" %in% names(formals(enrichplot::cnetplot))) {
    p <- enrichplot::cnetplot(
      obj,
      showCategory = show_n,
      circular = FALSE,
      colorEdge = TRUE,
      node_label = "category"
    )
  } else {
    p0 <- enrichplot::cnetplot(
      obj,
      showCategory = show_n,
      circular = FALSE,
      colorEdge = TRUE
    )
    p0$layers <- Filter(function(ly) {
      !inherits(ly$geom, "GeomText") && !inherits(ly$geom, "GeomLabel")
    }, p0$layers)
    p <- p0
  }
  
  p +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = .5, face = "bold", size = base_size),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2)
    )
}

# ---------- 8) dotplot ----------
dot_data <- function(obj, top_n = 12) {
  if (is.null(obj)) return(NULL)
  df <- tryCatch(as.data.frame(obj), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  if (!is.numeric(df$GeneRatio)) {
    rr <- strsplit(as.character(df$GeneRatio), "/")
    df$GeneRatio <- vapply(rr, function(x) {
      if (length(x) == 2) as.numeric(x[1]) / as.numeric(x[2]) else NA_real_
    }, 0.0)
  }
  
  df <- df[order(df$p.adjust, -df$GeneRatio, -df$Count), ]
  df <- df[seq_len(min(nrow(df), top_n)), , drop = FALSE]
  df$Description <- factor(df$Description, levels = rev(df$Description))
  df$neglogFDR   <- -log10(df$p.adjust)
  df
}

make_dot <- function(df, title, left_margin = 1) {
  if (is.null(df) || !nrow(df)) {
    return(
      ggplot() +
        annotate("text", x = .5, y = .5,
                 label = paste0(title, "\nNo significant terms"),
                 size = base_size / 3) +
        theme_void(base_size = base_size)
    )
  }
  
  ggplot(df, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = Count, colour = neglogFDR)) +
    scale_color_gradientn(colors = pal4, name = "-log10(FDR)") +
    scale_size_continuous(range = pt_range, name = "Gene count") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = y_wrap)) +
    labs(x = "Gene ratio", y = NULL, title = title) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title   = element_text(hjust = .5, face = "bold", size = base_size),
      axis.title.x = element_text(face = "bold"),
      axis.text.y  = element_text(color = "black"),
      axis.text.x  = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      plot.margin  = margin(6, 8, 6, left_margin)
    ) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.06)))
}

# ---------- 9) Main Figure 2 (Shared + LV-only only) ----------
df2A <- dot_data(ego_shared, top_n = top_n_go)
df2B <- dot_data(ego_LV,     top_n = top_n_go)
dfRV <- dot_data(ego_RV,     top_n = top_n_go)   # export only

g2A <- make_dot(df2A, "Shared DEGs-GO BP", left_margin = 2)
g2B <- make_dot(df2B, "LV-only DEGs-GO BP", left_margin = 2)

g2C <- make_cnet_terms_only(
  ego_shared_s,
  "Shared DEGs-GO BP cnet",
  base_size = base_size,
  show_n = show_n_cnet
)

g2D <- make_cnet_terms_only(
  ego_LV_s,
  "LV-only DEGs-GO BP cnet",
  base_size = base_size,
  show_n = show_n_cnet
)

fig2 <- ((g2A | g2B) / (g2C | g2D)) +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(
  file.path(out_dir, "Figure2_final.pdf"),
  plot = fig2,
  width = fig2_w,
  height = fig2_h
)

ggsave(
  file.path(out_dir, "Figure2_final.png"),
  plot = fig2,
  width = fig2_w,
  height = fig2_h,
  dpi = 300
)

message("✅ Main Figure 2 saved in: ", out_dir)

# ============================
# Single-panel export
# ============================
dot_w  <- 4.8
dot_h  <- 4.5
cnet_w <- 6.2
cnet_h <- 6.2

ggsave(
  file.path(out_dir, "Fig2A_Shared_DEGs_GO_BP_dot.pdf"),
  plot = g2A,
  width = dot_w,
  height = dot_h
)

ggsave(
  file.path(out_dir, "Fig2B_LVonly_DEGs_GO_BP_dot.pdf"),
  plot = g2B,
  width = dot_w,
  height = dot_h
)

ggsave(
  file.path(out_dir, "Fig2C_Shared_DEGs_GO_BP_cnet.pdf"),
  plot = g2C,
  width = cnet_w,
  height = cnet_h
)

ggsave(
  file.path(out_dir, "Fig2D_LVonly_DEGs_GO_BP_cnet.pdf"),
  plot = g2D,
  width = cnet_w,
  height = cnet_h
)

# ---------- 10) Export source data ----------
safe_df <- function(x) {
  if (is.null(x)) return(data.frame())
  out <- tryCatch(as.data.frame(x), error = function(e) data.frame())
  if (is.null(out)) data.frame() else out
}

wb <- createWorkbook()

addWorksheet(wb, "Shared_genes")
writeData(wb, "Shared_genes", data.frame(Gene = genes_shared))

addWorksheet(wb, "LV_only_genes")
writeData(wb, "LV_only_genes", data.frame(Gene = genes_LVonly))

addWorksheet(wb, "RV_only_genes")
writeData(wb, "RV_only_genes", data.frame(Gene = genes_RVonly))

addWorksheet(wb, "Shared_GO_BP")
writeData(wb, "Shared_GO_BP", safe_df(ego_shared))

addWorksheet(wb, "LVonly_GO_BP")
writeData(wb, "LVonly_GO_BP", safe_df(ego_LV))

addWorksheet(wb, "RVonly_GO_BP")
writeData(wb, "RVonly_GO_BP", safe_df(ego_RV))

addWorksheet(wb, "Shared_GO_BP_simplified")
writeData(wb, "Shared_GO_BP_simplified", safe_df(ego_shared_s))

addWorksheet(wb, "LVonly_GO_BP_simplified")
writeData(wb, "LVonly_GO_BP_simplified", safe_df(ego_LV_s))

addWorksheet(wb, "RVonly_GO_BP_simplified")
writeData(wb, "RVonly_GO_BP_simplified", safe_df(ego_RV_s))

if (!is.null(df2A)) {
  addWorksheet(wb, "Fig2A_dot_data")
  writeData(wb, "Fig2A_dot_data", df2A)
}

if (!is.null(df2B)) {
  addWorksheet(wb, "Fig2B_dot_data")
  writeData(wb, "Fig2B_dot_data", df2B)
}

if (!is.null(dfRV)) {
  addWorksheet(wb, "RVonly_dot_data")
  writeData(wb, "RVonly_dot_data", dfRV)
}

addWorksheet(wb, "ID_mapping_summary")
writeData(
  wb, "ID_mapping_summary",
  data.frame(
    Set = c("Shared", "LV-only", "RV-only"),
    Input_gene_count = c(length(genes_shared), length(genes_LVonly), length(genes_RVonly)),
    Mapped_ENTREZ_count = c(length(m_shared$entrez), length(m_LV$entrez), length(m_RV$entrez))
  )
)

addWorksheet(wb, "Plot_note")
writeData(
  wb, "Plot_note",
  data.frame(
    Item = c(
      "Main_figure_panel_A",
      "Main_figure_panel_B",
      "Main_figure_panel_C",
      "Main_figure_panel_D",
      "RV_only_plot_status"
    ),
    Value = c(
      "Shared GO BP dotplot",
      "LV-only GO BP dotplot",
      "Shared GO BP cnet",
      "LV-only GO BP cnet",
      "RV-only results exported only; not plotted in main figure because term number is limited"
    )
  )
)

addWorksheet(wb, "Meta")
writeData(
  wb, "Meta",
  data.frame(
    Parameter = c(
      "Input_file",
      "top_n_go",
      "show_n_cnet",
      "pval_cut",
      "qval_cut",
      "min_gs",
      "max_gs",
      "y_wrap",
      "pt_range_min",
      "pt_range_max",
      "fig2_w",
      "fig2_h",
      "dot_w",
      "dot_h"
    ),
    Value = c(
      fig1_file,
      top_n_go,
      show_n_cnet,
      pval_cut,
      qval_cut,
      min_gs,
      max_gs,
      y_wrap,
      pt_range[1],
      pt_range[2],
      fig2_w,
      fig2_h,
      dot_w,
      dot_h
    )
  )
)

addWorksheet(wb, "sessionInfo")
writeData(wb, "sessionInfo", capture.output(sessionInfo()))

saveWorkbook(
  wb,
  file.path(out_dir, "Figure2_source_data.xlsx"),
  overwrite = TRUE
)

message("✅ Figure 2 source data saved in: ", out_dir)

