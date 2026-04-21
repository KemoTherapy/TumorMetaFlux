#### Initialize ----
rm(list = ls())

suppressPackageStartupMessages({
  library(here); library(readxl); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(pheatmap)
  library(scales); library(stringr); library(ggrepel); library(grid)
})

#### Paths ----
MAPPING_FILE  <- here("data", "processed", "pancancer_metabolomics",
                       "data", "MasterMapping_MetImmune_03_16_2022_release.csv")
FLUX_DIR      <- here("results", "flux_calculations")
SUBSYSTEM_DIR <- here("results", "subsystems")
CORR_DIR      <- here("results", "correlations")
FIG_DIR       <- here("results", "figures", "flux_labeled")
dir.create(file.path(FIG_DIR, "pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "png"), recursive = TRUE, showWarnings = FALSE)

#### Constants ----
HAS_NORMAL <- c("BRCA1", "COAD", "GBM", "PDAC", "PRAD", "ccRCC3", "ccRCC4")

CANCER_ORDER <- c("BRCA1", "BRCA2", "COAD", "DLBCL", "GBM", "HCC", "ICC",
                  "HurthleCC", "OV", "PDAC", "PRAD", "ccRCC1", "ccRCC2", "ccRCC3", "ccRCC4")

CT_COLORS <- c(
  BRCA1 = "#C0392B", BRCA2 = "#E8A598", COAD = "#2980B9", DLBCL = "#1ABC9C",
  GBM = "#34495E", HCC = "#E59866", ICC = "#76B7A0", HurthleCC = "#7D6C8A",
  OV = "#922B21", PDAC = "#6D4C41", PRAD = "#9E9D8E", ccRCC1 = "#1A5276",
  ccRCC2 = "#2471A3", ccRCC3 = "#5499C7", ccRCC4 = "#A9CCE3"
)

EXCLUDE_SUBS <- c("Transport reactions", "Exchange/demand reactions",
                  "Drug metabolism", "Xenobiotics metabolism", "Isolated", "")

LOG2FC_THRESH <- 0.5
FDR_THRESH    <- 0.05

TAG_THEME <- theme(plot.tag = element_text(size = 13, face = "bold",
                                            color = "grey10", margin = margin(b=2, r=2)))

#### Helpers ----
norm_stem <- function(x) {
  x <- tools::file_path_sans_ext(basename(x))
  gsub("_+", "_", tolower(gsub("[.\\-\\s]+", "_", x)))
}

fix_id <- function(x) str_replace_all(x, "-", ".")

load_flux <- function(fp) {
  if (is.null(fp) || !file.exists(fp)) return(NULL)
  cat("  Reading:", basename(fp), "...")
  m <- tryCatch(read.csv(fp, row.names = 1, check.names = FALSE),
                error = function(e) { cat(" ERROR:", e$message, "\n"); NULL })
  if (!is.null(m)) cat("[", nrow(m), "x", ncol(m), "]\n")
  as.matrix(m)
}

match_samp <- function(fmat, ct_map) {
  if (is.null(fmat)) return(NULL)
  ids <- fix_id(ct_map$RNAID)
  ok  <- ids %in% colnames(fmat)
  if (!any(ok)) return(NULL)
  list(flux = fmat[, ids[ok], drop = FALSE], mapping = ct_map[ok, , drop = FALSE])
}

cbrt_norm <- function(x) sign(x) * abs(x)^(1/3)

save_fig <- function(p, name, w = 10, h = 8) {
  for (ext in c("pdf", "png")) {
    fp <- file.path(FIG_DIR, ext, paste0(name, ".", ext))
    tryCatch({
      if (inherits(p, "ggplot") || inherits(p, "patchwork"))
        ggsave(fp, p, width = w, height = h,
               dpi = if (ext == "png") 300 else 72, bg = "white")
      else {
        if (ext == "pdf") pdf(fp, w, h, bg = "white")
        else png(fp, w, h, units = "in", res = 300, bg = "white")
        print(p); dev.off()
      }
    }, error = function(e) cat("  Save error:", e$message, "\n"))
  }
  cat("  Saved:", name, "\n")
}

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
  theme(
    text             = element_text(family = "sans"),
    axis.text        = element_text(color = "grey20", size = base_size - 1),
    axis.title       = element_text(color = "grey10", size = base_size, face = "plain"),
    axis.line        = element_line(color = "grey40", linewidth = 0.35),
    axis.ticks       = element_line(color = "grey40", linewidth = 0.35),
    plot.title       = element_text(color = "grey10", size = base_size + 1,
                                    face = "bold", hjust = 0, margin = margin(b = 3)),
    plot.subtitle    = element_text(color = "grey45", size = base_size - 1,
                                    hjust = 0, margin = margin(b = 6)),
    strip.background = element_blank(),
    strip.text       = element_text(color = "grey15", size = base_size - 1, face = "bold"),
    panel.grid.major = element_line(color = "grey94", linewidth = 0.28),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(color = "grey15", size = base_size - 1),
    legend.text      = element_text(color = "grey25", size = base_size - 2),
    legend.key.size  = unit(0.85, "lines"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
}

# Mean absolute flux per subsystem per sample (rows = subsystems, cols = samples)
agg_sub <- function(fmat) {
  common <- intersect(rownames(fmat), rxn_annot$reaction)
  if (!length(common)) return(NULL)
  fsub <- fmat[common, , drop = FALSE]
  ann  <- rxn_annot$subsystem[match(common, rxn_annot$reaction)]
  result <- vapply(unique(ann), function(s) {
    colMeans(abs(fsub[which(ann == s), , drop = FALSE]), na.rm = TRUE)
  }, numeric(ncol(fsub)))
  t(result)
}

#### Load mapping and flux files ----
cat("Loading master mapping...\n")
mapping <- read.csv(MAPPING_FILE, stringsAsFactors = FALSE)
cat(" ", nrow(mapping), "samples,", length(unique(mapping$Dataset)), "datasets\n\n")

all_flux    <- list.files(FLUX_DIR, pattern = "\\.csv$", full.names = TRUE)
flux_lookup <- setNames(all_flux, norm_stem(all_flux))

find_file <- function(rna_file, lkp) {
  if (is.na(rna_file) || rna_file == "") return(NULL)
  s    <- norm_stem(rna_file)
  cand <- paste0("flux_results_", s)
  if (cand %in% names(lkp)) return(lkp[[cand]])
  if (s    %in% names(lkp)) return(lkp[[s]])
  hits <- names(lkp)[str_detect(names(lkp), fixed(s))]
  if (length(hits)) return(lkp[[hits[1]]])
  NULL
}

dataset_meta <- mapping %>%
  group_by(Dataset) %>%
  summarise(RNAFile = first(RNAFile), RNA_raw = first(RNA_raw_File), .groups = "drop") %>%
  rowwise() %>%
  mutate(f_tpm = list(find_file(RNAFile, flux_lookup)),
         f_raw = list(find_file(RNA_raw,  flux_lookup))) %>%
  ungroup()

cat("Dataset → flux file:\n")
for (i in seq_len(nrow(dataset_meta))) {
  ft <- dataset_meta$f_tpm[[i]]
  cat(" ", dataset_meta$Dataset[i], "→",
      if (!is.null(ft)) basename(ft) else "NOT FOUND", "\n")
}

#### Load subsystem annotation ----
cat("\nLoading subsystem annotation...\n")
sub_files <- list.files(SUBSYSTEM_DIR, pattern = "subsystem_analysis_.*\\.xlsx",
                         full.names = TRUE)
rxn_annot <- if (length(sub_files) > 0) {
  read_excel(sub_files[1], sheet = "Reactions_Annotated")
} else {
  gem_file <- here("human_GEM_equations.csv")
  if (!file.exists(gem_file))
    stop("No subsystem annotation found. Run 7_Subsystem.R or place human_GEM_equations.csv in project root.")
  read.csv(gem_file) %>% select(reaction = ID, subsystem = SUBSYSTEM)
}
rxn_annot <- filter(rxn_annot, !is.na(subsystem), subsystem != "",
                    !subsystem %in% EXCLUDE_SUBS)
cat(" ", n_distinct(rxn_annot$subsystem), "subsystems\n\n")

#### Load flux data ----
cat("Loading flux data...\n")
file_cache   <- list()
dataset_flux <- list()

for (i in seq_len(nrow(dataset_meta))) {
  ct <- dataset_meta$Dataset[i]
  ft <- dataset_meta$f_tpm[[i]]
  fr <- dataset_meta$f_raw[[i]]
  cat("\n[", ct, "]\n")

  fk <- if (!is.null(ft)) ft else "__null__"
  if (!fk %in% names(file_cache)) file_cache[[fk]] <- load_flux(ft)
  fmat <- file_cache[[fk]]
  if (is.null(fmat)) { cat("  SKIPPED\n"); next }

  # Apply cubic root normalization if files are un-normalized (max|v| > 10)
  if (max(abs(fmat), na.rm = TRUE) > 10) {
    warning(ct, ": flux appears un-normalized, applying cbrt now. Run 5_FluxCalculations.R to fix source files.")
    fmat <- cbrt_norm(fmat)
    file_cache[[fk]] <- fmat
  }

  rk <- if (!is.null(fr) && !identical(fr, ft)) fr else "__same__"
  if (rk != "__same__" && !rk %in% names(file_cache)) {
    fraw_loaded <- load_flux(fr)
    if (!is.null(fraw_loaded) && max(abs(fraw_loaded), na.rm = TRUE) > 10)
      fraw_loaded <- cbrt_norm(fraw_loaded)
    file_cache[[rk]] <- fraw_loaded
  }
  fraw <- if (rk != "__same__") file_cache[[rk]] else NULL

  ct_map  <- mapping[mapping$Dataset == ct, ]
  matched <- match_samp(fmat, ct_map)
  if (is.null(matched)) { cat("  SKIPPED: no sample match\n"); next }

  raw_m <- if (!is.null(fraw)) match_samp(fraw, ct_map) else NULL
  t_idx <- matched$mapping$TN == "Tumor"
  n_idx <- matched$mapping$TN == "Normal"
  cat("  Matched:", sum(t_idx), "tumor,", sum(n_idx), "normal\n")

  dataset_flux[[ct]] <- list(
    all_flux    = matched$flux,
    tumor_flux  = matched$flux[,  t_idx, drop = FALSE],
    normal_flux = if (sum(n_idx) >= 5) matched$flux[, n_idx, drop = FALSE] else NULL,
    flux_raw    = if (!is.null(raw_m)) raw_m$flux else NULL,
    mapping     = matched$mapping,
    has_normal  = sum(n_idx) >= 5
  )
}

cat("\nLoaded", length(dataset_flux), "datasets\n\n")
if (!length(dataset_flux)) stop("No datasets loaded. Check FLUX_DIR.")

cat("Aggregating to subsystem level...\n")
sub_list <- Filter(Negate(is.null),
  lapply(dataset_flux, function(d)
    tryCatch(agg_sub(d$all_flux), error = function(e) NULL)))
cat("  Done:", length(sub_list), "datasets\n\n")

#### Figure 3A — Pan-cancer subsystem heatmap ----
cat("Figure 3A: Pan-cancer heatmap...\n")
tryCatch({
  mean_sub <- lapply(sub_list, function(m) rowMeans(m, na.rm = TRUE))
  all_subs <- Reduce(intersect, lapply(mean_sub, names))
  hmat <- do.call(cbind, lapply(mean_sub, function(x) x[all_subs]))
  rownames(hmat) <- all_subs
  colnames(hmat) <- names(mean_sub)
  hmat <- hmat[, intersect(CANCER_ORDER, colnames(hmat)), drop = FALSE]

  rv   <- apply(hmat, 1, var, na.rm = TRUE)
  hmat <- hmat[names(sort(rv, decreasing = TRUE))[1:min(40, sum(rv > 0))],, drop = FALSE]
  hsc  <- t(scale(t(hmat)))
  hsc[!is.finite(hsc)] <- 0

  col_ann <- data.frame(Cohort = ifelse(colnames(hsc) %in% HAS_NORMAL,
                                         "T+N paired", "Tumor only"),
                        row.names = colnames(hsc))
  ann_col <- list(Cohort = c("T+N paired" = "#34495E", "Tumor only" = "#C0392B"))

  for (ext in c("pdf", "png")) {
    fp <- file.path(FIG_DIR, ext, paste0("Figure1A_pancancer_heatmap.", ext))
    if (ext == "pdf") pdf(fp, 13, 11, bg = "white")
    else              png(fp, 13, 11, units = "in", res = 300, bg = "white")
    ph <- pheatmap(hsc,
      color          = colorRampPalette(c("#3B6FA0", "#F5F5F5", "#A63025"))(100),
      cluster_rows   = TRUE, cluster_cols = TRUE,
      annotation_col = col_ann, annotation_colors = ann_col,
      labels_row     = str_trunc(rownames(hsc), 45),
      fontsize_row   = 7.5, fontsize_col = 9.5, fontsize = 9,
      border_color   = NA, na_col = "grey88",
      treeheight_row = 20, treeheight_col = 18, angle_col = 45,
      main = "A   Metabolic subsystem flux — pan-cancer comparison (row z-score)",
      silent = TRUE)
    grid.newpage(); grid.draw(ph$gtable); dev.off()
  }
  cat("  Saved: Figure1A_pancancer_heatmap\n")
}, error = function(e) cat("  ERROR:", e$message, "\n"))

#### Differential flux analysis (T vs N) ----
# Welch t-test per reaction, BH FDR within cohort
# Thresholds: FDR < 0.05, |log2FC| > 0.5
# log2FC = log2((mean_T + 0.01) / (mean_N + 0.01))
cat("Running T vs N differential analysis...\n")
tn_cts <- intersect(CANCER_ORDER,
  names(dataset_flux)[sapply(dataset_flux, `[[`, "has_normal")])
cat("  Cohorts:", paste(tn_cts, collapse = ", "), "\n")

da_results_named <- setNames(
  Filter(Negate(is.null), lapply(tn_cts, function(ct) {
    d <- dataset_flux[[ct]]
    if (is.null(d$normal_flux) || ncol(d$tumor_flux) < 3) return(NULL)
    n <- nrow(d$tumor_flux)
    pv <- lfc <- numeric(n)
    for (i in seq_len(n)) {
      tv <- d$tumor_flux[i, ]; nv <- d$normal_flux[i, ]
      pv[i]  <- tryCatch(t.test(tv, nv)$p.value, error = function(e) 1)
      lfc[i] <- log2((mean(tv, na.rm = TRUE) + 0.01) / (mean(nv, na.rm = TRUE) + 0.01))
    }
    padj <- p.adjust(pv, "BH")
    data.frame(reaction = rownames(d$tumor_flux), log2fc = lfc, padj = padj,
               up   = padj < FDR_THRESH & lfc >  LOG2FC_THRESH,
               down = padj < FDR_THRESH & lfc < -LOG2FC_THRESH,
               sig  = padj < FDR_THRESH, Dataset = ct, stringsAsFactors = FALSE)
  })),
  tn_cts[sapply(tn_cts, function(ct) {
    d <- dataset_flux[[ct]]
    !is.null(d$normal_flux) && ncol(d$tumor_flux) >= 3
  })]
)
cat("  Done:", length(da_results_named), "cohorts\n\n")

# DA and DF scores per subsystem
da_all <- bind_rows(da_results_named) %>%
  left_join(rxn_annot, by = "reaction") %>%
  filter(!is.na(subsystem), !subsystem %in% EXCLUDE_SUBS)

scores <- da_all %>%
  group_by(subsystem, Dataset) %>%
  summarise(n_total = n(),
            DA = (sum(up) - sum(down)) / n(),
            DF = sum(sig) / n(),
            .groups = "drop")

keep_s <- scores %>%
  group_by(subsystem) %>%
  filter(n_distinct(Dataset) >= min(3, length(da_results_named))) %>%
  pull(subsystem) %>% unique()
scores <- filter(scores, subsystem %in% keep_s)

da_w <- scores %>% select(subsystem, Dataset, DA) %>%
  pivot_wider(names_from = Dataset, values_from = DA, values_fill = 0)
da_m <- as.matrix(da_w[, -1]); rownames(da_m) <- da_w$subsystem
da_m <- da_m[, intersect(CANCER_ORDER, colnames(da_m)), drop = FALSE]
da_m <- da_m[order(abs(rowMeans(da_m, na.rm = TRUE)), decreasing = TRUE),, drop = FALSE]
if (nrow(da_m) > 50) da_m <- da_m[1:50, ]

df_w <- scores %>% select(subsystem, Dataset, DF) %>%
  pivot_wider(names_from = Dataset, values_from = DF, values_fill = 0)
df_m <- as.matrix(df_w[, -1]); rownames(df_m) <- df_w$subsystem
df_m <- df_m[rownames(da_m), intersect(CANCER_ORDER, colnames(df_m)), drop = FALSE]

mk_long <- function(mat, val)
  as.data.frame(mat) %>% tibble::rownames_to_column("sub") %>%
  pivot_longer(-sub, names_to = "Dataset", values_to = val) %>%
  mutate(Dataset = factor(Dataset, levels = colnames(mat)),
         sub     = factor(sub, levels = rev(rownames(mat))))
pd2 <- left_join(mk_long(da_m, "DA"), mk_long(df_m, "DF"), by = c("sub", "Dataset"))

#### Figure 6A — DA/DF dot matrix ----
cat("Figure 6A: DA/DF dot matrix...\n")
tryCatch({
  fig_da <- ggplot(pd2, aes(Dataset, sub, color = DA, size = DF)) +
    geom_point() +
    scale_color_gradient2(low = "#3B6FA0", mid = "#F5F5F5", high = "#A63025",
                          midpoint = 0, limits = c(-1, 1), name = "Flux\nDA score") +
    scale_size_continuous(range = c(0.3, 5), name = "DF score\n(frac. sig.)",
                          breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    scale_x_discrete(position = "top") +
    labs(tag = "A", x = NULL, y = NULL,
         title = "Metabolic flux changes: tumor vs. normal",
         subtitle = paste0("DA = (n\u2191\u2212n\u2193)/n reactions  \u00b7  ",
                           "dot size = fraction significantly changed (BH FDR<0.05, |log\u2082FC|>0.5)")) +
    theme_pub(9) + TAG_THEME +
    theme(axis.text.x  = element_text(angle = 40, hjust = 0, size = 8.5),
          axis.text.y  = element_text(size = 6.8),
          axis.line    = element_blank(), axis.ticks = element_blank(),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.25))
  save_fig(fig_da, "Figure2A_TvN_DA_score_matrix",
           w = max(7, ncol(da_m) * 1.15 + 2.5),
           h = max(9, nrow(da_m) * 0.20 + 3))
}, error = function(e) cat("  ERROR:", e$message, "\n"))

#### Supplemental — DA score correlation across cohorts ----
cat("Supplemental: DA correlation...\n")
tryCatch({
  if (ncol(da_m) >= 2) {
    dc <- cor(da_m, use = "pairwise.complete.obs", method = "spearman")
    for (ext in c("pdf", "png")) {
      fp <- file.path(FIG_DIR, ext, paste0("Figure3A_DA_correlation.", ext))
      if (ext == "pdf") pdf(fp, 5.5, 5, bg = "white")
      else              png(fp, 5.5, 5, units = "in", res = 300, bg = "white")
      ph <- pheatmap(dc,
        color           = colorRampPalette(c("#3B6FA0", "#F5F5F5", "#A63025"))(100),
        breaks          = seq(-1, 1, length.out = 101),
        display_numbers = TRUE, number_format = "%.2f",
        fontsize_number = 7.5, fontsize_row = 9, fontsize_col = 9,
        main            = "A   Spearman correlation of DA scores across T+N cohorts",
        border_color    = "white", treeheight_row = 12, treeheight_col = 12,
        silent = TRUE)
      grid.newpage(); grid.draw(ph$gtable); dev.off()
    }
    cat("  Saved: Figure3A_DA_correlation\n")
  }
}, error = function(e) cat("  ERROR:", e$message, "\n"))

#### Figure 4A-B — PCA ----
cat("Figure 4A-B: PCA...\n")
tryCatch({
  pieces <- Filter(Negate(is.null), lapply(names(sub_list), function(ct) {
    sm  <- sub_list[[ct]]; d <- dataset_flux[[ct]]
    if (is.null(sm) || nrow(sm) == 0) return(NULL)
    t_sm <- t(sm)
    lkp  <- setNames(d$mapping$TN, fix_id(d$mapping$RNAID))
    list(mat = t_sm, subs = rownames(sm), CT = rep(ct, nrow(t_sm)), TN = lkp[rownames(t_sm)])
  }))

  cs   <- Reduce(intersect, lapply(pieces, `[[`, "subs"))
  mats <- lapply(pieces, function(p) { idx <- match(cs, p$subs); p$mat[, idx, drop = FALSE] })
  big  <- do.call(rbind, mats)
  big[!is.finite(big)] <- 0
  big  <- big[, apply(big, 2, var) > 0]
  ct_v <- unlist(lapply(pieces, `[[`, "CT"))
  tn_v <- unlist(lapply(pieces, `[[`, "TN"))

  cat("  PCA:", nrow(big), "samples ×", ncol(big), "subsystems\n")
  pr <- prcomp(big, scale. = TRUE, center = TRUE)
  pv <- round(100 * pr$sdev^2 / sum(pr$sdev^2), 1)
  sc <- as.data.frame(pr$x[, 1:3])
  sc$CT <- factor(ct_v, levels = CANCER_ORDER); sc$TN <- tn_v

  p4a <- ggplot(sc, aes(PC1, PC2, color = CT, shape = TN)) +
    geom_point(alpha = 0.7, size = 1.4, stroke = 0.3) +
    scale_color_manual(values = CT_COLORS, name = NULL) +
    scale_shape_manual(values = c("Tumor" = 16, "Normal" = 1), name = NULL) +
    labs(tag = "A", x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"),
         title = "Cancer type") +
    theme_pub() + TAG_THEME +
    guides(color = guide_legend(ncol = 1, keyheight = 0.75,
                                override.aes = list(size = 2.5, alpha = 1)))

  p4b <- ggplot(sc, aes(PC1, PC2, color = TN)) +
    geom_point(alpha = 0.65, size = 1.3, stroke = 0.3) +
    scale_color_manual(values = c("Tumor" = "#A63025", "Normal" = "#3B6FA0"), name = NULL) +
    labs(tag = "B", x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"),
         title = "Tissue type") +
    theme_pub() + TAG_THEME

  fig4 <- (p4a + p4b) +
    plot_annotation(
      title    = "PCA of metabolic flux profiles",
      subtitle = paste0(nrow(sc), " samples \u00b7 ", ncol(big), " subsystems"),
      theme    = theme(plot.title    = element_text(size = 12, face = "bold", color = "grey10"),
                       plot.subtitle = element_text(size = 9,  color = "grey45"),
                       plot.background = element_rect(fill = "white", color = NA)))
  save_fig(fig4, "Figure4AB_PCA_flux", w = 13, h = 5.5)
}, error = function(e) cat("  ERROR:", e$message, "\n"))

#### Figure 5A-H — TPM vs RAW concordance ----
cat("Figure 5A-H: TPM vs RAW concordance...\n")
tryCatch({
  conc_cts <- intersect(CANCER_ORDER,
    names(dataset_flux)[sapply(names(dataset_flux),
      function(ct) !is.null(dataset_flux[[ct]]$flux_raw))])
  cat("  RNA-seq cohorts with both inputs:", paste(conc_cts, collapse = ", "), "\n")

  if (length(conc_cts) >= 2) {
    cps <- Filter(Negate(is.null), lapply(seq_along(conc_cts), function(k) {
      ct <- conc_cts[k]; d <- dataset_flux[[ct]]
      mt <- rowMeans(d$all_flux, na.rm = TRUE)
      mr <- rowMeans(d$flux_raw, na.rm = TRUE)
      co <- intersect(names(mt), names(mr))
      if (length(co) < 100) return(NULL)
      df <- data.frame(TPM = mt[co], RAW = mr[co])
      rv <- round(cor(df$TPM, df$RAW, use = "complete.obs"), 3)
      ggplot(df, aes(TPM, RAW)) +
        geom_point(alpha = 0.10, size = 0.35, color = CT_COLORS[ct]) +
        geom_abline(slope = 1, intercept = 0, color = "grey30",
                    linetype = "dashed", linewidth = 0.45) +
        annotate("text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.4,
                 label = paste0("r = ", rv), size = 2.8, color = "grey20") +
        labs(tag = LETTERS[k], title = ct,
             x = "Mean flux (TPM input)", y = "Mean flux (raw input)") +
        theme_pub(9) + TAG_THEME +
        theme(plot.title = element_text(color = CT_COLORS[ct], face = "bold", size = 9))
    }))
    if (length(cps)) {
      fig5 <- wrap_plots(cps, ncol = 4) +
        plot_annotation(
          title    = "Concordance between TPM-derived and raw-input flux predictions",
          subtitle = "Each point = one reaction (mean across samples)  \u00b7  dashed line = perfect concordance",
          theme    = theme(plot.title    = element_text(size = 11, face = "bold", color = "grey10"),
                           plot.subtitle = element_text(size = 8.5, color = "grey45"),
                           plot.background = element_rect(fill = "white", color = NA)))
      save_fig(fig5, "Figure5AH_TPM_RAW_concordance", w = 14, h = ceiling(length(cps) / 4) * 3.5)
    }
  }
}, error = function(e) cat("  ERROR:", e$message, "\n"))

#### Figure 7A-H — Subsystem violin plots ----
cat("Figure 7A-H: Subsystem violins...\n")
tryCatch({
  KEY_SUBS <- c(
    "Arginine and proline metabolism",
    "Fatty acid biosynthesis",
    "Fatty acid oxidation",
    "Glycolysis / Gluconeogenesis",
    "Nucleotide metabolism",
    "Oxidative phosphorylation",
    "Pentose phosphate pathway",
    "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism"
  )

  valid_tn <- tn_cts[sapply(tn_cts, function(ct) !is.null(sub_list[[ct]]))]
  vd_all <- bind_rows(Filter(Negate(is.null), lapply(valid_tn, function(ct) {
    sm <- sub_list[[ct]]; d <- dataset_flux[[ct]]
    av <- intersect(KEY_SUBS, rownames(sm))
    if (!length(av)) return(NULL)
    mat <- sm[match(av, rownames(sm)),, drop = FALSE]; rownames(mat) <- av
    df  <- as.data.frame(t(mat)); df$SID <- rownames(df)
    lkp <- setNames(d$mapping$TN, fix_id(d$mapping$RNAID))
    df$TN <- lkp[df$SID]; df$CT <- ct
    df %>% pivot_longer(all_of(av), names_to = "subsystem", values_to = "flux") %>%
      filter(!is.na(TN))
  })))
  vd_all$CT        <- factor(vd_all$CT, levels = CANCER_ORDER)
  vd_all$TN        <- factor(vd_all$TN, levels = c("Normal", "Tumor"))
  vd_all$subsystem <- factor(vd_all$subsystem, levels = KEY_SUBS)

  vplots <- Filter(Negate(is.null), lapply(seq_along(KEY_SUBS), function(k) {
    sub <- KEY_SUBS[k]
    d   <- filter(vd_all, subsystem == sub)
    if (nrow(d) == 0) return(NULL)
    short <- if (nchar(sub) > 35) paste0(substr(sub, 1, 32), "...") else sub
    ggplot(d, aes(CT, flux, fill = TN)) +
      geom_violin(scale = "width", alpha = 0.75,
                  position = position_dodge(0.82), linewidth = 0.15) +
      geom_boxplot(width = 0.06, outlier.shape = NA,
                   position = position_dodge(0.82),
                   fill = "white", color = "grey30", linewidth = 0.28) +
      scale_fill_manual(values = c("Tumor" = "#A63025", "Normal" = "#3B6FA0"), name = NULL) +
      labs(tag = LETTERS[k], x = NULL, y = "Mean norm. flux", title = short) +
      theme_pub(8.5) + TAG_THEME +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
            legend.position = "none",
            plot.title = element_text(size = 8, face = "bold"))
  }))

  fig6 <- wrap_plots(vplots, ncol = 3) +
    plot_annotation(
      title    = "Core metabolic pathway flux: tumor vs. normal",
      subtitle = "Cubic root-normalized flux \u00b7 7 cohorts with matched normal \u00b7 red = Tumor, blue = Normal",
      theme    = theme(plot.title    = element_text(size = 11, face = "bold", color = "grey10"),
                       plot.subtitle = element_text(size = 9,  color = "grey45"),
                       plot.background = element_rect(fill = "white", color = NA)))
  save_fig(fig6, "Figure6AH_subsystem_violin_TvN", w = 13, h = 12)
}, error = function(e) cat("  ERROR:", e$message, "\n"))

#### Figure 8A-G — Reaction-level volcano plots ----
cat("Figure 8A-G: Volcano plots...\n")
tryCatch({
  if (!length(da_results_named)) stop("No DA results")
  volc <- bind_rows(da_results_named) %>%
    mutate(negP = pmin(-log10(pmax(padj, 1e-30)), 30),
           cat  = case_when(up ~ "Up", down ~ "Down", TRUE ~ "ns"))

  cohort_order <- intersect(CANCER_ORDER, names(da_results_named))
  vps <- lapply(seq_along(cohort_order), function(k) {
    ct   <- cohort_order[k]
    d    <- volc[volc$Dataset == ct, ]
    n_u  <- sum(d$up,   na.rm = TRUE)
    n_d  <- sum(d$down, na.rm = TRUE)
    xrng <- max(abs(d$log2fc), na.rm = TRUE) * 0.72
    ymax <- max(d$negP, na.rm = TRUE)
    ggplot(d, aes(log2fc, negP, color = cat)) +
      geom_point(data = d[d$cat == "ns",],  alpha = 0.18, size = 0.38, show.legend = FALSE) +
      geom_point(data = d[d$cat != "ns",],  alpha = 0.80, size = 0.65, show.legend = FALSE) +
      geom_vline(xintercept = c(-LOG2FC_THRESH, LOG2FC_THRESH),
                 linetype = "dashed", color = "grey60", linewidth = 0.28) +
      geom_hline(yintercept = -log10(FDR_THRESH),
                 linetype = "dashed", color = "grey60", linewidth = 0.28) +
      scale_color_manual(values = c("Up" = "#A63025", "Down" = "#3B6FA0", "ns" = "grey78")) +
      annotate("text", x =  xrng, y = ymax * 0.93,
               label = paste0("\u2191", n_u), color = "#A63025", size = 2.6, fontface = "bold") +
      annotate("text", x = -xrng, y = ymax * 0.93,
               label = paste0("\u2193", n_d), color = "#3B6FA0", size = 2.6, fontface = "bold") +
      labs(tag = LETTERS[k], title = ct,
           x = expression(log[2]*"FC (tumor/normal)"),
           y = expression(-log[10]*"(FDR)")) +
      theme_pub(9) + TAG_THEME +
      theme(plot.title = element_text(color = CT_COLORS[ct], face = "bold", size = 9))
  })

  fig7 <- wrap_plots(vps, ncol = 3) +
    plot_annotation(
      title    = "Reaction-level differential metabolic flux: tumor vs. normal",
      subtitle = paste0("Welch t-test, BH FDR  \u00b7  thresholds: |log\u2082FC|>",
                        LOG2FC_THRESH, ", FDR<", FDR_THRESH),
      theme    = theme(plot.title    = element_text(size = 12, face = "bold", color = "grey10"),
                       plot.subtitle = element_text(size = 8.5, color = "grey45"),
                       plot.background = element_rect(fill = "white", color = NA)))
  save_fig(fig7, "Figure7AG_reaction_volcanoes", w = 13, h = ceiling(length(vps) / 3) * 4.2)
}, error = function(e) cat("  ERROR:", e$message, "\n"))

#### Figure 9A-B — Pan-cancer metabolite hubs ----
cat("Figure 9A-B: Metabolite hubs...\n")
tryCatch({
  corr_files <- list.files(CORR_DIR, pattern = "correlations_.*\\.xlsx", full.names = TRUE)
  if (!length(corr_files)) {
    cat("  SKIPPED: no correlation files found. Run 6_Correlations.R first.\n")
  } else {
    cat("  Found", length(corr_files), "correlation files\n")
    sig_all <- bind_rows(lapply(corr_files, function(f) {
      ct <- sub("correlations_", "", sub("\\.xlsx", "", basename(f)))
      tryCatch({
        df <- read_excel(f, sheet = "Significant_FDR005")
        if (nrow(df) == 0) return(NULL)
        df$CancerType <- ct; df
      }, error = function(e) NULL)
    }))

    if (is.null(sig_all) || nrow(sig_all) == 0) {
      cat("  SKIPPED: no significant pairs found\n")
    } else {
      sig_all$CancerType <- factor(sig_all$CancerType, levels = CANCER_ORDER)
      cat("  Significant pairs:", nrow(sig_all), "\n")

      hub <- sig_all %>%
        group_by(metabolite) %>%
        summarise(n_cohorts = n_distinct(CancerType),
                  max_abs_r = max(abs(correlation), na.rm = TRUE),
                  .groups   = "drop") %>%
        filter(n_cohorts >= 2) %>%
        arrange(desc(n_cohorts), desc(max_abs_r))
      top_hub <- head(hub, 30)

      pc3a <- ggplot(top_hub, aes(x = reorder(metabolite, n_cohorts),
                                   y = n_cohorts, fill = max_abs_r)) +
        geom_col(width = 0.72) +
        geom_text(aes(label = n_cohorts), hjust = -0.15, size = 2.7, color = "grey20") +
        scale_fill_gradient(low = "#F5E6CC", high = "#7D2B0A", name = "Max |r|", limits = c(0, 1)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = 1:15, labels = 1:15) +
        coord_flip() +
        labs(tag = "A", x = NULL, y = "Number of cancer cohorts",
             title = "Pan-cancer metabolite hubs",
             subtitle = "Metabolites with significant flux correlations in \u22652 cohorts") +
        theme_pub(9) + TAG_THEME

      hub_mat_df <- sig_all %>%
        filter(metabolite %in% top_hub$metabolite) %>%
        group_by(metabolite, CancerType) %>%
        summarise(max_r = max(abs(correlation), na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = CancerType, values_from = max_r, values_fill = 0) %>%
        arrange(desc(rowSums(select(., -metabolite) > 0)))

      hub_mat <- as.matrix(hub_mat_df[, -1])
      rownames(hub_mat) <- hub_mat_df$metabolite
      hub_mat <- hub_mat[, intersect(CANCER_ORDER, colnames(hub_mat)), drop = FALSE]

      pc3b_df <- as.data.frame(hub_mat) %>%
        tibble::rownames_to_column("metabolite") %>%
        pivot_longer(-metabolite, names_to = "CancerType", values_to = "max_r") %>%
        mutate(CancerType = factor(CancerType, levels = CANCER_ORDER),
               metabolite = factor(metabolite, levels = rev(rownames(hub_mat))))
      pc3b_df$max_r[pc3b_df$max_r == 0] <- NA

      pc3b <- ggplot(pc3b_df, aes(CancerType, metabolite, fill = max_r)) +
        geom_tile(color = "white", linewidth = 0.35) +
        scale_fill_gradient(low = "#F5E6CC", high = "#7D2B0A",
                            name = "Max |r|", na.value = "white") +
        scale_x_discrete(position = "top") +
        labs(tag = "B", x = NULL, y = NULL,
             title = "Cross-cohort metabolite significance map",
             subtitle = "Color intensity = maximum |r| among significant pairs in that cohort") +
        theme_pub(8.5) + TAG_THEME +
        theme(axis.text.x = element_text(angle = 40, hjust = 0, size = 7.5),
              axis.text.y = element_text(size = 7.5),
              axis.line   = element_blank(), axis.ticks = element_blank())

      fig_c3 <- pc3a + pc3b +
        plot_annotation(
          title = "Pan-cancer metabolite hub analysis",
          theme = theme(plot.title      = element_text(size = 12, face = "bold", color = "grey10"),
                        plot.background = element_rect(fill = "white", color = NA)))
      save_fig(fig_c3, "FigureC3AB_metabolite_hubs", w = 16, h = 9)
    }
  }
}, error = function(e) cat("  ERROR:", e$message, "\n"))

#### Write summary CSV ----
summ <- bind_rows(lapply(names(dataset_flux), function(ct) {
  d  <- dataset_flux[[ct]]
  da <- da_results_named[[ct]]
  data.frame(
    CancerType  = ct,
    N_Tumor     = ncol(d$tumor_flux),
    N_Normal    = if (!is.null(d$normal_flux)) ncol(d$normal_flux) else 0L,
    Has_RAW     = !is.null(d$flux_raw),
    N_Up        = if (!is.null(da)) sum(da$up,   na.rm = TRUE) else NA_integer_,
    N_Down      = if (!is.null(da)) sum(da$down, na.rm = TRUE) else NA_integer_,
    stringsAsFactors = FALSE
  )
}))
write.csv(summ, file.path(FIG_DIR, "flux_summary.csv"), row.names = FALSE)

cat("\nDone. Output:", FIG_DIR, "\n")
