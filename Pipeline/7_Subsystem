#### Initialize ----
rm(list = ls())

suppressPackageStartupMessages({
  library(here); library(tidyverse); library(readxl)
  library(openxlsx); library(pheatmap); library(RColorBrewer)
})

#### Configuration ----
HMR_FILE    <- here("human_GEM_equations.csv")
MAPPING_FILE <- here("data", "processed", "pancancer_metabolomics",
                     "data", "MasterMapping_MetImmune_03_16_2022_release.csv")
FLUX_DIR    <- here("results", "flux_calculations")

# outputDir MUST match SUBSYSTEM_DIR in 04_all_figures.R
OUTPUT_DIR  <- here("results", "subsystems")

# plots go under figures/ alongside all other outputs
PLOT_DIR    <- here("results", "figures", "subsystem_plots")
LOG_FILE    <- file.path(OUTPUT_DIR, "subsystem_log.txt")

# Subsystems excluded from analysis (non-metabolic / too generic)
EXCLUDE_SUBS <- c("Transport reactions", "Exchange/demand reactions",
                   "Drug metabolism", "Xenobiotics metabolism", "Isolated", "")

normaliseID <- function(x) gsub("-", ".", x, fixed=TRUE)

#### Validate Inputs ----
for (p in c(HMR_FILE, MAPPING_FILE, FLUX_DIR))
  if (!file.exists(p) && !dir.exists(p))
    stop("Required file not found: ", p)

dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(PLOT_DIR,   recursive=TRUE, showWarnings=FALSE)

#### Load HMR Annotations ----
hmr <- read.csv(HMR_FILE, stringsAsFactors=FALSE)
cat("=== Subsystem Analysis ===\n")
cat("HMR reactions loaded:", nrow(hmr), "\n")

rxnSub <- hmr %>%
  select(ID, SUBSYSTEM) %>%
  filter(!is.na(SUBSYSTEM), !SUBSYSTEM %in% EXCLUDE_SUBS) %>%
  rename(reaction=ID, subsystem=SUBSYSTEM)

cat("Reactions with valid subsystem annotation:", nrow(rxnSub), "\n")
cat("Unique subsystems:", n_distinct(rxnSub$subsystem), "\n\n")

#### Load Mapping ----
mapping     <- read.csv(MAPPING_FILE, stringsAsFactors=FALSE)
cancerTypes <- unique(mapping$Dataset)
cat("Cancer types:", length(cancerTypes), "\n\n")

#### Logging ----
logLines <- c(paste("=== Subsystem Log —", Sys.time(), "==="),
              paste("Subsystems:", n_distinct(rxnSub$subsystem)), "")
writeLines(logLines, LOG_FILE)
lg <- function(...) { line <- paste(...); cat(line, "\n"); write(line, LOG_FILE, append=TRUE) }

#### Subsystem Analysis Function ----
analyzeSubsystems <- function(ct) {

  lg(strrep("=", 42))
  lg("Analyzing:", ct)
  lg(strrep("=", 42))

  ctMap      <- mapping[mapping$Dataset == ct, ]
  rnaBase    <- tools::file_path_sans_ext(unique(ctMap$RNAFile)[1])
  fluxFile   <- file.path(FLUX_DIR, paste0("flux_results_", rnaBase, ".csv"))
  outputFile <- file.path(OUTPUT_DIR, paste0("subsystem_analysis_", ct, ".xlsx"))

  if (!file.exists(fluxFile)) { lg("  SKIP: flux file not found"); return(NULL) }
  if (file.exists(outputFile)){ lg("  SKIP: output exists");       return(NULL) }

  tryCatch({

    fluxData <- read.csv(fluxFile, row.names=1, check.names=FALSE)
    lg("  Flux:", nrow(fluxData), "reactions x", ncol(fluxData), "samples")

    # Annotate reactions with subsystem membership
    annotated <- data.frame(reaction=rownames(fluxData), stringsAsFactors=FALSE) %>%
      left_join(rxnSub, by="reaction")

    nMapped <- sum(!is.na(annotated$subsystem))
    lg("  Mapped:", nMapped, "/", nrow(fluxData),
       "(", round(100*nMapped/nrow(fluxData), 1), "%)")

    if (nMapped < 100) { lg("  SKIP: too few mapped reactions"); return(NULL) }

    fluxAnnot <- fluxData
    fluxAnnot$subsystem <- annotated$subsystem
    fluxAnnot$reaction  <- rownames(fluxData)

    fluxSub <- fluxAnnot %>%
      filter(!is.na(subsystem)) %>%
      select(reaction, subsystem, everything())

    # Subsystem-level flux aggregation
    # IMPORTANT: uses abs(flux) to match agg_sub() in 04_all_figures.R
    # Subsystem activity = mean absolute flux; directionality analysed separately
    subFlux <- fluxSub %>%
      pivot_longer(cols=-c(reaction, subsystem),
                   names_to="sample", values_to="flux") %>%
      group_by(subsystem, sample) %>%
      summarise(
        mean_flux   = mean(abs(flux), na.rm=TRUE),
        median_flux = median(abs(flux), na.rm=TRUE),
        sd_flux     = sd(flux, na.rm=TRUE),
        n_reactions = n(),
        .groups     = "drop"
      )

    # Subsystem × sample matrix (for heatmap and T vs N comparison)
    subMat <- subFlux %>%
      select(subsystem, sample, mean_flux) %>%
      pivot_wider(names_from=sample, values_from=mean_flux) %>%
      column_to_rownames("subsystem") %>%
      as.matrix()

    lg("  Subsystem matrix:", nrow(subMat), "x", ncol(subMat))

    # Cross-sample variability
    subVar <- subFlux %>%
      group_by(subsystem) %>%
      summarise(
        mean_across_samples = mean(mean_flux, na.rm=TRUE),
        sd_across_samples   = sd(mean_flux, na.rm=TRUE),
        cv                  = sd_across_samples / abs(mean_across_samples),
        n_samples           = n(),
        .groups             = "drop"
      ) %>%
      arrange(desc(mean_across_samples))

    lg("  Top subsystem:", subVar$subsystem[1],
       "| mean flux =", round(subVar$mean_across_samples[1], 4))

    # Tumor vs Normal comparison (only if ≥5 samples per group)
    tnComp <- NULL
    tMap <- ctMap %>% filter(TN == "Tumor")
    nMap <- ctMap %>% filter(TN == "Normal")

    if (nrow(tMap) >= 5 && nrow(nMap) >= 5) {
      lg("  Running T vs N comparison...")

      tIDs <- intersect(colnames(subMat), normaliseID(unique(tMap$RNAID)))
      nIDs <- intersect(colnames(subMat), normaliseID(unique(nMap$RNAID)))
      lg("  T samples:", length(tIDs), "| N samples:", length(nIDs))

      if (length(tIDs) >= 5 && length(nIDs) >= 5) {
        tnComp <- data.frame(
          subsystem   = rownames(subMat),
          tumor_mean  = rowMeans(subMat[, tIDs, drop=FALSE], na.rm=TRUE),
          normal_mean = rowMeans(subMat[, nIDs, drop=FALSE], na.rm=TRUE),
          stringsAsFactors = FALSE
        )
        tnComp$log2FC <- log2((tnComp$tumor_mean  + 0.01) /
                               (tnComp$normal_mean + 0.01))
        tnComp$pValue <- vapply(seq_len(nrow(tnComp)), function(i) {
          tv <- subMat[tnComp$subsystem[i], tIDs]
          nv <- subMat[tnComp$subsystem[i], nIDs]
          if (sum(!is.na(tv)) < 3 || sum(!is.na(nv)) < 3) return(NA_real_)
          tryCatch(t.test(tv, nv)$p.value, error=function(e) NA_real_)
        }, numeric(1))
        tnComp$pValueAdjusted <- p.adjust(tnComp$pValue, method="BH")
        tnComp <- tnComp %>% arrange(pValueAdjusted) %>%
          mutate(significant = pValueAdjusted < 0.05)
        lg("  Sig subsystems (FDR<0.05):", sum(tnComp$significant, na.rm=TRUE))
      }
    }

    # Save Excel workbook
    wb <- createWorkbook()

    addWorksheet(wb, "Subsystem_Summary"); writeData(wb, "Subsystem_Summary", subVar)

    subFluxWide <- subFlux %>% select(subsystem, sample, mean_flux) %>%
      pivot_wider(names_from=sample, values_from=mean_flux)
    addWorksheet(wb, "Subsystem_by_Sample"); writeData(wb, "Subsystem_by_Sample", subFluxWide)

    if (!is.null(tnComp)) {
      addWorksheet(wb, "Tumor_vs_Normal"); writeData(wb, "Tumor_vs_Normal", tnComp)
    }

    # Reactions_Annotated sheet — required by 04_all_figures.R to load subsystem membership
    addWorksheet(wb, "Reactions_Annotated")
    writeData(wb, "Reactions_Annotated",
              fluxSub %>% select(reaction, subsystem) %>% distinct())

    saveWorkbook(wb, outputFile, overwrite=TRUE)
    lg("  Saved:", basename(outputFile))

    # Plots
    if (nrow(subMat) > 2 && ncol(subMat) > 2) {
      topSubs <- subVar %>% filter(!is.infinite(cv)) %>%
        arrange(desc(cv)) %>% head(30) %>% pull(subsystem)

      if (length(topSubs) > 5) {
        scaled <- t(scale(t(subMat[topSubs, , drop=FALSE])))
        pf1 <- file.path(PLOT_DIR, paste0("heatmap_subsystems_", ct, ".png"))
        png(pf1, width=12, height=10, units="in", res=300)
        pheatmap(scaled,
                 color=colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
                 main=paste("Subsystem Activity —", ct),
                 clustering_distance_rows="correlation",
                 clustering_distance_cols="correlation",
                 show_colnames=FALSE, fontsize_row=8,
                 breaks=seq(-3, 3, length.out=101))
        dev.off()
        lg("  Saved heatmap")
      }
    }

    top20 <- subVar %>% arrange(desc(mean_across_samples)) %>% head(20)
    p2 <- ggplot(top20, aes(reorder(subsystem, mean_across_samples),
                             mean_across_samples)) +
      geom_bar(stat="identity", fill="steelblue") + coord_flip() +
      labs(title=paste("Top 20 Subsystems —", ct), x=NULL, y="Mean |Flux|") +
      theme_minimal() + theme(axis.text.y=element_text(size=9))
    ggsave(file.path(PLOT_DIR, paste0("barplot_subsystems_", ct, ".png")),
           p2, width=10, height=8, dpi=300)

    if (!is.null(tnComp)) {
      tnPlot <- tnComp %>% filter(!is.na(pValueAdjusted), !is.infinite(log2FC))
      if (nrow(tnPlot) > 10) {
        p3 <- ggplot(tnPlot, aes(log2FC, -log10(pValue))) +
          geom_point(aes(color=significant), alpha=0.6, size=2) +
          scale_color_manual(values=c("grey70","red")) +
          geom_text(data=tnPlot %>% filter(significant) %>%
                         arrange(pValueAdjusted) %>% head(10),
                    aes(label=subsystem), size=2.5, hjust=0, vjust=0,
                    check_overlap=TRUE) +
          geom_vline(xintercept=0,        linetype="dashed", color="grey40") +
          geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
          labs(title=paste("T vs N Subsystem Flux —", ct),
               x=expression(log[2]*"FC"),
               y=expression(-log[10]*"(p)"),
               color="FDR<0.05") +
          theme_minimal() + theme(legend.position="top")
        ggsave(file.path(PLOT_DIR, paste0("volcano_TvN_", ct, ".png")),
               p3, width=10, height=8, dpi=300)
      }
    }

    lg("  OK:", ct)
    return(list(ct=ct, nSub=nrow(subMat), nSamp=ncol(subMat),
                topSub=subVar$subsystem[1],
                topFlux=round(subVar$mean_across_samples[1], 4)))

  }, error=function(e) { lg("  ERROR:", ct, "--", e$message); return(NULL) })
}

#### Run All Cancer Types ----
cat("Processing all cancer types...\n\n")
results <- Filter(Negate(is.null), lapply(cancerTypes, analyzeSubsystems))
names(results) <- sapply(results, `[[`, "ct")

#### Pan-Cancer Summary ----
lg(strrep("=", 42))
lg("Pan-cancer subsystem comparison")
lg(strrep("=", 42))

allData <- bind_rows(Filter(Negate(is.null), lapply(names(results), function(ct) {
  f <- file.path(OUTPUT_DIR, paste0("subsystem_analysis_", ct, ".xlsx"))
  if (!file.exists(f)) return(NULL)
  d <- read.xlsx(f, sheet="Subsystem_Summary")
  d$cancerType <- ct; d
})))

if (nrow(allData) > 0) {
  compMat <- allData %>%
    select(subsystem, cancerType, mean_across_samples) %>%
    pivot_wider(names_from=cancerType, values_from=mean_across_samples) %>%
    column_to_rownames("subsystem") %>% as.matrix()

  coverage    <- rowSums(!is.na(compMat))
  commonSubs  <- names(coverage[coverage >= 5])

  if (length(commonSubs) > 10) {
    scaled <- t(scale(t(compMat[commonSubs, , drop=FALSE])))
    pf <- file.path(PLOT_DIR, "heatmap_pancancer_subsystems.png")
    png(pf, width=14, height=12, units="in", res=300)
    pheatmap(scaled,
             color=colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             main="Pan-Cancer Subsystem Activity",
             clustering_distance_rows="correlation",
             clustering_distance_cols="correlation",
             fontsize_row=7, fontsize_col=10,
             breaks=seq(-3, 3, length.out=101))
    dev.off()
    lg("Saved pan-cancer heatmap")
  }

  wb2 <- createWorkbook()
  addWorksheet(wb2, "All_Cancer_Types"); writeData(wb2, "All_Cancer_Types", allData)
  if (length(commonSubs) > 0) {
    cDF <- as.data.frame(compMat[commonSubs, , drop=FALSE])
    cDF$subsystem <- rownames(cDF)
    cDF <- cDF %>% select(subsystem, everything())
    addWorksheet(wb2, "Comparison_Matrix"); writeData(wb2, "Comparison_Matrix", cDF)
  }
  saveWorkbook(wb2, file.path(OUTPUT_DIR, "pancancer_subsystem_summary.xlsx"),
               overwrite=TRUE)
  lg("Saved pan-cancer summary")
}

#### Final Summary ----
cat("\n=== Done ===\n")
if (length(results) > 0) {
  sumDF <- bind_rows(results) %>% select(ct, nSub, nSamp, topSub, topFlux)
  print(sumDF); cat("\n")
}
cat("Results:", OUTPUT_DIR, "\nPlots:", PLOT_DIR, "\nLog:", LOG_FILE, "\n")
