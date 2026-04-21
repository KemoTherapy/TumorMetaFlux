#### Initialize ----
rm(list = ls())

suppressPackageStartupMessages({
  library(here); library(readxl); library(openxlsx); library(ggplot2)
})

#### Configuration ----
MAPPING_FILE     <- here("data", "processed", "pancancer_metabolomics",
                          "data", "MasterMapping_MetImmune_03_16_2022_release.csv")
METABOLOMICS_DIR <- here("data", "processed", "pancancer_metabolomics",
                          "data", "metabolomics_processed")
FLUX_DIR         <- here("results", "flux_calculations")
OUTPUT_DIR       <- here("results", "correlations")
PLOT_DIR         <- here("results", "figures", "correlation_plots")
LOG_FILE         <- file.path(OUTPUT_DIR, "correlations_log.txt")
MIN_SAMPLES      <- 10     # minimum matched samples to test a pair
FDR_THRESH       <- 0.05   # BH FDR threshold

# Normalise sample IDs: master mapping uses hyphens; R may convert to dots
normaliseID <- function(x) gsub("-", ".", x, fixed=TRUE)

#### Validate Inputs ----
for (p in c(MAPPING_FILE, FLUX_DIR, METABOLOMICS_DIR))
  if (!file.exists(p) && !dir.exists(p))
    stop("Required path not found: ", p,
         "\nRun earlier pipeline scripts first.")

dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(PLOT_DIR,   recursive=TRUE, showWarnings=FALSE)

#### Load Mapping ----
mapping <- read.csv(MAPPING_FILE, stringsAsFactors=FALSE)
cancerTypes <- unique(mapping$Dataset)
cat("=== Flux-Metabolite Correlations ===\n")
cat("Samples:", nrow(mapping), "| Cancer types:", length(cancerTypes), "\n\n")

#### Logging ----
logLines <- c(paste("=== Correlation Log —", Sys.time(), "==="),
              paste("MIN_SAMPLES:", MIN_SAMPLES),
              paste("FDR_THRESH:", FDR_THRESH), "")
writeLines(logLines, LOG_FILE)
lg <- function(...) { line <- paste(...); cat(line, "\n"); write(line, LOG_FILE, append=TRUE) }

#### Core Correlation Function ----
correlateCancerType <- function(ct) {

  lg(strrep("=", 42))
  lg("Processing:", ct)
  lg(strrep("=", 42))

  ctMap       <- mapping[mapping$Dataset == ct, ]
  metabFile   <- file.path(METABOLOMICS_DIR, unique(ctMap$MetabFile)[1])
  rnaBase     <- tools::file_path_sans_ext(unique(ctMap$RNAFile)[1])
  fluxFile    <- file.path(FLUX_DIR, paste0("flux_results_", rnaBase, ".csv"))
  outputFile  <- file.path(OUTPUT_DIR, paste0("correlations_", ct, ".xlsx"))

  # Input validation
  if (!file.exists(fluxFile)) {
    lg("  SKIP: flux file not found:", fluxFile); return(NULL)
  }
  if (!file.exists(metabFile)) {
    lg("  SKIP: metabolomics file not found:", metabFile); return(NULL)
  }
  if (file.exists(outputFile)) {
    lg("  SKIP: output already exists"); return(NULL)
  }

  tryCatch({

    # Load metabolomics (metabolites × samples)
    metabRaw    <- read_excel(metabFile, sheet="metabo_imputed_filtered_Tumor")
    metabNames  <- metabRaw[[1]]
    metabData   <- as.data.frame(metabRaw[, -1])
    rownames(metabData) <- metabNames
    lg("  Metabolomics:", nrow(metabData), "metabolites x", ncol(metabData), "samples")

    # Load flux (reactions × samples) — check.names=FALSE preserves hyphens
    fluxData <- read.csv(fluxFile, row.names=1, check.names=FALSE)
    lg("  Flux:", nrow(fluxData), "reactions x", ncol(fluxData), "samples")

    # Sample matching via CommonID bridge
    # Normalise both sides to dot-delimited IDs before lookup
    m2c <- setNames(ctMap$CommonID, ctMap$MetabID)
    r2c <- setNames(ctMap$CommonID, normaliseID(ctMap$RNAID))

    metabCIDs <- m2c[colnames(metabData)]
    fluxCIDs  <- r2c[normaliseID(colnames(fluxData))]

    commonIDs <- intersect(
      metabCIDs[!is.na(metabCIDs)],
      fluxCIDs[ !is.na(fluxCIDs)]
    )
    lg("  Matched samples:", length(commonIDs))

    if (length(commonIDs) < MIN_SAMPLES) {
      lg("  SKIP: <", MIN_SAMPLES, "matched samples"); return(NULL)
    }

    metabIdx <- match(names(m2c)[m2c %in% commonIDs], colnames(metabData))
    fluxIdx  <- match(names(r2c)[r2c %in% commonIDs],
                      normaliseID(colnames(fluxData)))

    valid <- !is.na(metabIdx) & !is.na(fluxIdx)
    metabMatched <- metabData[,  metabIdx[valid], drop=FALSE]
    fluxMatched  <- fluxData[,   fluxIdx[valid],  drop=FALSE]

    lg("  Aligned dims — Metab:", nrow(metabMatched), "x", ncol(metabMatched),
       "| Flux:", nrow(fluxMatched), "x", ncol(fluxMatched))

    # Pearson correlations — preallocated vectors, no rbind in loop
    nM <- nrow(metabMatched); nR <- nrow(fluxMatched)
    maxN <- nM * nR
    lg("  Testing up to", format(maxN, big.mark=","), "pairs...")

    mVec <- character(maxN); rVec <- character(maxN)
    cVec <- numeric(maxN);   pVec <- numeric(maxN)
    idx  <- 1L

    for (i in seq_len(nM)) {
      mv <- as.numeric(metabMatched[i, ])
      for (j in seq_len(nR)) {
        fv   <- as.numeric(fluxMatched[j, ])
        ok   <- !is.na(mv) & !is.na(fv)
        if (sum(ok) < MIN_SAMPLES) next
        ct_  <- tryCatch(
          cor.test(mv[ok], fv[ok], method="pearson"),
          warning=function(w) NULL, error=function(e) NULL
        )
        if (!is.null(ct_) && !is.na(ct_$estimate)) {
          mVec[idx] <- rownames(metabMatched)[i]
          rVec[idx] <- rownames(fluxMatched)[j]
          cVec[idx] <- as.numeric(ct_$estimate)
          pVec[idx] <- ct_$p.value
          idx <- idx + 1L
        }
      }
      if (i %% 50 == 0)
        lg("  Progress:", i, "/", nM, "metabolites —",
           format(idx-1L, big.mark=","), "valid pairs")
    }

    # Trim and assemble
    n_     <- idx - 1L
    corDF  <- data.frame(
      metabolite  = mVec[1:n_],
      reaction    = rVec[1:n_],
      correlation = cVec[1:n_],
      pValue      = pVec[1:n_],
      stringsAsFactors = FALSE
    )
    lg("  Total valid pairs:", format(nrow(corDF), big.mark=","))

    # BH FDR correction within cohort
    corDF$pValueAdjusted <- p.adjust(corDF$pValue, method="BH")
    sigDF <- corDF[!is.na(corDF$pValueAdjusted) &
                     corDF$pValueAdjusted < FDR_THRESH, ]
    lg("  Significant (FDR <", FDR_THRESH, "):", nrow(sigDF))

    # Top 100 by |r|
    top100 <- corDF[order(-abs(corDF$correlation)), ][1:min(100, nrow(corDF)), ]

    # Summary sheet
    summSheet <- data.frame(
      Metric = c("Cancer Type", "Matched Samples", "Metabolites",
                 "Flux Reactions", "Total Correlations",
                 "Significant (FDR < 0.05)",
                 "Median Correlation", "Mean |Correlation|"),
      Value  = c(ct, ncol(metabMatched), nrow(metabMatched),
                 nrow(fluxMatched), nrow(corDF), nrow(sigDF),
                 round(median(corDF$correlation, na.rm=TRUE), 3),
                 round(mean(abs(corDF$correlation), na.rm=TRUE), 3)),
      stringsAsFactors = FALSE
    )

    # Write workbook
    wb <- createWorkbook()
    addWorksheet(wb, "Significant_FDR005")
    writeData(wb, "Significant_FDR005",
              sigDF[order(-abs(sigDF$correlation)), ])
    addWorksheet(wb, "Top100")
    writeData(wb, "Top100", top100)
    addWorksheet(wb, "Summary")
    writeData(wb, "Summary", summSheet)
    saveWorkbook(wb, outputFile, overwrite=TRUE)
    lg("  Saved:", basename(outputFile))

    # Volcano plot
    if (nrow(sigDF) > 0) {
      plotDF <- corDF
      plotDF$sig    <- plotDF$pValueAdjusted < FDR_THRESH
      plotDF$negP   <- -log10(plotDF$pValue)
      topMets       <- names(sort(table(sigDF$metabolite), decreasing=TRUE)[1:10])
      plotDF$label  <- ifelse(plotDF$metabolite %in% topMets, plotDF$metabolite, "")

      p <- ggplot(plotDF, aes(correlation, negP)) +
        geom_point(aes(color=sig), alpha=0.4, size=0.5) +
        scale_color_manual(values=c("grey70","red")) +
        geom_text(data=plotDF[plotDF$label!="",],
                  aes(label=label), size=2.5, hjust=0, vjust=0, check_overlap=TRUE) +
        labs(title=paste("Flux-Metabolite Correlations —", ct),
             x="Pearson r", y=expression(-log[10]*"(p)"), color="FDR<0.05") +
        theme_minimal() + theme(legend.position="top")
      ggsave(file.path(PLOT_DIR, paste0("correlation_plot_", ct, ".png")),
             p, width=10, height=8, dpi=300)
      lg("  Plot saved")
    }

    lg("  OK:", ct)
    return(list(ct=ct, n=ncol(metabMatched),
                nCorr=nrow(corDF), nSig=nrow(sigDF)))

  }, error=function(e) {
    lg("  ERROR:", ct, "--", e$message); return(NULL)
  })
}

#### Run All Cancer Types ----
cat("\nProcessing all cancer types...\n\n")
results <- Filter(Negate(is.null), lapply(cancerTypes, correlateCancerType))
names(results) <- sapply(results, `[[`, "ct")

#### Summary ----
cat("\n=== Summary ===\n")
if (length(results) > 0) {
  sumDF <- do.call(rbind, lapply(results, function(x)
    data.frame(CancerType=x$ct, Samples=x$n,
               Correlations=x$nCorr, Significant=x$nSig,
               stringsAsFactors=FALSE)))
  print(sumDF); cat("\n")
  wb2 <- createWorkbook()
  addWorksheet(wb2, "Summary"); writeData(wb2, "Summary", sumDF)
  saveWorkbook(wb2, file.path(OUTPUT_DIR, "analysis_summary.xlsx"), overwrite=TRUE)
}
cat("Results:", OUTPUT_DIR, "\nPlots:", PLOT_DIR, "\nLog:", LOG_FILE, "\n")
