#### Initialize ----
rm(list = ls())

suppressPackageStartupMessages({
  library(METAFlux)
  library(here)
})

#### Configuration ----
TRANSCRIPTOMICS_DIR <- here("data", "processed", "pancancer_metabolomics",
                             "data", "transcriptomics_processed")
OUTPUT_DIR          <- here("results", "flux_calculations")
LOG_FILE            <- file.path(OUTPUT_DIR, "flux_calculations_log.txt")
MEDIUM              <- "human_blood"   # human_blood | cell_medium

# Sign-preserving cubic root normalization — paper formula: sign(v) * |v|^(1/3)
cubicRootTransform <- function(x) sign(x) * abs(x)^(1/3)

#### Validate Inputs ----
if (!dir.exists(TRANSCRIPTOMICS_DIR))
  stop("Transcriptomics directory not found: ", TRANSCRIPTOMICS_DIR,
       "\nRun 4_PreprocessingData.R first.")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

transcriptomicsFiles <- list.files(
  TRANSCRIPTOMICS_DIR, pattern = "\\.csv$", full.names = TRUE)

if (length(transcriptomicsFiles) == 0)
  stop("No .csv files in: ", TRANSCRIPTOMICS_DIR,
       "\nRun 4_PreprocessingData.R first.")

cat("=== Flux Calculations ===\n")
cat("Input: ", TRANSCRIPTOMICS_DIR, "\n")
cat("Output:", OUTPUT_DIR, "\n")
cat("Files: ", length(transcriptomicsFiles), "\n")
cat("Medium:", MEDIUM, "\n\n")

#### Load Medium and Initialize Log ----
data(list = MEDIUM)
mediumData <- get(MEDIUM)

logLines <- c(paste("=== Flux Log —", Sys.time(), "==="),
              paste("Medium:", MEDIUM),
              paste("Files found:", length(transcriptomicsFiles)), "")
writeLines(logLines, LOG_FILE)
lg <- function(...) { line <- paste(...); cat(line, "\n"); write(line, LOG_FILE, append=TRUE) }

#### Process Each File ----
results <- list()

for (f in transcriptomicsFiles) {

  baseName   <- tools::file_path_sans_ext(basename(f))
  outputPath <- file.path(OUTPUT_DIR, paste0("flux_results_", baseName, ".csv"))

  if (file.exists(outputPath)) { lg("SKIP:", baseName, "(exists)"); next }

  lg("START:", baseName)
  t0 <- proc.time()

  tryCatch({
    # check.names=FALSE preserves sample ID hyphens for downstream matching
    expr <- read.csv(f, row.names=1, check.names=FALSE)
    if (nrow(expr)==0 || ncol(expr)==0) stop("Empty expression matrix")
    lg("  Dim:", nrow(expr), "genes x", ncol(expr), "samples")

    # MRAS: AND-nodes = min across subunits; OR-nodes = max across isoforms
    lg("  Computing MRAS...")
    mras <- calculate_reaction_score(expr)

    lg("  Running FBA...")
    rawFlux  <- compute_flux(mras=mras, medium=mediumData)
    normFlux <- cubicRootTransform(rawFlux)

    # Sanity check: cubic root of anything >1 is smaller in magnitude
    if (max(abs(normFlux), na.rm=TRUE) > max(abs(rawFlux), na.rm=TRUE))
      warning(baseName, ": normalized max > raw max — check input normalization state")

    write.csv(normFlux, outputPath, row.names=TRUE)
    elapsed <- round((proc.time()-t0)["elapsed"], 1)
    lg("  OK:", nrow(normFlux), "reactions |", elapsed, "s")

    results[[baseName]] <- list(ok=TRUE, nR=nrow(normFlux), nS=ncol(normFlux), t=elapsed)

  }, error=function(e) {
    lg("  ERROR:", baseName, "--", e$message)
    results[[baseName]] <<- list(ok=FALSE, error=e$message)
  })
}

#### Summary ----
nOK  <- sum(sapply(results, `[[`, "ok"), na.rm=TRUE)
nErr <- sum(!sapply(results, `[[`, "ok"), na.rm=TRUE)
nSkip <- length(transcriptomicsFiles) - length(results)
cat("\n=== Done: OK =", nOK, "| Errors =", nErr, "| Skipped =", nSkip, "===\n")
if (nErr > 0) for (nm in names(results))
  if (!results[[nm]]$ok) cat(" FAILED:", nm, "--", results[[nm]]$error, "\n")
cat("Results:", OUTPUT_DIR, "\nLog:", LOG_FILE, "\n")
