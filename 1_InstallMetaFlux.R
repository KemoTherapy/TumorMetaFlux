#### Initialize ----
rm(list = ls())

#### Install Prerequisites ----
# METAFlux requires processx, devtools, and Matrix from CRAN

if (!require("processx", quietly = TRUE)) {
  install.packages("processx")
}

if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

if (!require("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}

#### Clear Package Locks ----
# remove any existing package lock directories to prevent installation errors
lockDirectory <- file.path(Sys.getenv("R_LIBS_USER"), "00LOCK")

if (dir.exists(lockDirectory)) {
  unlink(lockDirectory, recursive = TRUE, force = TRUE)
  cat("Cleared package lock directory\n")
}

#### Install METAFlux from GitHub ----
# source: https://github.com/KChen-lab/METAFlux
library(devtools)

devtools::install_github("KChen-lab/METAFlux")

library(METAFlux)

cat("\n✓ METAFlux installation complete!\n")
