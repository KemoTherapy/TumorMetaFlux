#### Initialize ----
rm(list = ls())

#### Load Required Packages ----
library(here)
library(tidyverse)
library(openxlsx)
library(SummarizedExperiment)
library(reshape2)
library(maplet)

#### Set Working Directory to Preprocessing Folder ----
setwd(here("scripts", "preprocessing"))

#### Copy Data Folder ----
sourceDataDir <- here("data", "processed", "pancancer_metabolomics", "data")
targetDataDir <- file.path(getwd(), "data")

if (dir.exists(targetDataDir)) {
  unlink(targetDataDir, recursive = TRUE)
}

file.copy(sourceDataDir, getwd(), recursive = TRUE)

#### Comment Out setwd Lines in Preprocessing Scripts ----
for (scriptNum in 1:5) {
  scriptFile <- list.files(pattern = paste0("^", scriptNum, "_.*\\.R$"))[1]
  
  if (!is.na(scriptFile) && file.exists(scriptFile)) {
    lines <- readLines(scriptFile)
    lines <- gsub("^setwd\\(", "# setwd(", lines)
    writeLines(lines, scriptFile)
  }
}

#### Run Preprocessing Scripts ----
source("1_PreprocessMetabo.R")
source("2_ImputeMetabo.R")
source("3_FilterMetabo.R")
source("4_FilterRNA.R")
source("5_Concordance.R")

cat("✓ Preprocessing complete\n")

setwd(here())
