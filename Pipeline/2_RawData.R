#### Initialize ----
rm(list = ls())

#### Load Required Packages ----
if (!require("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

#### Set Download Configuration ----
# increase timeout for large file downloads (~700MB total)
options(timeout = 10000)

rawDataDirectory       <- here("data", "raw")
processedDataDirectory <- here("data", "processed")

# create directories if they don't exist
if (!dir.exists(rawDataDirectory))       dir.create(rawDataDirectory,       recursive = TRUE)
if (!dir.exists(processedDataDirectory)) dir.create(processedDataDirectory, recursive = TRUE)

#### Download Molecular Data ----
# pancancer metabolomics dataset from Zenodo (record 7348648)
# ~500MB compressed | ~800 metabolites x ~900 samples
cat("\nDownloading pancancer metabolomics data...\n")

tempMolecularArchive <- tempfile()

download.file(
  url      = "https://zenodo.org/record/7348648/files/pancancer_metabolomics_v.0.3.2.tar.gz?download=1",
  destfile = tempMolecularArchive,
  mode     = "wb"
)

untar(tarfile = tempMolecularArchive, exdir = rawDataDirectory)
unlink(tempMolecularArchive)

cat("✓ Molecular data downloaded\n")

#### Download Precalculated Data ----
# precomputed analysis results from Zenodo (record 7352546)
# ~200MB compressed
cat("\nDownloading precalculated data...\n")

tempPrecalcArchive <- tempfile()

download.file(
  url      = "https://zenodo.org/record/7352546/files/data_for_scripts_v0.3.3.tar.gz?download=1",
  destfile = tempPrecalcArchive,
  mode     = "wb"
)

untar(tarfile = tempPrecalcArchive, exdir = rawDataDirectory)
unlink(tempPrecalcArchive)

cat("✓ Precalculated data downloaded\n")

#### Copy Data to Working Location ----
# keep raw data pristine, work on copies
cat("\nOrganizing data...\n")

downloadedDataSource   <- here("data", "raw",       "pancancer_metabolomics")
workingDataDestination <- here("data", "processed", "pancancer_metabolomics")

if (!dir.exists(downloadedDataSource)) {
  stop("ERROR: Downloaded data not found at expected location: ", downloadedDataSource)
}

if (!dir.exists(workingDataDestination)) {
  file.copy(
    from      = downloadedDataSource,
    to        = processedDataDirectory,
    recursive = TRUE,
    overwrite = FALSE
  )
  cat("✓ Data copied to working directory\n")
} else {
  cat("✓ Working copy already exists\n")
}

#### Verify ----
cat("\n=== Data structure ===\n")
cat("Raw data folders:\n")
print(list.dirs(rawDataDirectory,       recursive = FALSE, full.names = FALSE))

cat("\nProcessed data folders:\n")
print(list.dirs(workingDataDestination, recursive = FALSE, full.names = FALSE))

cat("\n✓ Download and setup complete | Total: ~700MB\n")
