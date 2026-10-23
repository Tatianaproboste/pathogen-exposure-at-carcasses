# -------------------------------------------------------------------------
# Interspecies Interaction Network - Data Preparation
# -------------------------------------------------------------------------
# Purpose: Create species-by-time-and-station groups for network analysis
# Author: Tatiana Proboste
# Date: 2026-05-15
# -------------------------------------------------------------------------

# Load libraries ----------------------------------------------------------
suppressMessages({
  library(data.table)
  library(spatsoc)
  library(dplyr)
  library(readr)
  library(lubridate)
})

rm(list = ls())

# Setup paths -------------------------------------------------------------
INPUT_FILE <- "./Raw_data/species_clean.csv"

OUTPUT_DIR <- "./Outputs/Networks"
OUTPUT_FILE <- file.path(OUTPUT_DIR, "outputs_full.rds")

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load data ---------------------------------------------------------------
dat <- read_csv(INPUT_FILE, show_col_types = FALSE)

#check years available
table(dat$year, dat$location)

# Create temporal groups with MULTIPLE thresholds -------------------------

setDT(dat)

# Define time windows
time_windows <- c("24 hours", "2 days", "3 days")

# Create a list to store results for each window
outputs <- list()

for(window in time_windows) {
  cat(paste0("\n  Processing ", window, " window...\n"))
  
  # Create a copy to avoid overwriting
  tmp_data <- copy(dat)
  
  # Group observations within the specified window
  tmp_data <- spatsoc::group_times(
    tmp_data, 
    datetime = "datetime", 
    threshold = window
  )
  
  cat(paste0("    Time groups created: ", uniqueN(tmp_data$timegroup), "\n"))
  
  # Create spatial-temporal interaction groups
  tmp_data[, group := .GRP, by = .(station, timegroup)]
  
  cat(paste0("    Interaction groups: ", uniqueN(tmp_data$group), "\n"))
  
  # Add ID column for network analysis
  tmp_data[, id := paste(species, location, station, year, season, sep = "_")]
  
  # Calculate summary statistics
  group_summary <- tmp_data[, .(
    n_observations = .N,
    n_species = uniqueN(species),
    mean_time_span = as.numeric(difftime(max(datetime), min(datetime), units = "hours"))
  ), by = group]
  # Store in list
  outputs[[window]] <- tmp_data
}

# Save all outputs --------------------------------------------------------
# Save individual files
saveRDS(outputs[["24 hours"]], file.path(OUTPUT_DIR, "outputs_24h.rds"))
saveRDS(outputs[["2 days"]], file.path(OUTPUT_DIR, "outputs_48h.rds"))
saveRDS(outputs[["3 days"]], file.path(OUTPUT_DIR, "outputs_72h.rds"))

# Also save as a combined list
saveRDS(outputs, file.path(OUTPUT_DIR, "outputs_all_windows.rds"))

# For backward compatibility, save 24h as the default
saveRDS(outputs[["24 hours"]], OUTPUT_FILE)

