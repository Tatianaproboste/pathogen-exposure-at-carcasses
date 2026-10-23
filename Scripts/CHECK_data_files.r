# CHECK: Are 24h and 72h files actually different?
# =========================================================================

library(data.table)
library(tidyverse)

INPUT_FOLDER <- "./Outputs/Networks/"

# Load both files
cat("Loading 24h data...\n")
df_24h <- readRDS(file.path(INPUT_FOLDER, "outputs_24h.rds")) %>%
  janitor::clean_names() %>%
  as.data.table()

cat("Loading 72h data...\n")
df_72h <- readRDS(file.path(INPUT_FOLDER, "outputs_72h.rds")) %>%
  janitor::clean_names() %>%
  as.data.table()

# Compare
cat("\n=== COMPARISON ===\n")
cat("24h rows:", nrow(df_24h), "\n")
cat("72h rows:", nrow(df_72h), "\n")
cat("24h groups:", uniqueN(df_24h$group), "\n")
cat("72h groups:", uniqueN(df_72h$group), "\n")
cat("24h species:", uniqueN(df_24h$species), "\n")
cat("72h species:", uniqueN(df_72h$species), "\n\n")

# Check if groups are identical
cat("Are group IDs identical?", identical(unique(df_24h$group), unique(df_72h$group)), "\n")
cat("Are entire datasets identical?", identical(df_24h, df_72h), "\n\n")

# Sample from each
cat("=== SAMPLE FROM 24h ===\n")
print(head(df_24h[, .(group, species, datetime, location)], 15))

cat("\n=== SAMPLE FROM 72h ===\n")
print(head(df_72h[, .(group, species, datetime, location)], 15))

# Check one location in detail
cat("\n=== KNP COMPARISON ===\n")
knp_24 <- df_24h[location == "KNP"]
knp_72 <- df_72h[location == "KNP"]

cat("KNP 24h: ", nrow(knp_24), "rows,", uniqueN(knp_24$group), "groups\n")
cat("KNP 72h: ", nrow(knp_72), "rows,", uniqueN(knp_72$group), "groups\n")

cat("\n=== DIAGNOSIS ===\n")
if(identical(df_24h, df_72h)) {
  cat("❌ PROBLEM: Both files are IDENTICAL!\n")
  cat("   The 24h and 72h files contain the exact same data.\n")
  cat("   You need to regenerate these files with proper time windows.\n")
} else if(uniqueN(df_24h$group) == uniqueN(df_72h$group)) {
  cat("⚠️  WARNING: Same number of groups in both files!\n")
  cat("   This is suspicious - 72h should have FEWER groups (more aggregated)\n")
  cat("   Check how the 'group' variable was created.\n")
} else {
  cat("✓ Files appear different.\n")
  cat("  The problem may be elsewhere (network construction or c-z calculation)\n")
}
