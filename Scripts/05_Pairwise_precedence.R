# =========================================================================
# PAIRWISE PRECEDENCE — with bootstrap CIs and minimum-n filter
# For every pair of species seen at the same carcass, score which arrived
# first (1 = A first, 0 = B first, 0.5 = tie). A species' precedence score
# is its mean score across all pairs; CIs come from resampling carcasses.
# =========================================================================

library(data.table)
library(ggplot2)
library(janitor)

source("Scripts/_functions_shared.R")

# SETTINGS ----------------------------------------------------------------
INPUT_FOLDER         <- "./Outputs/Networks/"
FIRST_ARRIVAL_DIR    <- "./Outputs/first_arrival/"
TIME_WINDOW          <- "24h"
TIE_SECONDS          <- 60     # arrivals this close count as a tie
COOCCUR_WINDOW       <- Inf    # max hours between arrivals to count as a pair
MIN_COOCCUR_PAIR     <- 1      # minimum co-occurrences to retain a SPECIES PAIR
MIN_COOCCUR_SPECIES  <- 1      # minimum total co-occurrences to plot a SPECIES
N_BOOT               <- 1000   # bootstrap replicates


dir.create(FIRST_ARRIVAL_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# DATA --------------------------------------------------------------------
# datetime is already POSIXct in the RDS; each station is one carcass
df <- load_network_rds(INPUT_FOLDER, TIME_WINDOW)
df <- df[!is.na(datetime) & !is.na(station) & !is.na(species)]
df[, carcass := station]

# Functions -----------------------------------------------------------------
# All ordered species pairs (A, B) at each carcass, with each species' first arrival
first_arrival_pairs <- function(dt) {
  arrivals <- dt[, .(first_time = min(datetime)), by = .(carcass, species)]
  pairs <- arrivals[arrivals, on = "carcass", allow.cartesian = TRUE]
  setnames(pairs, c("species", "first_time", "i.species", "i.first_time"), c("A", "t_A", "B", "t_B"))
  pairs[A != B]
}

# Score each pair: 1 if A arrived first, 0 if B did, 0.5 if within TIE_SECONDS
score_pairs <- function(dt) {
  pairs <- first_arrival_pairs(dt)
  pairs[, gap_sec := as.numeric(difftime(t_B, t_A, units = "secs"))]
  pairs <- pairs[abs(gap_sec) <= COOCCUR_WINDOW * 3600]
  pairs[, contribution := fifelse(abs(gap_sec) <= TIE_SECONDS, 0.5, fifelse(gap_sec > 0, 1, 0))]
  pairs[, .(carcass, A, B, contribution)]
}

# Mean score per species, with 95% CI from resampling carcasses with replacement
bootstrap_precedence_scores <- function(events) {
  carcasses <- unique(events$carcass)
  point_est <- events[, .(score = mean(contribution), total_n = .N), by = .(species = A)]

  boot_mat <- matrix(NA_real_, nrow = N_BOOT, ncol = nrow(point_est),
                     dimnames = list(NULL, point_est$species))
  for (b in seq_len(N_BOOT)) {
    samp <- events[.(sample(carcasses, replace = TRUE)), on = "carcass", allow.cartesian = TRUE]
    b_scores <- samp[, .(score = mean(contribution)), by = .(species = A)]
    boot_mat[b, b_scores$species] <- b_scores$score
  }

  # species missing from too many replicates get no CI
  min_boot_non_na <- max(50, round(N_BOOT * 0.05))
  cis <- apply(boot_mat, 2, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < min_boot_non_na) return(c(NA_real_, NA_real_))
    unname(quantile(x, c(0.025, 0.975)))
  })

  point_est[, `:=`(ci_lower = cis[1, ], ci_upper = cis[2, ])]
}

# DIAGNOSTICS
pair_counts <- first_arrival_pairs(df)[, .N, by = .(A, B)]$N
cat(sprintf("PAIRWISE: pairs=%d, mean=%.2f, median=%.1f, max=%d\n",
            length(pair_counts), mean(pair_counts), median(pair_counts), max(pair_counts)))

png(file.path(FIRST_ARRIVAL_DIR, "pre_diag_pair_counts.png"), width = 800, height = 400)
hist(pair_counts, breaks = seq(0, max(pair_counts) + 1, by = 1), col = "#2A1B5D",
     xlab = "n_cooccur", main = "Pairwise co-occurrence counts")
abline(v = MIN_COOCCUR_PAIR, col = "red", lty = 2)
dev.off()

ths <- 1:min(10, max(pair_counts))
cum_pairs <- sapply(ths, function(t) sum(pair_counts >= t))
cat("Pairs retained at thresholds 1..10: ", paste(ths, cum_pairs, sep = "=", collapse = ", "), "\n")

# RUN PER LOCATION --------------------------------------------------------
precedence_scores <- list()
precedence_pairs  <- list()

for (loc in unique(df$location)) {
  message("Processing: ", loc)
  events <- score_pairs(df[location == loc])
  if (nrow(events) == 0) next

  # per-pair summary
  precedence_pairs[[loc]] <- events[, .(
    n_cooccur    = .N,
    n_A_first    = sum(contribution),
    p_A_before_B = mean(contribution)
  ), by = .(A, B)][n_cooccur >= MIN_COOCCUR_PAIR][, location := loc]

  precedence_scores[[loc]] <- bootstrap_precedence_scores(events)[, location := loc]
}

precedence_score <- rbindlist(precedence_scores)
precedence_all   <- rbindlist(precedence_pairs)

precedence_score[, guild := assign_guild(species)]
precedence_score[, distinguishable := ci_lower > 0.5 | ci_upper < 0.5]  # CI excludes 0.5

# FIGURE ------------------------------------------------------------------
plot_dt <- precedence_score[total_n >= MIN_COOCCUR_SPECIES]
message("Species retained: ", nrow(plot_dt), " of ", nrow(precedence_score))

# order species within each location so the earliest arrivers sit at the top
plot_dt[, species_re := reorder_within(species, score, location)]

location_labels <- c(KNP = "Alpine", SD = "Desert", WV = "Temperate")

p_marg <- ggplot(plot_dt, aes(x = species_re, y = score)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#6E6E6E") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, linewidth = 0.6, alpha = 0.85) +
  geom_point(aes(size = total_n), color = "#2A1B5D", alpha = 0.95) +
  facet_wrap(~ location, scales = "free_y", ncol = 3,
             labeller = labeller(location = location_labels)) +
  coord_flip() +
  scale_size_continuous(name = "Co-occurrences", range = c(3, 8)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_reordered() +
  labs(x = NULL,
       y = "Relative arrival order among scavenger species (95% bootstrap CI)") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(
    axis.title   = element_text(size = 14, color = "grey20"),
    axis.text    = element_text(size = 13, color = "grey30"),
    axis.text.y  = element_text(face = "italic", size = 13),
    strip.text   = element_text(face = "bold", size = 14),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12)
  )

ggsave(file.path(FIRST_ARRIVAL_DIR, "precedence_score_by_species_CI.png"), p_marg,
       width = 14, height = 7, dpi = 600, bg = "white")

# SAVE --------------------------------------------------------------------
fwrite(precedence_score, file.path(FIRST_ARRIVAL_DIR, "precedence_score_with_CI.csv"))
fwrite(precedence_all,   file.path(FIRST_ARRIVAL_DIR, "pairwise_precedence.csv"))

message("Done — outputs written to: ", FIRST_ARRIVAL_DIR)
