# =========================================================================
# DISEASE EXPOSURE RISK MODELS - MULTIPLE DISEASES
# =========================================================================
#
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(data.table)
library(viridis)
library(patchwork)
library(tidytext)
library(scales)
library(ggplot2)

rm(list = ls())

# SOURCE SHARED FUNCTIONS -------------------------------------------------
# Keep shared objects in this script environment when sourced locally.
source("Scripts/_functions_shared.R", local = TRUE)

# Consistent risk gradients used across visualizations in this script
risk_palette <- unname(risk_colors[c("low", "medium", "high", "veryhigh")])

# Location labels: short code → full region name
location_labels <- c(KNP = "Alpine", SD = "Desert", WV = "Temperate")

# SETUP -------------------------------------------------------------------
INPUT_FILE <- "./Raw_data/species_clean.csv"
DISEASE_PARAMS_FILE <- "./Raw_data/disease_parameters.csv"
SPECIES_SUSCEPTIBILITY_FILE <- "./Raw_data/species_susceptibility.csv"
OUTPUT_DIR <- "./Outputs/Disease_Risk_Models"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# LOAD DATA ---------------------------------------------------------------

# Load disease parameters
disease_params <- read_csv(DISEASE_PARAMS_FILE, show_col_types = FALSE) %>%
  mutate(
    # Convert resistance period (λ) to hours for exponential decay calculation
    resistance_hours = resistance_days * 24
  )

# Load species data
df <- read_csv(INPUT_FILE, show_col_types = FALSE) %>%
  janitor::clean_names()

cat(paste0("  Data loaded: ", nrow(df), " observations\n"))
cat(paste0("  Species: ", length(unique(df$species)), "\n"))
cat(paste0("  Carcasses: ", length(unique(df$station)), "\n\n"))

# =========================================================================
# STEP 1: CALCULATE ARRIVAL TIMING & PATHOGEN EXPOSURE
# =========================================================================

# Get carcass deployment times (keep location so it flows through to risk summaries)
# If carcass_deploy_time is missing, fall back to the earliest detection at that station
carcass_deployment <- df %>%
  group_by(station) %>%
  summarise(
    deployment_time_raw = first(carcass_deploy_time),
    earliest_detection  = min(datetime, na.rm = TRUE),
    location            = first(location),
    .groups = "drop"
  ) %>%
  mutate(
    deployment_time = if_else(
      !is.na(deployment_time_raw),
      as.POSIXct(deployment_time_raw),
      as.POSIXct(earliest_detection)
    )
  ) %>%
  select(station, deployment_time, location)

n_with_deploy <- sum(!is.na(df$carcass_deploy_time))

# Calculate when each species arrived at each carcass
arrival_timing <- df %>%
  group_by(station, species) %>%
  summarise(
    first_arrival = min(datetime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(carcass_deployment, by = "station") %>%
  mutate(
    # Hours since carcass deployed (or first detection if deployment time missing)
    hours_since_deployment = as.numeric(difftime(first_arrival, deployment_time, units = "hours")),
    hours_since_deployment = pmax(0, hours_since_deployment)
  )

# Add pathogen availability for each disease dynamically
for (i in 1:nrow(disease_params)) {
  disease_code <- disease_params$disease[i]
  resistance_days <- disease_params$resistance_days[i]
  resistance_hours <- resistance_days * 24  # Convert resistance days to hours
  
  # Calculate pathogen availability with exponential decay
  # Pathogen decays exponentially: P(t) = exp(-t / λ) where λ = resistance period
  # Animals arriving early (small t) have high exposure (P ≈ 1)
  # Animals arriving late (large t) have low exposure (P → 0)
  arrival_timing <- arrival_timing %>%
    mutate(
      !!paste0("pathogen_", disease_code, "_raw") := exp(-hours_since_deployment / !!resistance_hours),
      !!paste0("pathogen_", disease_code) := ifelse(
        .data[[paste0("pathogen_", disease_code, "_raw")]] < 0.01, 
        0, 
        .data[[paste0("pathogen_", disease_code, "_raw")]]
      )
    )
}

arrival_timing <- arrival_timing %>%
  select(-ends_with("_raw"))  # Remove intermediate columns

cat(paste0("  Species-carcass combinations: ", nrow(arrival_timing), "\n"))
cat(paste0("  Arrival time range: ", round(min(arrival_timing$hours_since_deployment), 1), 
           " to ", round(max(arrival_timing$hours_since_deployment), 1), " hours\n\n"))

# Report pathogen resistance periods (λ = persistence time)
cat("  Pathogen resistance periods (λ):\n")
for (i in 1:nrow(disease_params)) {
  disease_name <- disease_params$disease[i]
  resist_hours <- disease_params$resistance_hours[i]
  resist_days <- disease_params$resistance_days[i]
  
  if (resist_days > 365) {
    cat(paste0("    ", disease_name, ": λ = ", round(resist_hours, 0), " hours (", 
               round(resist_days / 365, 1), " years)\n"))
  } else {
    cat(paste0("    ", disease_name, ": λ = ", round(resist_hours, 0), " hours (", 
               round(resist_days, 1), " days)\n"))
  }
}
cat("\n")

# =========================================================================
# STEP 2: DEFINE DISEASE SUSCEPTIBILITY
# =========================================================================

# Load species susceptibility from file
susceptibility_matrix <- read_csv(SPECIES_SUSCEPTIBILITY_FILE, show_col_types = FALSE) %>%
  rename(species = Species)  # Rename to match data column name

# Get all disease columns dynamically (exclude non-disease columns)
disease_cols <- setdiff(names(susceptibility_matrix), 
                       c("species", "functional_group", "notes"))

# Keep only needed columns
susceptibility_matrix <- susceptibility_matrix %>%
  select(species, functional_group, all_of(disease_cols))

# Validate that all diseases in parameters file have susceptibility columns
missing_diseases <- setdiff(disease_params$disease, disease_cols)

# Check for missing species
all_species <- unique(df$species)
missing_species <- setdiff(all_species, susceptibility_matrix$species)

if (length(missing_species) > 0) {
  cat(paste0("  ⚠ WARNING: ", length(missing_species), " species in data but not in susceptibility file:\n"))
  cat(paste0("    ", paste(missing_species, collapse = ", "), "\n"))
  cat("    Adding conservative default susceptibility values (low risk)...\n")
  cat("    ➜ RECOMMENDATION: Add these species to species_susceptibility.csv with proper values\n\n")}
  
  # Create default values dynamically - use conservative (low) defaults for all diseases
  default_values <- list(
    species = missing_species,
    functional_group = "Unknown"
  )
  
  # Add low default susceptibility for each disease (0.05 for all)
  for (disease in disease_cols) {
    default_values[[disease]] <- 0.05  # Conservative default
  }
  
  default_susceptibility <- as_tibble(default_values)
  
# =========================================================================
# STEP 3: CALCULATE RISK FOR EACH DISEASE
# =========================================================================

# Join susceptibility with arrival timing
risk_per_carcass <- arrival_timing %>%
  left_join(susceptibility_matrix, by = "species")

# Calculate risk for each disease using formula: Risk = S × P × T
# Where S = susceptibility, P = pathogen availability (exp decay), T = transmission rate

for (i in 1:nrow(disease_params)) {
  disease_code <- disease_params$disease[i]
  trans_rate <- disease_params$transmission_rate[i]
  
  risk_per_carcass <- risk_per_carcass %>%
    mutate(
      # Risk = Susceptibility (S) × Pathogen_availability (P) × Transmission_rate (T)
      !!paste0(disease_code, "_risk") := .data[[disease_code]] *                        # S
                                          .data[[paste0("pathogen_", disease_code)]] *  # P
                                          !!trans_rate                                  # T
    )
}

# Aggregate by species (mean risk across all carcass visits)
# Build the summarise call dynamically
risk_summary_expr <- rlang::exprs(
  n_carcasses = n(),
  mean_hours_since_deployment = mean(hours_since_deployment),
  median_arrival_hours = median(hours_since_deployment),
  min_hours = min(hours_since_deployment),
  max_hours = max(hours_since_deployment),
  locations = paste(sort(unique(location)), collapse = ", ")
)

# Add risk and susceptibility columns for each disease
for (i in 1:nrow(disease_params)) {
  disease_code <- disease_params$disease[i]
  
  # Add risk column
  risk_summary_expr[[paste0(disease_code, "_risk")]] <- 
    rlang::expr(mean(!!rlang::sym(paste0(disease_code, "_risk")), na.rm = TRUE))
  
  # Add susceptibility column
  risk_summary_expr[[paste0(disease_code, "_susceptibility")]] <- 
    rlang::expr(first(!!rlang::sym(disease_code)))
}

# Create the summary
risk_summary <- risk_per_carcass %>%
  group_by(species) %>%
  summarise(!!!risk_summary_expr, .groups = "drop")

# Sort by total risk across all diseases
risk_cols <- paste0(disease_params$disease, "_risk")
risk_summary <- risk_summary %>%
  mutate(total_risk = rowSums(across(all_of(risk_cols)))) %>%
  arrange(desc(total_risk)) %>%
  select(-total_risk)

cat(paste0("  Risk calculated for ", nrow(risk_summary), " species\n\n"))

# =========================================================================
# STEP 4: SAVE RESULTS
# =========================================================================


# Save detailed results
write_csv(risk_summary, file.path(OUTPUT_DIR, "multi_disease_risk_summary_SIMPLE.csv"))
write_csv(risk_per_carcass, file.path(OUTPUT_DIR, "multi_disease_risk_per_visit_SIMPLE.csv"))

# =========================================================================
# STEP 5: CREATE VISUALIZATIONS
# =========================================================================

# Comparative risk heatmap ---
risk_heatmap_data <- risk_summary %>%
  select(species, ASF_risk, FMD_risk, Anthrax_risk, Leptospirosis_risk, LSD_risk, Hydatidosis_risk) %>%
  pivot_longer(cols = ends_with("_risk"), names_to = "disease", values_to = "risk") %>%
  mutate(
    disease = str_remove(disease, "_risk"),
    risk_category = case_when(
      risk == 0 ~ "None",
      risk <= 0.3 ~ "Low",
      risk <= 0.6 ~ "Medium",
      TRUE ~ "High"
    ),
    risk_category = factor(risk_category, levels = c("None", "Low", "Medium", "High"))
  )

# Disease-specific bubble plots
# Function to create bubble plot for a specific disease
create_disease_bubble_plot <- function(risk_data, disease_name, disease_fullname,
                                       threshold_hours, color_high, size_range = c(10, 25),
                                       label_size = 3.5, species_below = NULL,
                                       fill_limits = c(0, NA)) {
  
  # Dynamic column names
  risk_col <- paste0(disease_name, "_risk")
  suscept_col <- paste0(disease_name, "_susceptibility")
  
  # Get top 5 species for this disease
  top5 <- risk_data %>%
    arrange(desc(.data[[risk_col]])) %>%
    head(5) %>%
    mutate(
      risk_label = paste0(round(.data[[risk_col]] * 100, 1), "%"),
      species_wrapped = str_wrap(species, width = 15),
      label_vjust = ifelse(species %in% species_below, 2.5, -1.8)
    )
  
  # Calculate appropriate x-axis limit based on actual arrival times
  max_arrival <- max(top5$max_hours, na.rm = TRUE)
  
  # Special handling for LSD - needs longer x-axis to show risk window
  if (disease_name == "LSD") {
    x_max <- max(max_arrival * 1.3, threshold_hours * 1.15)  # Extend to show risk window
    x_max <- min(x_max, 600)  # Cap at 600 hours for LSD
    x_breaks_interval <- 50
  } else if (threshold_hours > 1000) {
    # For very long thresholds (Anthrax, Brucella, Hydatidosis), just use max arrival + buffer
    x_max <- min(max_arrival * 1.2, 1000)
    x_breaks_interval <- 100
  } else {
    # For shorter thresholds, use threshold as reference
    x_max <- min(max_arrival * 1.15, threshold_hours * 1.2, 1000)
    x_breaks_interval <- 50
  }
  
  x_breaks_max <- floor(x_max / x_breaks_interval) * x_breaks_interval
  
  # Create plot
  p <- ggplot(top5, aes(x = mean_hours_since_deployment, y = .data[[suscept_col]])) +
    
    # Background risk zones - subtle, muted colors
    annotate("rect", xmin = 0, xmax = threshold_hours, ymin = -Inf, ymax = Inf, 
             fill = "#FEF3E2", alpha = 0.3) +  # Very light warm beige for risk window
    annotate("rect", xmin = threshold_hours, xmax = Inf, ymin = -Inf, ymax = Inf, 
             fill = "#E8E8E8", alpha = 0.2) +  # Very light grey for outside risk window
    
    # Error bars showing min/max arrival times - more subtle
    geom_errorbarh(aes(xmin = min_hours, xmax = max_hours), 
                   height = 0.02, color = "grey50", linewidth = 0.4, alpha = 0.5) +
    
    # Bubbles - cleaner stroke
    geom_point(aes(size = .data[[risk_col]], fill = .data[[risk_col]]), 
               shape = 21, alpha = 0.9, color = "grey20", stroke = 0.5) +
    
    # Species labels with risk percentage - cleaner font
    geom_text(aes(label = paste0(species_wrapped, "\n(", risk_label, ")"), vjust = label_vjust), 
              size = label_size, fontface = "plain", lineheight = 0.85, family = "sans") +
    
    # Risk threshold line - more subtle
    geom_vline(xintercept = threshold_hours, linetype = "dashed", 
               color = "#E76F51", linewidth = 0.6, alpha = 0.8) +
    annotate("text", x = threshold_hours, y = 1.05,
             label = paste0("Risk window\n(", threshold_hours, " h)"),
             hjust = 0.5, size = 3, color = "#E76F51", fontface = "plain", family = "sans") +
    
    # Color scale - journal-style gradient (purple-navy → teal → yellow)
    scale_fill_gradientn(colors = color_palette(100),
                         name = "Mean risk",
                         labels = scales::label_percent(accuracy = 1),
                         limits = fill_limits) +
    scale_size_continuous(range = size_range, guide = "none") +
    
    # Axes
    scale_x_continuous(limits = c(0, x_max), 
                       breaks = seq(0, x_breaks_max, x_breaks_interval),
                       expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(-0.05, 1.1), breaks = seq(0, 1, 0.2),
                       labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
                       expand = c(0.01, 0)) +
    
    # Labels - clean, sentence-case styling
    labs(
      title = disease_fullname,
      x = "Mean arrival time (hours since carcass deployment)",
      y = "Species biological susceptibility"
    ) +
    
    # Theme - journal style
    theme_minimal(base_size = 11, base_family = "sans") +
    theme(
      # Clean white background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),

      # Minimal grid
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),

      # Thin or no panel border
      panel.border = element_blank(),
      axis.line = element_line(color = "grey30", linewidth = 0.4),

      # Typography - sans-serif, restrained
      plot.title = element_text(face = "bold", size = 16, margin = margin(b = 8), hjust = 0),
      axis.title = element_text(face = "plain", size = 16, color = "grey20"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      axis.text = element_text(size = 9, color = "grey30"),
      axis.ticks = element_line(color = "grey30", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),

      # Legend - compact and unobtrusive
      legend.position = "right",
      legend.background = element_blank(),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.4, "cm"),
      legend.title = element_text(face = "plain", size = 9, color = "grey20"),
      legend.text = element_text(size = 8, color = "grey30"),
      legend.margin = margin(l = 5),

      # Tight margins
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    )
  
  return(p)
}

# Generate and save bubble plot for each disease
for (i in 1:nrow(disease_params)) {
  disease <- disease_params$disease[i]
  fullname <- disease_params$disease_fullname[i]
  threshold <- disease_params$resistance_hours[i]  # Use resistance period (λ in hours)
  color <- disease_params$color_gradient[i]
  
 
  # Define species to place below bubbles (for better label spacing)
  species_below <- if(disease == "Anthrax") {
    c("Dingo", "Red Fox")  # Move these two species' labels below their bubbles
  } else if (disease == "Hydatidosis") {
    # For Hydatidosis some labels overlap at high susceptibility; move these below bubbles
    c("Eastern Grey Kangaroo", "Dingo", "Red Fox")
  } else {
    NULL  # Default: all labels above
  }
  
  # Create plot
  p <- create_disease_bubble_plot(
    risk_data = risk_summary,
    disease_name = disease,
    disease_fullname = fullname,
    threshold_hours = threshold,
    color_high = color,
    species_below = species_below
  )
  
  # Save plot
  filename <- paste0(disease, "_risk_top5_bubble.png")
  ggsave(file.path(OUTPUT_DIR, filename), p,
         width = 10, height = 6, dpi = 600, bg = "white")  # Standard width with wrapped names
}

# add all the bubble plots together in a single figure in 2 columns
library(patchwork)

# Shared fill scale so every panel's "Mean risk" legend is identical and can
# be collected into a single legend instead of repeating per panel
global_risk_max <- risk_summary %>%
  select(all_of(risk_cols)) %>%
  unlist() %>%
  max(na.rm = TRUE)
shared_fill_limits <- c(0, global_risk_max)

# Create plots (Figure 4) with panel labels (A, B, C, etc.)
plot_list <- lapply(seq_along(disease_params$disease), function(i) {
  disease <- disease_params$disease[i]

  # Define species to place below bubbles (for better label spacing)
  species_below <- if(disease == "Anthrax") {
    c("Dingo", "Red Fox")  # Move these two species' labels below their bubbles
  } else if (disease == "Hydatidosis") {
    # Mirror single-plot logic for combined figure
    c("Eastern Grey Kangaroo", "Cat", "Red Fox")
  } else {
    NULL  # Default: all labels above
  }

  disease_fullname <- disease_params$disease_fullname[disease_params$disease == disease]

  p <- create_disease_bubble_plot(
    risk_data = risk_summary,
    disease_name = disease,
    disease_fullname = disease_fullname,
    threshold_hours = disease_params$resistance_hours[disease_params$disease == disease],  # Use λ
    color_high = disease_params$color_gradient[disease_params$disease == disease],
    size_range = c(4, 12),  # Smaller bubbles for combined plot
    label_size = 2.2,  # Smaller text labels for combined plot
    species_below = species_below,
    fill_limits = shared_fill_limits  # Same scale on every panel → one shared legend
  )

  # Folded into the title (instead of a floating plot.tag) so it never
  # overlaps the in-panel "Risk window" annotation, which for short-λ
  # diseases can sit near the same top-left corner as a floating tag.
  panel_label <- LETTERS[i]
  p <- p +
    labs(title = paste0(panel_label, ") ", disease_fullname))

  return(p)
})

combined_plot <- wrap_plots(plotlist = plot_list, ncol = 2) +
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(
    theme = theme(plot.margin = margin(5, 5, 5, 5))
  ) &
  theme(legend.position = "right")
  
ggsave(file.path(OUTPUT_DIR, "combined_bubble_plot.png"), combined_plot,
       width = 16, height = 20, dpi = 600, bg = "white")
