# Disease Transmission Risk: Species Centrality Analysis
# Identifies superspreader species across network contexts
# =========================================================================

library(igraph)
library(tidyverse)
library(viridis)
library(patchwork)
library(RColorBrewer)

# SOURCE SHARED FUNCTIONS -------------------------------------------------
source("Scripts/_functions_shared.R", local = TRUE)
gradient_palette <- color_palette

# Location labels: short code → full region name (consistent across all scripts)
location_labels <- c(KNP = "Alpine", SD = "Desert", WV = "Temperate")
OUTPUT_DIR <- "./Outputs/Networks/centrality measures"

# =========================================================================
# 1. LOAD DATA & DEFINE FUNCTIONAL GROUPS
# =========================================================================

# Load networks (structure: list(cooccurrence = list(...), temporal = list(...)))
networks_file <- readRDS("./Outputs/Networks/structure_analysis/networks_v2.rds")

# Extract co-occurrence networks (undirected - best for centrality analysis)
all_networks <- networks_file$cooccurrence

# =========================================================================
# 2. CALCULATE CENTRALITY METRICS
# =========================================================================

calculate_centrality <- function(network, network_name) {
  
  species <- V(network)$name
  
  # Raw metrics
  degree_cent <- degree(network)
  betweenness_cent <- betweenness(network, normalized = FALSE)
  closeness_cent <- closeness(network, normalized = FALSE)
  eigen_cent <- eigen_centrality(network)$vector
  
  # Normalize to 0-1 scale
  normalize_01 <- function(x) {
    if(max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) return(rep(0.5, length(x)))
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  
  degree_norm <- normalize_01(degree_cent)
  betweenness_norm <- normalize_01(betweenness_cent)
  closeness_norm <- normalize_01(closeness_cent)
  eigen_norm <- normalize_01(eigen_cent)
  
  data.frame(
    network = network_name,
    species = species,
    functional_group = sapply(species, assign_guild),
    degree_norm = degree_norm,
    betweenness_norm = betweenness_norm,
    closeness_norm = closeness_norm,
    eigenvector_norm = eigen_norm,
    stringsAsFactors = FALSE
  )
}

# Process all networks
centrality_list <- lapply(names(all_networks), function(network_name) {
  network <- all_networks[[network_name]]
  if (vcount(network) >= 3) {
    calculate_centrality(network, network_name)
  } else {
    NULL
  }
})

centrality_df <- bind_rows(centrality_list) %>%
  mutate(
    net_parts = strsplit(network, "_"),
    time_window = vapply(net_parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1)),
    year = vapply(
      net_parts,
      function(x) {
        if (length(x) >= 3 && grepl("^[0-9]{4}$", x[2])) x[2] else NA_character_
      },
      character(1)
    ),
    location = vapply(
      net_parts,
      function(x) {
        if (length(x) == 2) {
          x[2]
        } else if (length(x) >= 3 && grepl("^[0-9]{4}$", x[2])) {
          x[3]
        } else if (length(x) >= 3) {
          x[2]
        } else {
          NA_character_
        }
      },
      character(1)
    )
  ) %>%
  select(-net_parts) %>%
  mutate(location_year = ifelse(is.na(year), location, paste(location, year, sep = "_")))

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Save full results
write.csv(centrality_df, 
          file.path(OUTPUT_DIR, "species_centrality_all.csv"),
          row.names = FALSE)

# =========================================================================
# 3. IDENTIFY CONSISTENTLY IMPORTANT SPECIES (BY DEGREE)
# =========================================================================

# Species in top 10 for ≥50% of contexts (using degree as primary metric)
consistency <- centrality_df %>%
  group_by(network) %>%
  mutate(is_top10 = rank(-degree_norm, ties.method = "min") <= 10) %>%
  group_by(species, functional_group) %>%
  summarise(
    times_in_top10 = sum(is_top10),
    pct_contexts = 100 * mean(is_top10),
    mean_degree = mean(degree_norm),
    mean_betweenness = mean(betweenness_norm),
    mean_closeness = mean(closeness_norm),
    mean_eigenvector = mean(eigenvector_norm),
    .groups = "drop"
  ) %>%
  filter(pct_contexts >= 50) %>%
  arrange(desc(mean_degree))

important_species <- consistency$species

print(consistency)

# Save top species by functional group
summary_by_group <- centrality_df %>%
  group_by(functional_group, species) %>%
  summarise(
    mean_degree = mean(degree_norm),
    mean_betweenness = mean(betweenness_norm),
    mean_closeness = mean(closeness_norm),
    mean_eigenvector = mean(eigenvector_norm),
    n_contexts = n(),
    .groups = "drop"
  ) %>%
  group_by(functional_group) %>%
  slice_max(mean_degree, n = 3) %>%
  ungroup()

write.csv(summary_by_group,
          file.path(OUTPUT_DIR, "top_species_by_functional_group.csv"),
          row.names = FALSE)

# =========================================================================
# 4. COMBINED HEATMAPS (3 LOCATIONS TOGETHER)
# =========================================================================

metrics <- list(
  "Degree" = "degree_norm",
  "Betweenness" = "betweenness_norm",
  "Closeness" = "closeness_norm",
  "Eigenvector" = "eigenvector_norm"
)

create_combined_centrality_heatmap <- function(data, metric_name, metric_col) {
  
  # Filter to important species and exclude combined 'AllLocations'
  combined_data <- data %>%
    filter(species %in% important_species, .data$location != "AllLocations")
  
  # Aggregate network contexts and prepare data for all locations (location-level means)
  heatmap_data <- combined_data %>%
    group_by(location, species, functional_group) %>%
    summarise(mean_value = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop")
  
  # Calculate species ordering within functional groups based on mean value
  species_order <- heatmap_data %>%
    group_by(species, functional_group) %>%
    summarise(avg_value = mean(mean_value, na.rm = TRUE), .groups = "drop") %>%
    arrange(functional_group, desc(avg_value)) %>%
    mutate(species_order = paste(functional_group, sprintf("%03d", row_number()), sep = "_"))
  
  # Join back to data
  heatmap_data <- heatmap_data %>%
    left_join(species_order %>% select(species, species_order), by = "species")
  
  # Metric-specific subtitles
  subtitle_text <- case_when(
    metric_name == "Degree" ~ "Number of direct contact partners (key for direct transmission)",
    metric_name == "Betweenness" ~ "Bridging role between groups (key for disease spread)",
    metric_name == "Closeness" ~ "Proximity to all species (reach within network)",
    metric_name == "Eigenvector" ~ "Connection to other central species (amplification potential)",
    TRUE ~ ""
  )
  
  # Create plot with a clean journal-style layout and location facets
  # Use location as the x-axis (one column per location)
  p <- ggplot(heatmap_data, aes(x = location, y = reorder(species, species_order), fill = mean_value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", mean_value)), 
              size = 4, color = "white", fontface = "plain") +
    
    # Color palette
    scale_fill_gradientn(
      colors = gradient_palette(100),
      name = "Normalized\nScore",
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c("Low", "Medium", "High"),
      na.value = "grey90"
    ) +
    scale_x_discrete(labels = location_labels) +
    
    # Facet by functional group (rows); location shown on the x-axis
    facet_grid(functional_group ~ ., scales = "free_y", space = "free_y") +
    
    labs(
      title = paste(metric_name, "centrality across locations"),
      subtitle = subtitle_text,
      x = "Location",
      y = NULL
    ) +
    
    theme_minimal(base_size = 14, base_family = "sans") +
    theme(
      # White background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Title styling - sentence case, bold
      plot.title = element_text(face = "bold", size = 16, hjust = 0,
                                margin = margin(b = 4)),
      plot.subtitle = element_text(size = 13, hjust = 0, color = "grey40", 
                                   margin = margin(b = 10)),
      
      # Axes
      axis.text.x = element_text(size = 13, face = "plain", angle = 0, hjust = 0.5,
                                color = "grey20"),
      axis.text.y = element_text(size = 13, face = "italic", color = "grey30"),
      axis.title.x = element_text(face = "plain", size = 14, margin = margin(t = 8),
                                 color = "grey20"),
      axis.title.y = element_text(face = "plain", size = 14, margin = margin(r = 8),
                                 color = "grey20"),
      
      # Facet strips - clean, minimal
      strip.text = element_text(face = "bold", size = 13, color = "grey20"),
      strip.text.y = element_text(angle = 0, hjust = 0),
      strip.background = element_rect(fill = "grey95", color = NA),
      
      # Grid and spacing
      panel.grid = element_blank(),
      panel.spacing = unit(0.3, "lines"),
      panel.border = element_rect(color = "grey80", linewidth = 0.3),
      
      # Legend
      legend.position = "right",
      legend.title = element_text(face = "plain", size = 13, color = "grey20"),
      legend.text = element_text(size = 12, color = "grey30"),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.4, "cm"),
      
      # Margins
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(p)
}

# Generate combined heatmaps for the 3 locations (exclude AllLocations)
centrality_df_filtered <- centrality_df %>%
  filter(location != "AllLocations")

for(metric_name in names(metrics)) {
  p <- create_combined_centrality_heatmap(centrality_df_filtered, metric_name, metrics[[metric_name]])
  
  if(!is.null(p)) {
    safe_metric <- tolower(metric_name)
    
    # Calculate height based on number of species
    n_species <- centrality_df_filtered %>%
      filter(species %in% important_species) %>%
      pull(species) %>%
      n_distinct()
    
    plot_height <- max(10, n_species * 0.35)
    plot_width <- 16
    
    ggsave(file.path(OUTPUT_DIR, paste0("heatmap_", safe_metric, "_3locations.png")),
           plot = p, 
           width = plot_width, 
           height = plot_height, 
           dpi = 600, 
           bg = "white")
    
  }
}
