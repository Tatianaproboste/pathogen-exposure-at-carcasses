# Creates UNDIRECTED: Co-occurrence (who feeds together)
# =========================================================================

# LIBRARIES ---------------------------------------------------------------
invisible(lapply(c("data.table", "igraph", "tidygraph",
 "ggraph", "spatsoc", "asnipe", "tidyverse", "ggrepel"), 
 library, character.only = TRUE))

# SOURCE SHARED FUNCTIONS -------------------------------------------------
# Keep shared objects (e.g., guild_colors) in the current script environment
# so they are available even when this script is sourced with local = TRUE.
source("Scripts/_functions_shared.R", local = TRUE)

# Location labels: short code → full region name used in all plot titles
location_labels <- c(KNP = "Alpine", SD = "Desert", WV = "Temperate")

# SETTINGS ----------------------------------------------------------------
INPUT_FOLDER <- "./Outputs/Networks/"
OUTPUT_FOLDER <- "./Outputs/Networks/Structure_Analysis/"
dir.create(OUTPUT_FOLDER, showWarnings = FALSE, recursive = TRUE)

TIME_WINDOW <- "24h"   # Options: "24h", "48h", "72h"
MIN_SPECIES <- 3       # Skip contexts with fewer species
MIN_OBS <- 10          # Skip contexts with fewer observations

file_path <- file.path(INPUT_FOLDER, paste0("outputs_", TIME_WINDOW, ".rds"))
df <- readRDS(file_path)
df <- janitor::clean_names(df)
# Ensure df is a data.table so data.table special symbols like .N work
data.table::setDT(df)

# =========================================================================
# INDEPENDENT FEEDING EVENTS
# =========================================================================
# An "independent feeding event" = one station visited within one time
# window (the `group` column created by spatsoc::group_times() in
# 01_Network_construction.R). Records sharing the same `group` value were
# observed at the same carcass within the same TIME_WINDOW.

feeding_events_by_location <- df[, list(
  n_feeding_events = data.table::uniqueN(group),
  n_observations   = .N,
  n_species        = data.table::uniqueN(species),
  n_stations       = data.table::uniqueN(station)
), by = location]

feeding_events_total <- df[, list(
  n_feeding_events = data.table::uniqueN(group),
  n_observations   = .N,
  n_species        = data.table::uniqueN(species),
  n_stations       = data.table::uniqueN(station)
)]

cat("\nIndependent feeding events (", TIME_WINDOW, " window):\n", sep = "")
print(feeding_events_by_location)
cat("Total across all locations:", feeding_events_total$n_feeding_events, "\n\n")

fwrite(feeding_events_by_location, file.path(OUTPUT_FOLDER, "feeding_events_by_location_V2.csv"))

# =========================================================================
# GUILD CLASSIFICATION
# =========================================================================

assign_guild <- function(species_name) {  
  case_when(
    # Carnivores (mammalian predators)
    species_name %in% c("Dingo", "Red Fox", "Cat") ~ "Carnivore",
    # Herbivores (large herbivorous mammals)
    species_name %in% c("Red Kangaroo", "Eastern Grey Kangaroo", "Deer", "Camel") ~ "Herbivore",
   # Omnivores (opportunistic feeders)
    species_name %in% c("Feral Pig", "Brushtail Possum") ~ "Omnivore",
    # Raptors (birds of prey)
    species_name %in% c("Wedge-Tailed Eagle", "Little Eagle", "Black Kite", 
                       "Brown Falcon", "Nankeen Kestrel", "Nankeen", "Grey Goshawk") ~ "Raptor",
   # Corvids (intelligent scavenger birds)
    species_name %in% c("Corvid") ~ "Corvid",
    # Other birds (generalist/opportunistic birds)
    species_name %in% c("Kookaburra", "Grey Butcherbird", "Pied Butcherbird", 
                       "Willie Wag-Tail", "White Face Heron", "White-Winged Chough","Black Faced Woodswallow", 
                       "Singing Honeyeater", "Indian Miner", "Magpie Lark", "White-Throated Tree Creeper", "Magpie", "Currawong") ~ "Other bird",
    # Reptiles (scavenging lizards)
    species_name %in% c("Lace Goanna", "Sand Goanna", "Bearded Dragon", 
                       "Millitary Dragon", "Central Netted Dragon", "Blue-Tongued Lizard") ~ "Reptile",
    # Default
    TRUE ~ "Other"
  )
}

# =========================================================================
# STEP 2: CALCULATE NETWORK METRICS
# =========================================================================

calculate_metrics <- function(df, context_name) {
  # Build GBI matrix (groups × species)
   
  gbi <- tryCatch({
    get_gbi(df, group = 'group', id = 'species')
  }, error = function(e) {
    cat("    ⚠ Error creating GBI\n")
    return(NULL)
  })
  
  if (is.null(gbi) || ncol(gbi) < MIN_SPECIES) {
    cat(paste0("    ⚠ Skipped: only ", ncol(gbi), " species\n"))
    return(NULL)
  }
  
  cat(paste0("    GBI: ", nrow(gbi), " carcasses × ", ncol(gbi), " species\n"))
  
  # -----------------------------------------------------------------------
  # NETWORK 1: Co-occurrence (UNDIRECTED)
  # -----------------------------------------------------------------------
  
    # Calculate HWI (association strength)
  hwi_matrix <- get_network(gbi, data_format = "GBI", association_index = "HWI")
  
  # Create network
  net_cooccur <- graph_from_adjacency_matrix(
    hwi_matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  # -----------------------------------------------------------------------
  # GLOBAL METRICS (network-level)
  # -----------------------------------------------------------------------
  
  global <- data.table::data.table(
    context = context_name,
    
  # Structure
  n_stations = data.table::uniqueN(df$station),  # Number of unique carcasses/stations
  n_carcass_days = nrow(gbi),  # Number of carcass × time window combinations
    n_species = vcount(net_cooccur),
    n_edges = ecount(net_cooccur),
    density = edge_density(net_cooccur),
    
    # Connectivity
    n_components = components(net_cooccur)$no,
    largest_comp = max(components(net_cooccur)$csize),
    
    # Centralization (hub dominance)
    central_degree = centr_degree(net_cooccur)$centralization,
    central_between = centr_betw(net_cooccur)$centralization,
    
    # Clustering
    clustering = transitivity(net_cooccur, type = "global"),
    
    # Assortativity (do similar species connect?)
    # Positive = similar species interact, Negative = dissimilar species interact
    assortativity_degree = tryCatch(
      assortativity_degree(net_cooccur), 
      error = function(e) NA
    ),
    
    # Path lengths
    diameter = tryCatch(diameter(net_cooccur), error = function(e) NA),
    avg_path = tryCatch(mean_distance(net_cooccur), error = function(e) NA),
    
    # Degree
    mean_degree = mean(degree(net_cooccur)),
    max_degree = max(degree(net_cooccur))
  )
  
  # -----------------------------------------------------------------------
  # NODE METRICS (species-level)
  # -----------------------------------------------------------------------
  
  nodes <- data.table(
    context = context_name,
    species = V(net_cooccur)$name,
    
    # Connectivity
    degree = degree(net_cooccur),
    strength = strength(net_cooccur),  # Weighted degree
    
    # Centrality
    betweenness = betweenness(net_cooccur, normalized = TRUE),
    closeness = closeness(net_cooccur, normalized = TRUE),
    eigenvector = eigen_centrality(net_cooccur)$vector,
    
    # Clustering
    clustering = transitivity(net_cooccur, type = "local", isolates = "zero")
  )
  
  # Add temporal metrics if available (safe check if variable exists)
  if (exists("temporal_stats") && !is.null(temporal_stats)) {
    nodes <- merge(nodes, temporal_stats, by = c("context", "species"), all.x = TRUE)
  }
  
  # -----------------------------------------------------------------------
  # RETURN
  # -----------------------------------------------------------------------
    return(list(
    global = global,
    nodes = nodes,
    net_cooccur = net_cooccur 
  ))
}

# =========================================================================
# STEP 3: VISUALIZE NETWORKS
# =========================================================================

plot_network <- function(net, context, save_path, type = "cooccur", panel_label = NULL) {
  
  # Prepare for plotting and add guild classification
  tg <- as_tbl_graph(net) %>%
    activate("nodes") %>%
    mutate(
      deg = degree(net),
      guild = assign_guild(name)
    )
  
  if (type == "temporal") {
    # DIRECTED network (arrows show early → late) with a clean journal-style layout
    p <- ggraph(tg, layout = "fr") +
      # Edges with arrows - subtle grey
      geom_edge_link(
        aes(width = weight, alpha = weight),
        arrow = arrow(length = unit(3, 'mm'), type = "closed"),
        end_cap = circle(4, 'mm'),
        color = "grey50",
        show.legend = FALSE
      ) +
      # Nodes with guild colors - cleaner strokes
      geom_node_point(aes(size = deg, fill = guild), 
                      shape = 21, color = "grey30", stroke = 0.5, alpha = 0.9) +
      # Labels - cleaner style
      geom_node_text(aes(label = name), repel = TRUE, size = 3.5, 
                     fontface = "italic", color = "grey20",
                     bg.color = "white", bg.r = 0.12, segment.color = "grey70") +
      # Scales
      scale_edge_width(range = c(0.3, 2.0), guide = "none") +
      scale_edge_alpha(range = c(0.3, 0.7), guide = "none") +
      scale_size_continuous(range = c(3, 14), name = "Network degree",
                           breaks = c(4, 8, 12)) +
      scale_fill_manual(values = guild_colors, name = "Functional group") +
      # Labels - clear, clean style
      labs(
        title = paste("Temporal precedence network of scavenger arrivals —", context),
        subtitle = paste(vcount(net), "species •", ecount(net), "edges"),
        tag = panel_label
      ) +
      # Theme - clean journal-style layout
      theme_graph(base_family = "sans", base_size = 14) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(face = "bold", hjust = 0, size = 16, 
                                  margin = margin(b = 4)),
        plot.subtitle = element_text(hjust = 0, size = 13, color = "grey40", 
                                     margin = margin(b = 15)),
        plot.tag = element_text(face = "bold", size = 18, hjust = 0),
        plot.tag.position = c(0.02, 0.91),  # Left side, below title
        legend.position = "right",
        legend.title = element_text(face = "plain", size = 13, color = "grey20"),
        legend.text = element_text(size = 12, color = "grey30"),
        legend.key.size = unit(0.6, "cm"),
        legend.spacing.y = unit(0.2, "cm"),
        plot.margin = margin(8, 8, 8, 8)
      ) +
      guides(
        fill = guide_legend(override.aes = list(size = 6), order = 1),
        size = guide_legend(order = 2)
      )
    
  } else {
    # UNDIRECTED network (co-occurrence) with a clean journal-style layout
    p <- ggraph(tg, layout = "fr") +
      # Edges - subtle grey with transparency
      geom_edge_link(aes(width = weight, alpha = weight), 
                     color = "grey60", show.legend = FALSE) +
      # Nodes with guild colors - cleaner appearance
      geom_node_point(aes(size = deg, fill = guild), 
                      shape = 21, color = "grey30", stroke = 0.5, alpha = 0.9) +
      # Labels - refined typography
      geom_node_text(aes(label = name), repel = TRUE, size = 3.5, force = 4,
                    force_pull = 0.2,  box.padding = 0.8,  point.padding = 1.2,
                     fontface = "italic", color = "grey20", bg.color = "white", bg.r = 0.05, segment.color = "grey70") +
      # Scales
      scale_edge_width(range = c(0.4, 2.5), guide = "none") +
      scale_edge_alpha(range = c(0.3, 0.8), guide = "none") +
      scale_size_continuous(range = c(3, 14), name = "Network degree",
                           breaks = c(4, 8, 12, 16)) +
      scale_fill_manual(values = guild_colors, name = "Functional group") +
      # Labels - sentence case, clean styling
      labs(
        title = paste("Co-occurrence network of scavenger feeding associations —", context),
        subtitle = paste(vcount(net), "species •", ecount(net), "edges"),
        tag = panel_label
      ) +
      # Theme - clean journal-style layout
      theme_graph(base_family = "sans", base_size = 14) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(face = "bold", hjust = 0, size = 16, 
                                  margin = margin(b = 4)),
        plot.subtitle = element_text(hjust = 0, size = 13, color = "grey40", 
                                     margin = margin(b = 15)),
        plot.tag = element_text(face = "bold", size = 18, hjust = 0),
        plot.tag.position = c(0.02, 0.91),  # Left side, below title
        legend.position = "right",
        legend.title = element_text(face = "plain", size = 13, color = "grey20"),
        legend.text = element_text(size = 12, color = "grey30"),
        legend.key.size = unit(0.6, "cm"),
        legend.spacing.y = unit(0.2, "cm"),
        plot.margin = margin(8, 8, 8, 8)
      ) +
      guides(
        fill = guide_legend(override.aes = list(size = 5), order = 1),
        size = guide_legend(order = 2)
      )
  }
  
  ggsave(save_path, p, width = 14, height = 11, dpi = 600, bg = "white")
  return(p)
}

# =========================================================================
# STEP 4: IDENTIFY CONTEXTS
# =========================================================================

contexts <- df[, list(
  n_obs = .N,
  n_species = data.table::uniqueN(species),
  n_groups = data.table::uniqueN(group)
), by = "location"] # one context per location (no year/season)

# =========================================================================
# STEP 5: ANALYZE EACH CONTEXT
# =========================================================================

all_global <- list()
all_nodes <- list()
all_nets_cooccur <- list()

for (i in seq_len(nrow(contexts))) {

  # Get context (location only)
  loc <- contexts$location[i]

  context_name <- paste(TIME_WINDOW, loc, sep = "_")

  # Filter data for this location (all years/seasons)
  df_subset <- df[location == loc]
  
  # Calculate metrics
  result <- calculate_metrics(df_subset, context_name)
  
  if (!is.null(result)) {
    # Store results
    all_global[[context_name]] <- result$global
    all_nodes[[context_name]] <- result$nodes
    all_nets_cooccur[[context_name]] <- result$net_cooccur
  
    # No seasonal panel labels for location-level networks
    panel_label <- NULL

    # Display label (location only) - remove TIME_WINDOW from figure titles
    display_label <- ifelse(loc %in% names(location_labels), location_labels[[loc]], loc)

    # Plot co-occurrence network (versioned) - filename uses location only
    plot_path <- file.path(OUTPUT_FOLDER, paste0("cooccur_", loc, "_V2.png"))
    plot_network(result$net_cooccur, display_label, plot_path, type = "cooccur", panel_label = panel_label)
  }
 }

# =========================================================================
# STEP 8: COMBINE RESULTS
# =========================================================================

# Combine all global metrics (EXCLUDING overview networks)
global_metrics <- rbindlist(all_global[!grepl("_overview$", names(all_global))])

# Parse context into components
global_metrics <- global_metrics %>%
  separate(context, into = c("time_window", "location"), 
           sep = "_", remove = FALSE, fill = "right") %>%
  as.data.table()

# Combine all node metrics (EXCLUDING overview networks)
node_metrics <- rbindlist(all_nodes[!grepl("_overview$", names(all_nodes))])

# =========================================================================
# STEP 9: SAVE RESULTS
# =========================================================================

# Save data
fwrite(global_metrics, file.path(OUTPUT_FOLDER, "global_metrics_V2.csv"))
fwrite(node_metrics, file.path(OUTPUT_FOLDER, "node_metrics_V2.csv"))

# Save networks
saveRDS(list(
  cooccurrence = all_nets_cooccur
), file.path(OUTPUT_FOLDER, "networks_V2.rds"))

# Summary
summary <- global_metrics[, list(
  mean_density = mean(density),
  mean_clustering = mean(clustering),
  mean_centralization = mean(central_degree),
  mean_assortativity = mean(assortativity_degree, na.rm = TRUE)
), by = "location"]

fwrite(summary, file.path(OUTPUT_FOLDER, "summary_V2.csv"))

