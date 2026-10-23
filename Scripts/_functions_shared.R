# =========================================================================
# SHARED FUNCTIONS & COLOR PALETTES
# =========================================================================
# Source this file at the top of all analysis scripts to ensure consistency
# across guild assignments, color schemes, and helper functions.
# =========================================================================

# Load required libraries (silently, for sourcing)
suppressPackageStartupMessages({
  require(dplyr, quietly = TRUE)
  require(ggplot2, quietly = TRUE)
  require(data.table, quietly = TRUE)
  require(janitor, quietly = TRUE)
})


# COLOR PALETTES ------------------------------------------

# Sequential gradient palette: purple-navy → teal → yellow-green
color_palette <- function(n = 100) {
  colorRampPalette(c("#2A1B5D", "#1F4788", "#1E6B8C", "#2A9D8F", "#8FBC8F", "#E9C46A", "#F4E76E"))(n)
}

# Guild-specific categorical colors (used for species visualization)
guild_colors <- c(
  "Carnivore"    = "#D62828",  # Deep red - dominant predators
  "Herbivore"    = "#E9C46A",  # Gold/yellow - large herbivores
  "Omnivore"     = "#2A9D8F",  # Teal - generalists
  "Raptor"       = "#2A1B5D",  # Deep purple-navy - birds of prey
  "Corvid"       = "#8B4789",  # Purple - intelligent birds
  "Other bird"   = "#4A90C6",  # Blue - general avian
  "Reptile"      = "#8FBC8F"   # Sage green - ectotherms
)

# Risk-specific colors (for disease/threat visualization)
risk_colors <- c(
  low = "#8FBC8F",
  medium = "#F4A261",
  high = "#E76F51",
  veryhigh = "#D62828"
)

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
                        "Willie Wag-Tail", "White Face Heron", "White-Winged Chough",
                        "Black Faced Woodswallow", "Singing Honeyeater", "Indian Miner", 
                        "Magpie Lark", "White-Throated Tree Creeper", "Magpie", "Currawong") ~ "Other bird",
    
    # Reptiles (scavenging lizards)
    species_name %in% c("Lace Goanna", "Sand Goanna", "Bearded Dragon", 
                        "Millitary Dragon", "Central Netted Dragon", "Blue-Tongued Lizard") ~ "Reptile",
    
    # Default
    TRUE ~ "Other"
  )
}

# =========================================================================
# FUNCTIONS FOR PLOTTING
# =========================================================================

# Reorder factor within facets (preserves within-group ordering)
# Used for ggplot2 faceted plots where you want species ordered within each facet
reorder_within <- function(x, by, within, fun = mean, sep = "___") {
  new_x <- paste(as.character(x), within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

# Scale for removing the within-group suffix added by reorder_within
scale_x_reordered <- function(..., sep = "___") {
  ggplot2::scale_x_discrete(labels = function(x) gsub(paste0(sep, ".+$"), "", x), ...)
}

# =========================================================================
# DATA LOADING 
# =========================================================================

# Load and prepare network RDS file with standard cleaning
load_network_rds <- function(input_folder = "./Outputs/Networks/", time_window = "24h") {
  file_path <- file.path(input_folder, paste0("outputs_", time_window, ".rds"))
  if (!file.exists(file_path)) {
    stop("RDS file not found: ", file_path)
  }
  
  df <- readRDS(file_path)
  df <- janitor::clean_names(df)
  data.table::setDT(df)
  
  return(df)
}
