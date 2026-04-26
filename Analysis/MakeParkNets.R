# Purpose: Build park-level plant-pollinator networks for Broward County,
# calculate network metrics with bipartite, and write out a park-level
# table of network metrics for downstream analysis with honeybee abundance.
#
# Also makes a Broward study-area map where parks are colored by
# relative abundance of Apis mellifera.

# Libraries ####
library(tidyverse)
library(bipartite)
library(sf)
library(ggplot2)
library(tigris)
library(cowplot)
library(scales)

options(tigris_use_cache = TRUE)

# 1. Read data ####
data <- read.csv("../Data/interactions_data_4_4_2025.csv", stringsAsFactors = FALSE)

# Clean confidence column
# Original code had: FLW_ID_Conf != c(NA, "<80%"), which does not work as intended.
data <- data %>%
  filter(
    !is.na(FLW_ID_Conf),
    FLW_ID_Conf != "<80%"
  )

# Inspect names first
names(data)
glimpse(data)

# 2. Define the key columns ####
park_col  <- "Park.name"
plant_col <- "Flower_species"
poll_col  <- "Taxon.name"

# Define honeybee taxon name
apis_name <- "Apis mellifera"

# 3. Basic cleaning ####
dat_clean <- data %>%
  filter(
    !is.na(.data[[park_col]]),
    !is.na(.data[[plant_col]]),
    !is.na(.data[[poll_col]])
  ) %>%
  mutate(
    across(all_of(c(park_col, plant_col, poll_col)), trimws)
  ) %>%
  filter(
    .data[[park_col]] != "",
    .data[[plant_col]] != "",
    .data[[poll_col]] != ""
  )

# Optional: drop honeybees from network metrics if desired.
# Leave commented if you want honeybees included in network structure.
# dat_clean_no_apis <- dat_clean %>%
#   filter(.data[[poll_col]] != apis_name)

# 4. Function to make a bipartite interaction matrix for one park ####
make_web <- function(df, plant_col, poll_col, count_col = NULL) {
  
  if (is.null(count_col)) {
    web_df <- df %>%
      count(.data[[plant_col]], .data[[poll_col]], name = "freq")
  } else {
    web_df <- df %>%
      group_by(.data[[plant_col]], .data[[poll_col]]) %>%
      summarise(freq = sum(.data[[count_col]], na.rm = TRUE), .groups = "drop")
  }
  
  web <- web_df %>%
    pivot_wider(
      names_from = all_of(poll_col),
      values_from = freq,
      values_fill = 0
    ) %>%
    as.data.frame()
  
  rownames(web) <- web[[plant_col]]
  web[[plant_col]] <- NULL
  
  web <- as.matrix(web)
  
  return(web)
}

# 5. Function to calculate network metrics for one park ####
calc_metrics <- function(web) {
  
  # Require at least 2 plants and 2 pollinators for many metrics to behave well
  if (nrow(web) < 2 || ncol(web) < 2 || sum(web) == 0) {
    return(tibble(
      n_plants = nrow(web),
      n_pollinators = ncol(web),
      n_links = sum(web > 0),
      interaction_sum = sum(web),
      connectance = NA_real_,
      web_asymmetry = NA_real_,
      nestedness = NA_real_,
      modularity = NA_real_,
      H2 = NA_real_,
      linkage_density = NA_real_,
      interaction_evenness = NA_real_
    ))
  }
  
  out <- tryCatch({
    
    nl <- bipartite::networklevel(
      web,
      index = c(
        "connectance",
        "web asymmetry",
        "nestedness",
        "H2",
        "linkage density",
        "interaction evenness"
      ),
      ISAmethod = "Bluethgen"
    )
    
    # Modularity is separate / sometimes finicky
    mod <- tryCatch({
      bipartite::computeModules(web)@likelihood
    }, error = function(e) NA_real_)
    
    tibble(
      n_plants = nrow(web),
      n_pollinators = ncol(web),
      n_links = sum(web > 0),
      interaction_sum = sum(web),
      connectance = unname(nl["connectance"]),
      web_asymmetry = unname(nl["web asymmetry"]),
      nestedness = unname(nl["nestedness"]),
      modularity = mod,
      H2 = unname(nl["H2"]),
      linkage_density = unname(nl["linkage density"]),
      interaction_evenness = unname(nl["interaction evenness"])
    )
    
  }, error = function(e) {
    
    tibble(
      n_plants = nrow(web),
      n_pollinators = ncol(web),
      n_links = sum(web > 0),
      interaction_sum = sum(web),
      connectance = NA_real_,
      web_asymmetry = NA_real_,
      nestedness = NA_real_,
      modularity = NA_real_,
      H2 = NA_real_,
      linkage_density = NA_real_,
      interaction_evenness = NA_real_
    )
  })
  
  return(out)
}

# 6. Split into park-level networks and calculate metrics ####
park_list <- split(dat_clean, dat_clean[[park_col]])

park_metrics <- imap_dfr(park_list, function(df_park, park_name) {
  
  web <- make_web(
    df = df_park,
    plant_col = plant_col,
    poll_col = poll_col
  )
  
  metrics <- calc_metrics(web)
  
  metrics %>%
    mutate(
      park_name = park_name,
      n_records = nrow(df_park)
    ) %>%
    relocate(park_name)
})

# 7. Calculate Apis relative abundance per park ####
apis_per_park <- dat_clean %>%
  mutate(is_apis = .data[[poll_col]] == apis_name) %>%
  group_by(park_name = .data[[park_col]]) %>%
  summarise(
    total_interactions = n(),
    apis_interactions = sum(is_apis, na.rm = TRUE),
    apis_relative_abundance = apis_interactions / total_interactions,
    .groups = "drop"
  )

# Add Apis data to network metric table
park_metrics <- park_metrics %>%
  left_join(apis_per_park, by = "park_name")

# 8. View output ####
print(park_metrics)
glimpse(park_metrics)

# 9. Write output ####
write.csv(
  park_metrics,
  "../outputs/park_level_network_metrics.csv",
  row.names = FALSE
)

