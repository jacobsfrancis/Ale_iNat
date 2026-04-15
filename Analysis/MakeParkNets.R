# Purpose: Build park-level plant-pollinator networks for Broward County,
# calculate network metrics with bipartite, and write out a park-level
# table of network metrics for downstream analysis with honeybee abundance.

# Libraries ####
library(tidyverse)
library(bipartite)

# 1. Read data ####
data <- read.csv("../Data/interactions_data_4_4_2025.csv", stringsAsFactors = FALSE)
data<- data %>% filter(FLW_ID_Conf!=c(NA,"<80%"))
# Inspect names first
names(data)
glimpse(data)

# 2. Define the key columns ####
# CHANGE THESE to match your actual column names
park_col      <- "Park.name"
plant_col     <- "Flower_species"
poll_col      <- "Taxon.name"



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

# Optional: drop honeybees, casual IDs, unknowns, etc. here if needed
# dat_clean <- dat_clean %>%
#    filter(Taxon.name != "Apis mellifera")

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
  
  # networklevel returns a named vector/list depending on metric
  out <- tryCatch({
    nl <- networklevel(
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
      computeModules(web)@likelihood
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

# 7. View output ####
print(park_metrics)
glimpse(park_metrics)

# 8. Write output ####
write.csv(
  park_metrics,
  "../Output/park_level_network_metrics.csv",
  row.names = FALSE
)
