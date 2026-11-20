# Purpose: Test whether plant-pollinator interaction composition varies
# with park-level honey bee relative abundance.
#
# Question:
#   Do parks with higher Apis mellifera relative abundance have different
#   plant-pollinator interaction composition?
#
# Approach:
#   Treat each plant x pollinator pair as an "interaction species".
#   Build a park x interaction matrix.
#   Test composition using Bray-Curtis dissimilarity and PERMANOVA.
#
# This script runs two versions:
#   1. Weighted/count composition:
#        cells = number of observations of each plant x pollinator interaction
#   2. Binary/presence-absence composition:
#        cells = 1 if interaction was observed in a park, 0 otherwise
#
# Filter:
#   Only parks with at least 10 unique plant-pollinator links are included.
#
# Inputs:
#   ../Data/interactions_data_4_4_2025.csv
#   ../Output/ApisRelAbund.csv
#   ../Output/park_level_network_metrics.csv
#
# Outputs:
#   ../Output/interaction_composition_matrix.csv
#   ../Output/interaction_composition_matrix_filtered_ge10_links_weighted.csv
#   ../Output/interaction_composition_matrix_filtered_ge10_links_binary.csv
#   ../Output/interaction_composition_metadata_filtered_ge10_links.csv
#   ../Output/interaction_composition_permanova_weighted.csv
#   ../Output/interaction_composition_permanova_binary.csv
#   ../Output/interaction_composition_dispersion_test_weighted.txt
#   ../Output/interaction_composition_dispersion_test_binary.txt
#   ../Output/interaction_composition_envfit_weighted.txt
#   ../Output/interaction_composition_envfit_binary.txt
#   ../Output/interaction_composition_envfit_vectors_weighted.csv
#   ../Output/interaction_composition_envfit_vectors_binary.csv
#   ../Output/interaction_composition_nmds_weighted_poster.png
#   ../Output/interaction_composition_nmds_binary_poster.png
#   ../Output/interaction_composition_links_correlated_with_honeybee_weighted.csv
#   ../Output/interaction_composition_links_correlated_with_honeybee_binary.csv

# 0. Packages ####

library(tidyverse)
library(vegan)
library(broom)
library(stringr)

# 0.1 Plot style and colors ####

fau_blue <- "#003366"
fau_red <- "#CC0000"
fau_gray <- "#CCCCCC"
fau_dark_gray <- "#4D4C55"

hb_low <- "#E8E8E8"
hb_mid <- "#E6A3A3"
hb_high <- fau_red
hb_ramp <- c(hb_low, hb_mid, hb_high)

theme_hb <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
      axis.title = element_text(color = "grey15"),
      axis.text = element_text(color = "grey25"),
      plot.title = element_text(face = "bold", color = "grey10"),
      plot.subtitle = element_text(color = "grey25"),
      legend.title = element_text(color = "grey15"),
      legend.text = element_text(color = "grey25")
    )
}

# 1. Read data ####

interactions <- read.csv(
  "../Data/interactions_data_4_4_2025.csv",
  stringsAsFactors = FALSE
)

relabund <- read.csv(
  "../Output/ApisRelAbund.csv",
  stringsAsFactors = FALSE
)

netMets <- read.csv(
  "../Output/park_level_network_metrics.csv",
  stringsAsFactors = FALSE
)

# 2. Define columns ####

park_col <- "Park.name"
plant_genus_col <- "Flower_Genus"
plant_species_col <- "Flower_species"
poll_col <- "Taxon.name"

apis_name <- "Apis mellifera"

# 3. Clean interaction data ####

interactions_clean <- interactions %>%
  filter(
    !is.na(FLW_ID_Conf),
    FLW_ID_Conf != "<80%",
    !is.na(.data[[park_col]]),
    !is.na(.data[[plant_genus_col]]),
    !is.na(.data[[plant_species_col]]),
    !is.na(.data[[poll_col]])
  ) %>%
  mutate(
    across(
      all_of(c(park_col, plant_genus_col, plant_species_col, poll_col)),
      ~ str_squish(as.character(.x))
    )
  ) %>%
  filter(
    .data[[park_col]] != "",
    .data[[plant_genus_col]] != "",
    .data[[plant_species_col]] != "",
    .data[[poll_col]] != ""
  ) %>%
  mutate(
    park_name = .data[[park_col]],
    park_join = str_trim(str_to_lower(park_name)),
    plant_taxon = str_squish(
      paste(.data[[plant_genus_col]], .data[[plant_species_col]])
    ),
    pollinator_taxon = .data[[poll_col]],
    interaction_id = paste(plant_taxon, pollinator_taxon, sep = " __ ")
  )

# Quick check of interaction IDs
interactions_clean %>%
  distinct(interaction_id) %>%
  arrange(interaction_id)

# 4. Build park x interaction matrix ####
# Each column is a plant-pollinator interaction.
# Each row is a park.
# Values are counts of that interaction in that park.

interaction_matrix_df <- interactions_clean %>%
  count(park_join, interaction_id, name = "interaction_count") %>%
  pivot_wider(
    names_from = interaction_id,
    values_from = interaction_count,
    values_fill = 0
  )

interaction_matrix <- interaction_matrix_df %>%
  column_to_rownames("park_join") %>%
  as.matrix()

# Drop parks with zero interactions, just in case
interaction_matrix <- interaction_matrix[
  rowSums(interaction_matrix) > 0,
  ,
  drop = FALSE
]

# Drop interactions occurring in zero parks, just in case
interaction_matrix <- interaction_matrix[
  ,
  colSums(interaction_matrix) > 0,
  drop = FALSE
]

cat("\n--- INTERACTION MATRIX SUMMARY BEFORE PARK FILTER ---\n")
cat("Parks:", nrow(interaction_matrix), "\n")
cat("Unique plant-pollinator interactions:", ncol(interaction_matrix), "\n")

write.csv(
  as.data.frame(interaction_matrix) %>%
    rownames_to_column("park_join"),
  "../Output/interaction_composition_matrix.csv",
  row.names = FALSE
)

# 5. Prepare park-level metadata ####

relabund_clean <- relabund %>%
  mutate(
    park_join = str_trim(str_to_lower(as.character(poly_id))),
    hb_percent = as.numeric(percentageHB),
    hb_raw = hb_percent / 100
  ) %>%
  select(
    park_join,
    n_apidae,
    n_hb,
    hb_percent,
    hb_raw
  )

netMets_clean <- netMets %>%
  mutate(
    park_join = str_trim(str_to_lower(as.character(park_name))),
    interaction_sum = as.numeric(interaction_sum),
    log_interaction_sum = log1p(interaction_sum)
  ) %>%
  select(
    park_join,
    interaction_sum,
    log_interaction_sum,
    n_plants,
    n_pollinators
  )

park_labels <- interactions_clean %>%
  distinct(park_join, park_name)

meta <- tibble(
  park_join = rownames(interaction_matrix)
) %>%
  left_join(park_labels, by = "park_join") %>%
  left_join(relabund_clean, by = "park_join") %>%
  left_join(netMets_clean, by = "park_join")

cat("\n--- METADATA SUMMARY BEFORE PARK FILTER ---\n")
print(
  meta %>%
    summarise(
      n_parks = n(),
      missing_hb = sum(is.na(hb_raw)),
      missing_interaction_sum = sum(is.na(log_interaction_sum))
    )
)

cat("\nParks missing honey bee relative abundance:\n")
print(
  meta %>%
    filter(is.na(hb_raw)) %>%
    select(park_name, park_join)
)

cat("\nParks missing interaction_sum:\n")
print(
  meta %>%
    filter(is.na(log_interaction_sum)) %>%
    select(park_name, park_join)
)

# 6. Filter to sufficiently sampled networks ####
# Keep only parks with at least 10 unique plant-pollinator links.
# Here, links = unique plant x pollinator interaction pairs observed in that park.

park_link_counts <- interactions_clean %>%
  distinct(park_join, interaction_id) %>%
  count(park_join, name = "n_unique_links")

meta <- meta %>%
  left_join(park_link_counts, by = "park_join")

cat("\n--- LINK COUNT SUMMARY ---\n")
print(summary(meta$n_unique_links))

cat("\nParks with fewer than 10 unique links:\n")
print(
  meta %>%
    filter(is.na(n_unique_links) | n_unique_links < 10) %>%
    select(park_name, park_join, n_unique_links, hb_percent, interaction_sum)
)

# 7. Restrict to parks with complete metadata and at least 10 links ####

complete_parks <- meta %>%
  filter(
    !is.na(hb_raw),
    !is.na(log_interaction_sum),
    !is.na(n_unique_links),
    n_unique_links >= 10
  ) %>%
  pull(park_join)

interaction_matrix_complete <- interaction_matrix[
  complete_parks,
  ,
  drop = FALSE
]

meta_complete <- meta %>%
  filter(park_join %in% complete_parks) %>%
  arrange(match(park_join, rownames(interaction_matrix_complete)))

# Confirm row order
stopifnot(all(meta_complete$park_join == rownames(interaction_matrix_complete)))

cat("\n--- COMPLETE DATA FOR COMPOSITION ANALYSIS ---\n")
cat("Parks with >=10 unique links:", nrow(interaction_matrix_complete), "\n")
cat("Interactions before dropping absent columns:", ncol(interaction_matrix_complete), "\n")

# 8. Remove interactions absent after park filtering ####

interaction_matrix_complete <- interaction_matrix_complete[
  ,
  colSums(interaction_matrix_complete) > 0,
  drop = FALSE
]

cat("Interactions after dropping absent columns:", ncol(interaction_matrix_complete), "\n")

# Weighted/count matrix
interaction_matrix_weighted <- interaction_matrix_complete

# Binary/presence-absence matrix
interaction_matrix_binary <- (interaction_matrix_complete > 0) * 1

write.csv(
  as.data.frame(interaction_matrix_weighted) %>%
    rownames_to_column("park_join"),
  "../Output/interaction_composition_matrix_filtered_ge10_links_weighted.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(interaction_matrix_binary) %>%
    rownames_to_column("park_join"),
  "../Output/interaction_composition_matrix_filtered_ge10_links_binary.csv",
  row.names = FALSE
)

write.csv(
  meta_complete,
  "../Output/interaction_composition_metadata_filtered_ge10_links.csv",
  row.names = FALSE
)

# 9. Helper function to run composition analysis ####

run_interaction_composition_analysis <- function(
    matrix_input,
    meta_input,
    analysis_label = "weighted",
    matrix_type_label = "Weighted interaction composition"
) {
  
  cat("\n========================================\n")
  cat("RUNNING ANALYSIS:", analysis_label, "\n")
  cat("========================================\n")
  
  # 9.1 PERMANOVA ####
  
  set.seed(123)
  
  permanova_interactions <- vegan::adonis2(
    matrix_input ~ hb_raw + log_interaction_sum,
    data = meta_input,
    method = "bray",
    permutations = 999,
    by = "margin"
  )
  
  permanova_results <- as.data.frame(permanova_interactions) %>%
    rownames_to_column("term")
  
  cat("\n--- PERMANOVA:", analysis_label, "---\n")
  print(permanova_results)
  
  write.csv(
    permanova_results,
    paste0("../Output/interaction_composition_permanova_", analysis_label, ".csv"),
    row.names = FALSE
  )
  
  # Extract honey bee PERMANOVA P-value if available
  hb_p <- permanova_results %>%
    filter(term == "hb_raw") %>%
    pull(`Pr(>F)`)
  
  if (length(hb_p) == 0) {
    hb_p <- NA_real_
  }
  
  # 9.2 Dispersion check ####
  
  meta_for_disp <- meta_input %>%
    mutate(
      hb_group = if_else(
        hb_raw >= median(hb_raw, na.rm = TRUE),
        "High honey bee",
        "Low honey bee"
      ),
      hb_group = factor(
        hb_group,
        levels = c("Low honey bee", "High honey bee")
      )
    )
  
  bray_dist <- vegdist(matrix_input, method = "bray")
  
  disp <- betadisper(bray_dist, meta_for_disp$hb_group)
  disp_test <- permutest(disp, permutations = 999)
  
  cat("\n--- DISPERSION TEST:", analysis_label, "---\n")
  print(disp_test)
  
  capture.output(
    disp_test,
    file = paste0("../Output/interaction_composition_dispersion_test_", analysis_label, ".txt")
  )
  
  # 9.3 NMDS ordination ####
  
  set.seed(123)
  
  nmds <- metaMDS(
    matrix_input,
    distance = "bray",
    k = 2,
    trymax = 100,
    autotransform = FALSE,
    trace = FALSE
  )
  
  cat("\n--- NMDS STRESS:", analysis_label, "---\n")
  print(nmds$stress)
  
  nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
    rownames_to_column("park_join") %>%
    left_join(meta_for_disp, by = "park_join")
  
  # 9.4 Envfit vectors ####
  
  envfit_hb <- envfit(
    nmds,
    meta_for_disp %>%
      select(hb_raw, log_interaction_sum),
    permutations = 999
  )
  
  cat("\n--- ENVFIT:", analysis_label, "---\n")
  print(envfit_hb)
  
  capture.output(
    envfit_hb,
    file = paste0("../Output/interaction_composition_envfit_", analysis_label, ".txt")
  )
  
  envfit_vectors <- as.data.frame(scores(envfit_hb, display = "vectors")) %>%
    rownames_to_column("variable") %>%
    mutate(
      r2 = envfit_hb$vectors$r,
      p.value = envfit_hb$vectors$pvals,
      label = case_when(
        variable == "hb_raw" ~ "% honey bee",
        variable == "log_interaction_sum" ~ "Sampling effort",
        TRUE ~ variable
      )
    )
  
  arrow_scale <- 0.75 * min(
    diff(range(nmds_scores$NMDS1, na.rm = TRUE)),
    diff(range(nmds_scores$NMDS2, na.rm = TRUE))
  )
  
  envfit_vectors <- envfit_vectors %>%
    mutate(
      NMDS1_end = NMDS1 * arrow_scale,
      NMDS2_end = NMDS2 * arrow_scale,
      label_x = NMDS1_end * 1.12,
      label_y = NMDS2_end * 1.12,
      arrow_color = case_when(
        label == "% honey bee" ~ fau_red,
        label == "Sampling effort" ~ fau_dark_gray,
        TRUE ~ fau_dark_gray
      )
    )
  
  print(envfit_vectors)
  
  write.csv(
    envfit_vectors,
    paste0("../Output/interaction_composition_envfit_vectors_", analysis_label, ".csv"),
    row.names = FALSE
  )
  
  # 9.5 Poster NMDS plot ####
  
  nmds_plot <- ggplot(
    nmds_scores,
    aes(x = NMDS1, y = NMDS2)
  ) +
    geom_hline(
      yintercept = 0,
      color = "grey90",
      linewidth = 0.35
    ) +
    geom_vline(
      xintercept = 0,
      color = "grey90",
      linewidth = 0.35
    ) +
    geom_point(
      aes(
        fill = hb_percent,
        size = interaction_sum
      ),
      shape = 21,
      color = fau_dark_gray,
      stroke = 0.6,
      alpha = 0.85
    ) +
    geom_segment(
      data = envfit_vectors,
      aes(
        x = 0,
        y = 0,
        xend = NMDS1_end,
        yend = NMDS2_end,
        color = label
      ),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.18, "cm")),
      linewidth = 0.8
    ) +
    geom_text(
      data = envfit_vectors,
      aes(
        x = label_x,
        y = label_y,
        label = label,
        color = label
      ),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4
    ) +
    scale_color_manual(
      values = c(
        "% honey bee" = fau_red,
        "Sampling effort" = fau_dark_gray
      ),
      guide = "none"
    ) +
    scale_fill_gradientn(
      colors = hb_ramp,
      limits = c(0, 100),
      guide = "none"
    ) +
    scale_size_continuous(
      range = c(2.5, 7),
      guide = "none"
    ) +
    labs(
      title = "Interaction composition",
      subtitle = paste0(
        matrix_type_label,
        "; honey bee P = ",
        signif(hb_p, 2),
        "; stress = ",
        round(nmds$stress, 2)
      ),
      x = "NMDS1",
      y = "NMDS2"
    ) +
    theme_hb(base_size = 20) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 24, face = "plain"),
      plot.subtitle = element_text(size = 20),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 25, 10, 10)
    )
  
  print(nmds_plot)
  
  ggsave(
    paste0("../Output/interaction_composition_nmds_", analysis_label, "_poster.png"),
    nmds_plot,
    width = 7.5,
    height = 6,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  # 9.6 Exploratory interaction-level correlations ####
  
  interaction_correlations <- map_dfr(
    colnames(matrix_input),
    function(interaction_id) {
      x <- matrix_input[, interaction_id]
      
      if (sd(x, na.rm = TRUE) == 0) {
        return(tibble())
      }
      
      ct <- suppressWarnings(
        cor.test(x, meta_for_disp$hb_raw, method = "spearman")
      )
      
      tibble(
        interaction_id = interaction_id,
        rho = unname(ct$estimate),
        p.value = ct$p.value,
        total_count_or_presence = sum(x),
        n_parks_present = sum(x > 0)
      )
    }
  ) %>%
    separate(
      interaction_id,
      into = c("plant_taxon", "pollinator_taxon"),
      sep = " __ ",
      remove = FALSE
    ) %>%
    mutate(
      p_adj_fdr = p.adjust(p.value, method = "fdr")
    ) %>%
    arrange(p.value)
  
  write.csv(
    interaction_correlations,
    paste0("../Output/interaction_composition_links_correlated_with_honeybee_", analysis_label, ".csv"),
    row.names = FALSE
  )
  
  cat("\n--- TOP INTERACTIONS ASSOCIATED WITH HONEY BEE RELATIVE ABUNDANCE:", analysis_label, "---\n")
  print(interaction_correlations %>% slice_head(n = 20))
  
  # Return objects in case you want to inspect them interactively
  return(
    list(
      permanova = permanova_results,
      dispersion = disp_test,
      nmds = nmds,
      nmds_scores = nmds_scores,
      envfit = envfit_hb,
      envfit_vectors = envfit_vectors,
      plot = nmds_plot,
      interaction_correlations = interaction_correlations
    )
  )
}

# 10. Run weighted/count analysis ####
# This asks whether parks differ in the composition and relative frequency
# of plant-pollinator interactions.

weighted_results <- run_interaction_composition_analysis(
  matrix_input = interaction_matrix_weighted,
  meta_input = meta_complete,
  analysis_label = "weighted",
  matrix_type_label = "Weighted/count composition"
)

# 11. Run binary/presence-absence analysis ####
# This asks whether parks differ in which interaction links are present,
# ignoring how frequently each link was observed.

binary_results <- run_interaction_composition_analysis(
  matrix_input = interaction_matrix_binary,
  meta_input = meta_complete,
  analysis_label = "binary",
  matrix_type_label = "Presence-absence composition"
)

# 12. Print brief comparison ####

cat("\n========================================\n")
cat("WEIGHTED VS BINARY SUMMARY\n")
cat("========================================\n")

cat("\nWeighted PERMANOVA:\n")
print(weighted_results$permanova)

cat("\nBinary PERMANOVA:\n")
print(binary_results$permanova)

cat("\nWeighted NMDS stress:", weighted_results$nmds$stress, "\n")
cat("Binary NMDS stress  :", binary_results$nmds$stress, "\n")