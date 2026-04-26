# Purpose: Test whether non-Apis pollinator specialization varies with
# park-level honey bee relative abundance.
#
# Question:
#   Are non-Apis pollinators more specialized in parks where Apis mellifera
#   makes up a larger proportion of bee observations?
#
# Inputs:
#   ../Data/interactions_data_4_4_2025.csv
#   ../Output/ApisRelAbund.csv
#   ../Output/park_level_network_metrics.csv
#
# Outputs:
#   ../Output/non_apis_specieslevel_specialization.csv
#   ../Output/non_apis_specialization_models.csv
#   ../Output/non_apis_specialization_vs_honeybee.png

# 0. Packages ####

library(tidyverse)
library(bipartite)
library(broom)
library(stringr)
library(lme4)
library(broom.mixed)

# 0.1 Plot style and colors ####

fau_blue <- "#003366"
fau_red  <- "#CC0000"
fau_gray <- "#CCCCCC"
fau_dark_gray <- "#4D4C55"

hb_low  <- "#E8E8E8"
hb_mid  <- "#E6A3A3"
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

park_col  <- "Park.name"
plant_col <- "Flower_species"
poll_col  <- "Taxon.name"

apis_name <- "Apis mellifera"

# Inspect if needed
names(interactions)
names(relabund)
names(netMets)

# 3. Clean interaction data ####

interactions_clean <- interactions %>%
  filter(
    !is.na(FLW_ID_Conf),
    FLW_ID_Conf != "<80%",
    !is.na(.data[[park_col]]),
    !is.na(.data[[plant_col]]),
    !is.na(.data[[poll_col]])
  ) %>%
  mutate(
    across(
      all_of(c(park_col, plant_col, poll_col)),
      ~ str_squish(as.character(.x))
    )
  ) %>%
  filter(
    .data[[park_col]] != "",
    .data[[plant_col]] != "",
    .data[[poll_col]] != ""
  ) %>%
  mutate(
    park_join = str_trim(str_to_lower(.data[[park_col]])),
    is_apis = .data[[poll_col]] == apis_name
  )

# 4. Function to make plant x pollinator web ####

make_web <- function(df, plant_col, poll_col) {
  
  web_df <- df %>%
    count(.data[[plant_col]], .data[[poll_col]], name = "freq")
  
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

# 5. Calculate species-level specialization for each park ####
# For bipartite matrices:
#   rows = plants = lower level
#   columns = pollinators = higher level
#
# We use specieslevel(..., level = "higher", index = "d")
# because pollinators are columns.
#
# d' is a species-level specialization metric.
# Higher d' means a pollinator is more specialized relative to available resources.

park_list <- split(interactions_clean, interactions_clean[[park_col]])

non_apis_specialization <- imap_dfr(park_list, function(df_park, park_name) {
  
  web <- make_web(
    df = df_park,
    plant_col = plant_col,
    poll_col = poll_col
  )
  
  # Need at least 2 plants and 2 pollinators for species-level d' to behave
  if (nrow(web) < 2 || ncol(web) < 2 || sum(web) == 0) {
    return(tibble())
  }
  
  sp <- tryCatch({
    bipartite::specieslevel(
      web,
      level = "higher",
      index = "d"
    )
  }, error = function(e) {
    NULL
  })
  
  if (is.null(sp)) {
    return(tibble())
  }
  
  # specieslevel can return a matrix/data.frame with pollinator taxa as rownames
  sp_df <- as.data.frame(sp) %>%
    rownames_to_column(var = "pollinator_taxon")
  
  # The d' column name can vary slightly depending on bipartite version.
  # Detect it robustly.
  d_col <- names(sp_df)[names(sp_df) %in% c("d", "d'", "d.prime", "dprime")]
  
  if (length(d_col) == 0) {
    d_col <- names(sp_df)[str_detect(names(sp_df), "^d")]
  }
  
  if (length(d_col) == 0) {
    stop("Could not find species-level d' column in specieslevel() output.")
  }
  
  d_col <- d_col[1]
  
  poll_counts <- df_park %>%
    count(.data[[poll_col]], name = "pollinator_interactions") %>%
    rename(pollinator_taxon = all_of(poll_col))
  
  sp_df %>%
    transmute(
      park_name = park_name,
      park_join = str_trim(str_to_lower(park_name)),
      pollinator_taxon = pollinator_taxon,
      d_prime = as.numeric(.data[[d_col]])
    ) %>%
    left_join(poll_counts, by = "pollinator_taxon") %>%
    filter(
      pollinator_taxon != apis_name,
      !is.na(d_prime)
    )
})

# Inspect
glimpse(non_apis_specialization)

# 6. Prepare park-level predictors ####

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

# 7. Join specialization to park-level predictors ####

spec_dat <- non_apis_specialization %>%
  left_join(relabund_clean, by = "park_join") %>%
  left_join(netMets_clean, by = "park_join") %>%
  mutate(
    log_pollinator_interactions = log1p(pollinator_interactions)
  )

# Inspect joined data
glimpse(spec_dat)

cat("\n--- SPECIALIZATION DATA SUMMARY ---\n")
cat("Rows:", nrow(spec_dat), "\n")
cat("Parks:", n_distinct(spec_dat$park_name), "\n")
cat("Non-Apis pollinator taxa:", n_distinct(spec_dat$pollinator_taxon), "\n")

cat("\nRows missing honey bee relative abundance:\n")
print(spec_dat %>% filter(is.na(hb_raw)) %>% distinct(park_name, park_join))

cat("\nRows missing interaction_sum:\n")
print(spec_dat %>% filter(is.na(interaction_sum)) %>% distinct(park_name, park_join))

# Save species-level specialization table
write.csv(
  spec_dat,
  "../Output/non_apis_specieslevel_specialization.csv",
  row.names = FALSE
)

# 8. Filter to usable model data ####
# Optional but recommended: exclude pollinator taxa observed only once in a park.
# d' for singletons can be very noisy.
#
# Start with >= 2 interactions per pollinator within a park.
# You can change this to >= 1 if sample size gets too small.

spec_model_dat <- spec_dat %>%
  filter(
    !is.na(d_prime),
    !is.na(hb_raw),
    !is.na(log_interaction_sum),
    !is.na(log_pollinator_interactions),
    pollinator_interactions >= 5
  )

cat("\n--- MODEL DATA SUMMARY ---\n")
cat("Rows:", nrow(spec_model_dat), "\n")
cat("Parks:", n_distinct(spec_model_dat$park_name), "\n")
cat("Non-Apis pollinator taxa:", n_distinct(spec_model_dat$pollinator_taxon), "\n")

# 9. Models ####

# Model A:
# Simple park-level predictor model.
# This asks whether non-Apis d' changes with honey bee relative abundance,
# controlling for the total number of interactions observed in that park.

mod_lm_effort <- lm(
  d_prime ~ hb_raw + log_interaction_sum,
  data = spec_model_dat
)

# Model B:
# Adds the number of observations of each pollinator taxon within each park.
# This is important because species-level specialization can be noisy for rare taxa.

mod_lm_effort_taxonobs <- lm(
  d_prime ~ hb_raw + log_interaction_sum + log_pollinator_interactions,
  data = spec_model_dat
)

# Model C:
# Mixed model with pollinator taxon as a random intercept.
# This asks whether the same pollinator taxa tend to be more specialized
# in high-honey-bee parks.
#
# This may be singular if many taxa occur in only one park.

mod_mixed_taxon <- tryCatch({
  lmer(
    d_prime ~ hb_raw + log_interaction_sum + log_pollinator_interactions +
      (1 | pollinator_taxon),
    data = spec_model_dat,
    REML = FALSE
  )
}, error = function(e) {
  message("Mixed model failed: ", e$message)
  NULL
})

# 10. Summarize models ####

lm_effort_results <- tidy(mod_lm_effort) %>%
  mutate(
    model = "lm_effort",
    n = nobs(mod_lm_effort),
    r.squared = glance(mod_lm_effort)$r.squared,
    adj.r.squared = glance(mod_lm_effort)$adj.r.squared,
    AIC = glance(mod_lm_effort)$AIC
  )

lm_effort_taxonobs_results <- tidy(mod_lm_effort_taxonobs) %>%
  mutate(
    model = "lm_effort_taxonobs",
    n = nobs(mod_lm_effort_taxonobs),
    r.squared = glance(mod_lm_effort_taxonobs)$r.squared,
    adj.r.squared = glance(mod_lm_effort_taxonobs)$adj.r.squared,
    AIC = glance(mod_lm_effort_taxonobs)$AIC
  )

if (!is.null(mod_mixed_taxon)) {
  mixed_taxon_results <- broom.mixed::tidy(
    mod_mixed_taxon,
    effects = "fixed"
  ) %>%
    mutate(
      model = "mixed_taxon",
      n = nobs(mod_mixed_taxon),
      r.squared = NA_real_,
      adj.r.squared = NA_real_,
      AIC = AIC(mod_mixed_taxon)
    )
} else {
  mixed_taxon_results <- tibble()
}

model_results <- bind_rows(
  lm_effort_results,
  lm_effort_taxonobs_results,
  mixed_taxon_results
) %>%
  relocate(model, n)

cat("\n--- NON-APIS SPECIALIZATION MODEL RESULTS ---\n")
print(model_results)

write.csv(
  model_results,
  "../Output/non_apis_specialization_models.csv",
  row.names = FALSE
)

# 11. Build adjusted prediction for main plot ####
# Use Model B for plotting because it controls for:
#   1. total park interaction effort
#   2. pollinator taxon observation count
#
# Predictions hold both controls at their median values.

pred_dat <- tibble(
  hb_raw = seq(
    min(spec_model_dat$hb_raw, na.rm = TRUE),
    max(spec_model_dat$hb_raw, na.rm = TRUE),
    length.out = 100
  ),
  log_interaction_sum = median(spec_model_dat$log_interaction_sum, na.rm = TRUE),
  log_pollinator_interactions = median(spec_model_dat$log_pollinator_interactions, na.rm = TRUE)
)

preds <- predict(
  mod_lm_effort_taxonobs,
  newdata = pred_dat,
  interval = "confidence"
) %>%
  as.data.frame()

pred_dat <- bind_cols(pred_dat, preds) %>%
  mutate(
    hb_percent = hb_raw * 100
  )

main_p <- tidy(mod_lm_effort_taxonobs) %>%
  filter(term == "hb_raw") %>%
  pull(p.value)

main_r2 <- glance(mod_lm_effort_taxonobs)$r.squared

annotation_text <- paste0(
  "Adjusted model",
  "\nR\u00B2 = ", round(main_r2, 2),
  "\nP_honey bee = ", signif(main_p, 2),
  "\nN = ", nrow(spec_model_dat), " taxon-park pairs"
)

# 12. Plot non-Apis specialization vs honey bee relative abundance ####

non_apis_specialization_plot <- ggplot(
  spec_model_dat,
  aes(
    x = hb_raw,
    y = d_prime
  )
) +
  geom_ribbon(
    data = pred_dat,
    aes(
      x = hb_raw,
      ymin = lwr,
      ymax = upr
    ),
    inherit.aes = FALSE,
    fill = fau_gray,
    alpha = 0.45
  ) +
  geom_line(
    data = pred_dat,
    aes(
      x = hb_raw,
      y = fit
    ),
    inherit.aes = FALSE,
    color = fau_red,
    linewidth = 1.2
  ) +
  geom_point(
    aes(
      fill = hb_percent,
      size = pollinator_interactions
    ),
    shape = 21,
    color = fau_dark_gray,
    stroke = 0.7,
    alpha = 0.85
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = annotation_text,
    hjust = 1.08,
    vjust = 1.25,
    size = 4.5,
    color = "grey20"
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.03, 0.07))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  scale_fill_gradientn(
    colors = hb_ramp,
    limits = c(0, 100),
    name = "% honey bee"
  ) +
  scale_size_continuous(
    range = c(2.5, 8),
    name = "Pollinator\nobservations"
  ) +
  labs(
    x = expression(italic("Apis mellifera") ~ "relative abundance"),
    y = "Non-Apis pollinator specialization (d')",
    title = "Non-Apis pollinator specialization vs. honey bee relative abundance",
    subtitle = "Line controls for total park interactions and pollinator observation count"
  ) +
  theme_hb(base_size = 17) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 22),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

non_apis_specialization_plot

ggsave(
  "../Output/non_apis_specialization_vs_honeybee.png",
  non_apis_specialization_plot,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

# 13. Optional park-level summary plot ####
# Average non-Apis specialization per park.
# This is less statistically detailed but easier to explain on a poster.

park_spec_summary <- spec_model_dat %>%
  group_by(park_name, park_join) %>%
  summarise(
    mean_d_prime = mean(d_prime, na.rm = TRUE),
    median_d_prime = median(d_prime, na.rm = TRUE),
    n_non_apis_taxa = n_distinct(pollinator_taxon),
    n_taxon_park_pairs = n(),
    hb_raw = first(hb_raw),
    hb_percent = first(hb_percent),
    interaction_sum = first(interaction_sum),
    log_interaction_sum = first(log_interaction_sum),
    .groups = "drop"
  )

park_mod <- lm(
  mean_d_prime ~ hb_raw + log_interaction_sum,
  data = park_spec_summary
)

park_p <- tidy(park_mod) %>%
  filter(term == "hb_raw") %>%
  pull(p.value)

park_r2 <- glance(park_mod)$r.squared

park_annotation_text <- paste0(
  "Park-level model",
  "\nR\u00B2 = ", round(park_r2, 2),
  "\nP_honey bee = ", signif(park_p, 2),
  "\nN = ", nrow(park_spec_summary), " parks"
)

park_pred <- tibble(
  hb_raw = seq(
    min(park_spec_summary$hb_raw, na.rm = TRUE),
    max(park_spec_summary$hb_raw, na.rm = TRUE),
    length.out = 100
  ),
  log_interaction_sum = median(park_spec_summary$log_interaction_sum, na.rm = TRUE)
)

park_preds <- predict(
  park_mod,
  newdata = park_pred,
  interval = "confidence"
) %>%
  as.data.frame()

park_pred <- bind_cols(park_pred, park_preds)

park_specialization_plot <- ggplot(
  park_spec_summary,
  aes(
    x = hb_raw,
    y = mean_d_prime
  )
) +
  geom_ribbon(
    data = park_pred,
    aes(
      x = hb_raw,
      ymin = lwr,
      ymax = upr
    ),
    inherit.aes = FALSE,
    fill = fau_gray,
    alpha = 0.45
  ) +
  geom_line(
    data = park_pred,
    aes(
      x = hb_raw,
      y = fit
    ),
    inherit.aes = FALSE,
    color = fau_red,
    linewidth = 1.2
  ) +
  geom_point(
    aes(
      fill = hb_percent,
      size = n_non_apis_taxa
    ),
    shape = 21,
    color = fau_dark_gray,
    stroke = 0.8,
    alpha = 0.9
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = park_annotation_text,
    hjust = 1.08,
    vjust = 1.25,
    size = 4.5,
    color = "grey20"
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.03, 0.07))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  scale_fill_gradientn(
    colors = hb_ramp,
    limits = c(0, 100),
    name = "% honey bee"
  ) +
  scale_size_continuous(
    range = c(3, 8),
    name = "Non-Apis\ntaxa"
  ) +
  labs(
    x = expression(italic("Apis mellifera") ~ "relative abundance"),
    y = "Mean non-Apis specialization (d')",
    title = "Mean non-Apis specialization by park",
    subtitle = "Line controls for total park interactions observed"
  ) +
  theme_hb(base_size = 17) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 22),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

park_specialization_plot

ggsave(
  "../Output/non_apis_mean_specialization_by_park_vs_honeybee.png",
  park_specialization_plot,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

write.csv(
  park_spec_summary,
  "../Output/non_apis_park_specialization_summary.csv",
  row.names = FALSE
)

