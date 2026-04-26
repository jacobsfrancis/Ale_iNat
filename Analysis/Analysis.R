# Purpose: Test whether relative abundance of honey bees in a park
# predicts park-level plant-pollinator network structure while controlling
# for total number of plant-pollinator interactions observed.
#
# Inputs:
#   ../Output/ApisRelAbund.csv
#   ../Output/park_level_network_metrics.csv
#
# Outputs:
#   Model summaries, adjusted model summaries, correlation matrix,
#   and styled regression plots.

# Libraries ####
library(tidyverse)
library(broom)
library(stringr)

# 0. Plot style and colors ####

# FAU approved / inspired colors
fau_blue <- "#003366"
fau_red  <- "#CC0000"
fau_gray <- "#CCCCCC"
fau_dark_gray <- "#4D4C55"

# Honey bee alarm color ramp:
# low = pale gray, mid = muted pink/red, high = FAU red
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

clean_metric_label <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace("H2", "H\u00B2") %>%
    str_to_sentence()
}

# 1. Read in data ####

relabund <- read.csv("../Output/ApisRelAbund.csv")
netMets  <- read.csv("../Output/park_level_network_metrics.csv", stringsAsFactors = FALSE)

# Inspect ####
names(relabund)
names(netMets)
glimpse(relabund)
glimpse(netMets)

# 2. Define join column ####

park_col_relabund <- "poly_id"
park_col_netmets  <- "park_name"

# 3. Clean park names before join ####

relabund <- relabund %>%
  mutate(
    park_join = str_trim(str_to_lower(.data[[park_col_relabund]]))
  )

netMets <- netMets %>%
  mutate(
    park_join = str_trim(str_to_lower(.data[[park_col_netmets]]))
  )

# 4. Check overlap between datasets ####

parks_relabund <- relabund %>% distinct(park_join)
parks_netmets  <- netMets %>% distinct(park_join)

parks_in_both <- inner_join(parks_relabund, parks_netmets, by = "park_join")
parks_only_relabund <- anti_join(parks_relabund, parks_netmets, by = "park_join")
parks_only_netmets  <- anti_join(parks_netmets, parks_relabund, by = "park_join")

cat("\n--- PARK OVERLAP SUMMARY ---\n")
cat("Parks in relabund:", nrow(parks_relabund), "\n")
cat("Parks in netMets :", nrow(parks_netmets), "\n")
cat("Parks in both    :", nrow(parks_in_both), "\n")
cat("Only in relabund :", nrow(parks_only_relabund), "\n")
cat("Only in netMets  :", nrow(parks_only_netmets), "\n")

if (nrow(parks_only_relabund) > 0) {
  cat("\nParks only in relabund:\n")
  print(parks_only_relabund)
}

if (nrow(parks_only_netmets) > 0) {
  cat("\nParks only in netMets:\n")
  print(parks_only_netmets)
}

# 5. Join datasets ####

dat <- relabund %>%
  inner_join(netMets, by = "park_join", suffix = c("_relabund", "_netmets"))

cat("\nRows in joined dataset:", nrow(dat), "\n")

# 6. Define honey bee relative abundance column ####
# IMPORTANT:
# Use percentageHB, not n_hb.
# n_hb is the number of honey bee observations.
# percentageHB is the relative abundance metric.

hb_col <- "percentageHB"

# 7. Create honey bee relative abundance variables ####

dat <- dat %>%
  mutate(
    hb_percent = as.numeric(.data[[hb_col]]),
    hb_raw = hb_percent / 100
  )

# Bound to (0,1) for logit transform
eps <- 0.001

dat <- dat %>%
  mutate(
    hb_bounded = pmin(pmax(hb_raw, eps), 1 - eps),
    hb_logit = qlogis(hb_bounded),
    hb_arcsine = asin(sqrt(hb_raw))
  )

summary(dat$hb_raw)
summary(dat$hb_logit)

# 7.1 Define sampling effort control ####
# interaction_sum comes from the network metric file and represents the
# total number of plant-pollinator interactions used to build each park network.

if (!"interaction_sum" %in% names(dat)) {
  stop("interaction_sum not found in dat. Check names(netMets) for the interaction-count column.")
}

dat <- dat %>%
  mutate(
    interaction_sum = as.numeric(interaction_sum),
    log_interaction_sum = log1p(interaction_sum)
  )

summary(dat$interaction_sum)
summary(dat$log_interaction_sum)

# 8. Choose network metrics to model ####

candidate_metrics <- c(
  "connectance",
  "web_asymmetry",
  "nestedness",
  "modularity",
  "H2",
  "linkage_density",
  "interaction_evenness"
)

response_vars <- candidate_metrics[candidate_metrics %in% names(dat)]

cat("\nResponse variables being modeled:\n")
print(response_vars)

# 9. Unadjusted models: metric ~ honey bee relative abundance ####
# These are useful as a baseline but do NOT control for sampling effort.

fit_one_model <- function(response, predictor = "hb_raw", data = dat) {
  
  model_df <- data %>%
    select(all_of(c(response, predictor))) %>%
    drop_na()
  
  if (nrow(model_df) < 5) {
    return(tibble(
      response = response,
      predictor = predictor,
      n = nrow(model_df),
      term = predictor,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      r.squared = NA_real_,
      adj.r.squared = NA_real_,
      AIC = NA_real_
    ))
  }
  
  if (sd(model_df[[response]], na.rm = TRUE) == 0) {
    return(tibble(
      response = response,
      predictor = predictor,
      n = nrow(model_df),
      term = predictor,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      r.squared = NA_real_,
      adj.r.squared = NA_real_,
      AIC = NA_real_
    ))
  }
  
  form <- as.formula(paste(response, "~", predictor))
  mod <- lm(form, data = model_df)
  
  coef_tab <- tidy(mod) %>%
    filter(term == predictor)
  
  g <- glance(mod)
  
  tibble(
    response = response,
    predictor = predictor,
    n = nrow(model_df),
    term = predictor,
    estimate = coef_tab$estimate,
    std.error = coef_tab$std.error,
    statistic = coef_tab$statistic,
    p.value = coef_tab$p.value,
    r.squared = g$r.squared,
    adj.r.squared = g$adj.r.squared,
    AIC = g$AIC
  )
}

model_results_raw <- map_dfr(response_vars, fit_one_model, predictor = "hb_raw", data = dat)
model_results_logit <- map_dfr(response_vars, fit_one_model, predictor = "hb_logit", data = dat)

model_results_raw <- model_results_raw %>%
  mutate(p_adj_fdr = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.value)

model_results_logit <- model_results_logit %>%
  mutate(p_adj_fdr = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.value)

cat("\n--- RESULTS: UNADJUSTED RAW HONEY BEE RELATIVE ABUNDANCE ---\n")
print(model_results_raw)

cat("\n--- RESULTS: UNADJUSTED LOGIT HONEY BEE RELATIVE ABUNDANCE ---\n")
print(model_results_logit)

write.csv(
  model_results_raw,
  "../Output/honeybee_network_metric_models_raw.csv",
  row.names = FALSE
)

write.csv(
  model_results_logit,
  "../Output/honeybee_network_metric_models_logit.csv",
  row.names = FALSE
)

# 10. Adjusted models: metric ~ honey bee relative abundance + interaction effort ####
# Main model:
#   metric ~ honey bee relative abundance + log(total plant-pollinator interactions observed)

fit_one_adjusted_model <- function(
    response,
    predictor = "hb_raw",
    effort = "log_interaction_sum",
    data = dat
) {
  
  model_df <- data %>%
    select(all_of(c(response, predictor, effort))) %>%
    drop_na()
  
  if (nrow(model_df) < 5) {
    return(tibble(
      response = response,
      predictor = predictor,
      effort = effort,
      n = nrow(model_df),
      term = predictor,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      r.squared = NA_real_,
      adj.r.squared = NA_real_,
      AIC = NA_real_
    ))
  }
  
  if (sd(model_df[[response]], na.rm = TRUE) == 0) {
    return(tibble(
      response = response,
      predictor = predictor,
      effort = effort,
      n = nrow(model_df),
      term = predictor,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      r.squared = NA_real_,
      adj.r.squared = NA_real_,
      AIC = NA_real_
    ))
  }
  
  form <- as.formula(paste(response, "~", predictor, "+", effort))
  mod <- lm(form, data = model_df)
  
  coef_tab <- tidy(mod) %>%
    filter(term == predictor)
  
  g <- glance(mod)
  
  tibble(
    response = response,
    predictor = predictor,
    effort = effort,
    n = nrow(model_df),
    term = predictor,
    estimate = coef_tab$estimate,
    std.error = coef_tab$std.error,
    statistic = coef_tab$statistic,
    p.value = coef_tab$p.value,
    r.squared = g$r.squared,
    adj.r.squared = g$adj.r.squared,
    AIC = g$AIC
  )
}

model_results_adjusted_raw <- map_dfr(
  response_vars,
  fit_one_adjusted_model,
  predictor = "hb_raw",
  effort = "log_interaction_sum",
  data = dat
)

model_results_adjusted_logit <- map_dfr(
  response_vars,
  fit_one_adjusted_model,
  predictor = "hb_logit",
  effort = "log_interaction_sum",
  data = dat
)

model_results_adjusted_raw <- model_results_adjusted_raw %>%
  mutate(p_adj_fdr = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.value)

model_results_adjusted_logit <- model_results_adjusted_logit %>%
  mutate(p_adj_fdr = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.value)

cat("\n--- RESULTS: ADJUSTED RAW HONEY BEE RELATIVE ABUNDANCE ---\n")
print(model_results_adjusted_raw)

cat("\n--- RESULTS: ADJUSTED LOGIT HONEY BEE RELATIVE ABUNDANCE ---\n")
print(model_results_adjusted_logit)

# 10.1 Save full adjusted model coefficient tables ####
# This saves all terms, not just the honey bee term.

full_adjusted_model_coefficients <- map_dfr(response_vars, function(resp) {
  
  model_df <- dat %>%
    select(all_of(c(resp, "hb_raw", "log_interaction_sum"))) %>%
    drop_na()
  
  if (nrow(model_df) < 5 || sd(model_df[[resp]], na.rm = TRUE) == 0) {
    return(tibble())
  }
  
  mod <- lm(
    as.formula(paste(resp, "~ hb_raw + log_interaction_sum")),
    data = model_df
  )
  
  tidy(mod) %>%
    mutate(
      response = resp,
      n = nrow(model_df),
      r.squared = glance(mod)$r.squared,
      adj.r.squared = glance(mod)$adj.r.squared,
      AIC = glance(mod)$AIC
    ) %>%
    relocate(response, n)
})

cat("\n--- FULL ADJUSTED MODEL COEFFICIENTS ---\n")
print(full_adjusted_model_coefficients)

write.csv(
  model_results_adjusted_raw,
  "../Output/honeybee_network_metric_models_adjusted_raw.csv",
  row.names = FALSE
)

write.csv(
  model_results_adjusted_logit,
  "../Output/honeybee_network_metric_models_adjusted_logit.csv",
  row.names = FALSE
)

write.csv(
  full_adjusted_model_coefficients,
  "../Output/honeybee_network_metric_models_adjusted_full_coefficients.csv",
  row.names = FALSE
)

# 11. Quick unadjusted styled plots for each response ####

plot_dir_unadjusted <- "../Output/network_metric_regression_plots"
if (!dir.exists(plot_dir_unadjusted)) dir.create(plot_dir_unadjusted, recursive = TRUE)

for (resp in response_vars) {
  
  pdat <- dat %>%
    select(all_of(c(resp, "hb_raw", "hb_percent", "n_apidae", "n_hb"))) %>%
    drop_na()
  
  if (nrow(pdat) < 5) next
  if (sd(pdat[[resp]], na.rm = TRUE) == 0) next
  
  resp_label <- clean_metric_label(resp)
  
  mod <- lm(as.formula(paste(resp, "~ hb_raw")), data = pdat)
  mod_glance <- glance(mod)
  mod_tidy <- tidy(mod)
  
  p_val <- mod_tidy %>%
    filter(term == "hb_raw") %>%
    pull(p.value)
  
  r2_val <- mod_glance$r.squared
  
  annotation_text <- paste0(
    "Unadjusted model",
    "\nR\u00B2 = ", round(r2_val, 2),
    "\nP = ", signif(p_val, 2),
    "\nN = ", nrow(pdat), " parks"
  )
  
  p <- ggplot(
    pdat,
    aes(
      x = hb_raw,
      y = .data[[resp]]
    )
  ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = fau_red,
      fill = fau_gray,
      linewidth = 1.1,
      alpha = 0.45
    ) +
    geom_point(
      aes(
        fill = hb_percent,
        size = n_apidae
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
      label = annotation_text,
      hjust = 1.08,
      vjust = 1.25,
      size = 4.2,
      color = "grey20"
    ) +
    scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0.03, 0.07))
    ) +
    scale_fill_gradientn(
      colors = hb_ramp,
      limits = c(0, 100),
      name = "% honey bee"
    ) +
    scale_size_continuous(
      range = c(2.5, 7),
      name = "Number of\nAnthophila\nobservations"
    ) +
    labs(
      x = expression(italic("Apis mellifera") ~ "relative abundance"),
      y = resp_label,
      title = paste(resp_label, "vs. honey bee relative abundance"),
      subtitle = "Each point represents one park"
    ) +
    theme_hb(base_size = 15) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  print(p)
  
  ggsave(
    filename = file.path(plot_dir_unadjusted, paste0(resp, "_vs_honeybee_unadjusted.png")),
    plot = p,
    width = 7,
    height = 5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

# 12. Adjusted plots for each response ####
# These are not simple geom_smooth plots.
# They show the predicted relationship between honey bee relative abundance
# and the network metric while holding sampling effort at its median value.

plot_dir_adjusted <- "../Output/network_metric_regression_plots_adjusted"
if (!dir.exists(plot_dir_adjusted)) dir.create(plot_dir_adjusted, recursive = TRUE)

for (resp in response_vars) {
  
  pdat <- dat %>%
    select(all_of(c(resp, "hb_raw", "hb_percent", "interaction_sum", "log_interaction_sum"))) %>%
    drop_na()
  
  if (nrow(pdat) < 5) next
  if (sd(pdat[[resp]], na.rm = TRUE) == 0) next
  
  resp_label <- clean_metric_label(resp)
  
  mod <- lm(
    as.formula(paste(resp, "~ hb_raw + log_interaction_sum")),
    data = pdat
  )
  
  mod_glance <- glance(mod)
  mod_tidy <- tidy(mod)
  
  p_val <- mod_tidy %>%
    filter(term == "hb_raw") %>%
    pull(p.value)
  
  r2_val <- mod_glance$r.squared
  
  annotation_text <- paste0(
    "Adjusted model",
    "\nR\u00B2 = ", round(r2_val, 2),
    "\nP_honey bee = ", signif(p_val, 2),
    "\nN = ", nrow(pdat), " parks"
  )
  
  pred_dat <- tibble(
    hb_raw = seq(
      min(pdat$hb_raw, na.rm = TRUE),
      max(pdat$hb_raw, na.rm = TRUE),
      length.out = 100
    ),
    log_interaction_sum = median(pdat$log_interaction_sum, na.rm = TRUE)
  )
  
  preds <- predict(
    mod,
    newdata = pred_dat,
    interval = "confidence"
  ) %>%
    as.data.frame()
  
  pred_dat <- bind_cols(pred_dat, preds)
  
  p <- ggplot(
    pdat,
    aes(
      x = hb_raw,
      y = .data[[resp]]
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
        size = interaction_sum
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
      label = annotation_text,
      hjust = 1.08,
      vjust = 1.25,
      size = 4.2,
      color = "grey20"
    ) +
    scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0.03, 0.07))
    ) +
    scale_fill_gradientn(
      colors = hb_ramp,
      limits = c(0, 100),
      name = "% honey bee"
    ) +
    scale_size_continuous(
      range = c(2.5, 7),
      name = "Plant-pollinator\ninteractions"
    ) +
    labs(
      x = expression(italic("Apis mellifera") ~ "relative abundance"),
      y = resp_label,
      title = paste(resp_label, "vs. honey bee relative abundance"),
      subtitle = "Line controls for total plant-pollinator interactions observed"
    ) +
    theme_hb(base_size = 15) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  print(p)
  
  ggsave(
    filename = file.path(plot_dir_adjusted, paste0(resp, "_vs_honeybee_adjusted.png")),
    plot = p,
    width = 7,
    height = 5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

# 13. Main poster-style adjusted nestedness plot ####

if ("nestedness" %in% names(dat)) {
  
  nested_dat <- dat %>%
    select(nestedness, hb_raw, hb_percent, interaction_sum, log_interaction_sum) %>%
    drop_na()
  
  nested_mod <- lm(
    nestedness ~ hb_raw + log_interaction_sum,
    data = nested_dat
  )
  
  nested_glance <- glance(nested_mod)
  nested_tidy <- tidy(nested_mod)
  
  nested_p <- nested_tidy %>%
    filter(term == "hb_raw") %>%
    pull(p.value)
  
  nested_r2 <- nested_glance$r.squared
  
  annotation_text <- paste0(
    "Adjusted model",
    "\nR\u00B2 = ", round(nested_r2, 2),
    "\nP_honey bee = ", signif(nested_p, 2),
    "\nN = ", nrow(nested_dat), " parks"
  )
  
  nested_pred <- tibble(
    hb_raw = seq(
      min(nested_dat$hb_raw, na.rm = TRUE),
      max(nested_dat$hb_raw, na.rm = TRUE),
      length.out = 100
    ),
    log_interaction_sum = median(nested_dat$log_interaction_sum, na.rm = TRUE)
  )
  
  nested_preds <- predict(
    nested_mod,
    newdata = nested_pred,
    interval = "confidence"
  ) %>%
    as.data.frame()
  
  nested_pred <- bind_cols(nested_pred, nested_preds)
  
  nestedness_plot <- ggplot(
    nested_dat,
    aes(
      x = hb_raw,
      y = nestedness
    )
  ) +
    geom_ribbon(
      data = nested_pred,
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
      data = nested_pred,
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
        size = interaction_sum
      ),
      shape = 21,
      color = fau_dark_gray,
      stroke = 0.9,
      alpha = 0.9
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = annotation_text,
      hjust = 1.08,
      vjust = 1.25,
      size = 4.8,
      color = "grey20"
    ) +
    scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0.03, 0.07))
    ) +
    scale_fill_gradientn(
      colors = hb_ramp,
      limits = c(0, 100),
      name = "% honey bee"
    ) +
    scale_size_continuous(
      range = c(3, 8),
      name = "Plant-pollinator\ninteractions"
    ) +
    labs(
      x = expression(italic("Apis mellifera") ~ "relative abundance"),
      y = "Nestedness",
      title = "Nestedness vs. honey bee relative abundance",
      subtitle = "Regression controls for total plant-pollinator interactions observed"
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
  
  nestedness_plot
  
  ggsave(
    "../Output/nestedness_vs_honeybee_relative_abundance_adjusted_poster.png",
    nestedness_plot,
    width = 8,
    height = 5.5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

# 14. Correlation matrix ####
# Use numeric columns only; include honey bee abundance, effort, and network metrics.

corr_vars <- c("hb_raw", "hb_percent", "interaction_sum", "log_interaction_sum", response_vars)
corr_vars <- corr_vars[corr_vars %in% names(dat)]

corr_df <- dat %>%
  select(all_of(corr_vars)) %>%
  select(where(is.numeric))

corr_mat <- cor(corr_df, use = "pairwise.complete.obs", method = "pearson")

print(corr_mat)

write.csv(
  corr_mat,
  "../Output/honeybee_network_metric_correlation_matrix.csv"
)

# 15. Draw correlation matrix as a heatmap ####

corr_long <- as.data.frame(as.table(corr_mat)) %>%
  rename(var1 = Var1, var2 = Var2, correlation = Freq) %>%
  mutate(
    var1 = clean_metric_label(as.character(var1)),
    var2 = clean_metric_label(as.character(var2)),
    var1 = str_replace(var1, "Hb raw", "Honey bee relative abundance"),
    var2 = str_replace(var2, "Hb raw", "Honey bee relative abundance"),
    var1 = str_replace(var1, "Hb percent", "% honey bee"),
    var2 = str_replace(var2, "Hb percent", "% honey bee"),
    var1 = str_replace(var1, "Interaction sum", "Plant-pollinator interactions"),
    var2 = str_replace(var2, "Interaction sum", "Plant-pollinator interactions"),
    var1 = str_replace(var1, "Log interaction sum", "Log interactions"),
    var2 = str_replace(var2, "Log interaction sum", "Log interactions")
  )

corr_plot <- ggplot(corr_long, aes(x = var1, y = var2, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = round(correlation, 2)),
    size = 4,
    color = "grey15"
  ) +
  scale_fill_gradient2(
    low = fau_blue,
    mid = "white",
    high = fau_red,
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson\ncorrelation"
  ) +
  coord_equal() +
  theme_hb(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(size = 18)
  ) +
  labs(
    title = "Correlation matrix: honey bee abundance, effort, and network metrics",
    x = NULL,
    y = NULL
  )

corr_plot

ggsave(
  "../Output/honeybee_network_metric_correlation_matrix.png",
  corr_plot,
  width = 8,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

# 16. Optional: simple model comparison table ####
# This puts unadjusted and adjusted honey bee effects side-by-side.

model_comparison <- model_results_raw %>%
  select(
    response,
    n_unadjusted = n,
    estimate_unadjusted = estimate,
    p_unadjusted = p.value,
    r2_unadjusted = r.squared,
    p_adj_fdr_unadjusted = p_adj_fdr
  ) %>%
  left_join(
    model_results_adjusted_raw %>%
      select(
        response,
        n_adjusted = n,
        estimate_adjusted = estimate,
        p_adjusted = p.value,
        r2_adjusted = r.squared,
        adj_r2_adjusted = adj.r.squared,
        p_adj_fdr_adjusted = p_adj_fdr
      ),
    by = "response"
  ) %>%
  arrange(p_adjusted)

cat("\n--- MODEL COMPARISON: UNADJUSTED VS ADJUSTED ---\n")
print(model_comparison)

write.csv(
  model_comparison,
  "../Output/honeybee_network_metric_model_comparison.csv",
  row.names = FALSE
)