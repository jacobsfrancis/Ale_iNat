# File for exploring iNat data that Ale downloaded
# Authors: JSF + AA
# 8 April 2026
#
# Purpose:
#   1. Read Broward County Anthophila iNaturalist observations
#   2. Spatially join bee observations to Broward urban park polygons
#   3. Calculate the percentage of Apidae observations that are honey bees
#      within each park
#   4. Plot park-level honey bee relative abundance using FAU-inspired colors

# Clear workspace if you want ####
rm(list = ls())

# 0.1 Packages ####

library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(tigris)
library(scales)
library(stringr)

options(tigris_use_cache = TRUE)

# 0.2 FAU approved / inspired colors ####

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

# 0.3 Read in data ####

# Apid observations for Broward County
bees <- read.delim("../Data/Broward anthophila obs unfiltered.csv")

# Shapefile for Broward parks
parks <- st_read("../Data/urban_parks/urban_parks.shp")

# Inspect column names if needed
names(bees)
names(parks)

# 1.0 Geospatial join ####

# Make bee points in sf format
bees_sf <- bees %>%
  filter(
    !is.na(decimalLongitude),
    !is.na(decimalLatitude)
  ) %>%
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326
  )

# Transform bees to same CRS as parks
bees_sf <- st_transform(bees_sf, st_crs(parks))

# Spatial join: attach park attributes to each bee point
bees_by_park <- st_join(bees_sf, parks, join = st_within)

# Look at joined data
head(bees_by_park)

# 2.0 Remove observations not in parks ####

bees_in_parks <- bees_by_park %>%
  filter(!is.na(poly_id))

# 3.0 Calculate percentage of Apidae observations that are honey bees ####

Apid_Summary <- bees_in_parks %>%
  group_by(poly_id) %>%
  summarise(
    n_apidae = n(),
    n_hb = sum(genus == "Apis", na.rm = TRUE),
    percentageHB = 100 * n_hb / n_apidae,
    .groups = "drop"
  )

# Keep only parks with more than 5 Apidae observations
Apid_Summary_20 <- Apid_Summary %>%
  filter(n_apidae > 5)

# View summary table
as.data.frame(Apid_Summary_20)

# Write summary output
write.csv(
  as.data.frame(Apid_Summary_20)[, 1:4],
  "../Output/ApisRelAbund.csv",
  row.names = FALSE
)

# 4.0 Exploratory ranked dot plot ####

# Clean park names for plotting
Apid_Summary_20_plot <- Apid_Summary_20 %>%
  mutate(
    park_label = as.character(poly_id),
    park_label = str_replace_all(park_label, "_", " "),
    park_label = str_replace(park_label, "^Broward County ", ""),
    park_label = str_replace(park_label, "^Broward ", ""),
    park_label = str_replace(park_label, " Preserve Park$", " Preserve"),
    park_label = str_replace(park_label, " Nature Area$", " Natural Area"),
    park_label = str_squish(park_label)
  )

# 4.0 Exploratory ranked dot plot ####

# Clean park names for plotting
Apid_Summary_20_plot <- Apid_Summary_20 %>%
  mutate(
    park_label = as.character(poly_id),
    park_label = str_replace_all(park_label, "_", " "),
    park_label = str_replace(park_label, "^Broward County ", ""),
    park_label = str_replace(park_label, "^Broward ", ""),
    park_label = str_replace(park_label, " Preserve Park$", " Preserve"),
    park_label = str_replace(park_label, " Nature Area$", " Natural Area"),
    park_label = str_squish(park_label)
  )

explore_Plot1 <- ggplot(
  data = Apid_Summary_20_plot,
  aes(
    x = reorder(park_label, percentageHB),
    y = percentageHB,
    color = percentageHB
  )
) +
  geom_point(
    aes(size = n_apidae),
    alpha = 0.9
  ) +
  geom_text(
    aes(
      y = 106,
      label = n_apidae
    ),
    color = "grey25",
    size = 4.5,
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = length(unique(Apid_Summary_20_plot$park_label)) + 0.8,
    y = 106,
    label = "N bee obs.",
    hjust = 0.5,
    size = 4.8,
    fontface = "bold",
    color = "grey20"
  ) +
  coord_flip(clip = "off") +
  scale_color_gradientn(
    colors = hb_ramp,
    name = "% honey bee"
  ) +
  scale_y_continuous(
    limits = c(0, 110),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_size_continuous(
    range = c(3, 9),
    guide = "none"
  ) +
  theme_minimal(base_size = 17) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 18),
    legend.position = "none",
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 15),
    plot.title = element_text(face = "bold", size = 22),
    plot.subtitle = element_text(size = 17),
    plot.margin = margin(10, 35, 10, 10)
  ) +
  labs(
    title = expression(italic("Apis") ~ " relative abundance by park"),
    subtitle = "Parks with more than 5 Anthophila observations",
    x = NULL,
    y = "Honey bee observations"
  )

explore_Plot1

ggsave(
  "../Output/apis_relative_abundance_dotplot_poster.png",
  plot = explore_Plot1,
  width = 14,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)


# 5.0 Map honey bee relative abundance by park ####

# Only plot in urban portion of Broward
broward_urban <- st_read(
  "../Data/urban_parks/Broward_urban_border.shp",
  quiet = TRUE
)

broward_urban <- st_transform(broward_urban, st_crs(parks))

# Make sure poly_id joins cleanly.
# Using as.character() avoids problems if one object treats poly_id as numeric
# and the other treats it as character.
parks_clean <- parks %>%
  mutate(
    poly_id = as.character(poly_id)
  )

apis_map_data <- Apid_Summary_20 %>%
  st_drop_geometry() %>%
  mutate(
    poly_id = as.character(poly_id)
  ) %>%
  select(
    poly_id,
    n_apidae,
    n_hb,
    percentageHB
  )

# Join honey bee relative abundance back onto the park polygons
parks_apis <- parks_clean %>%
  left_join(apis_map_data, by = "poly_id")

# Make sure everything is in the same CRS
parks_apis <- st_transform(parks_apis, st_crs(broward_urban))

# Keep only park polygons that intersect the urban Broward boundary
parks_apis_urban <- parks_apis %>%
  st_filter(broward_urban, .predicate = st_intersects)

# Optional: clip park polygons to the urban Broward boundary.
# This is visually cleaner if any park polygon crosses the boundary.
parks_apis_urban <- st_intersection(
  st_make_valid(parks_apis_urban),
  st_make_valid(broward_urban)
)

# Diagnostic: summary rows that did not join to park polygons
missing_from_parks <- anti_join(
  apis_map_data %>% distinct(poly_id),
  parks_clean %>% st_drop_geometry() %>% distinct(poly_id),
  by = "poly_id"
)

missing_from_parks

# Diagnostic: urban park polygons with no honey bee summary data
parks_without_data <- parks_apis_urban %>%
  st_drop_geometry() %>%
  filter(is.na(percentageHB)) %>%
  select(poly_id)

parks_without_data

# Main park-level honey bee relative abundance map
apis_map <- ggplot() +
  geom_sf(
    data = broward_urban,
    fill = "grey95",
    color = fau_dark_gray,
    linewidth = 0.45
  ) +
  geom_sf(
    data = parks_apis_urban,
    aes(fill = percentageHB),
    color = NA,
    linewidth = 0
  ) +
  scale_fill_gradientn(
    colors = hb_ramp,
    labels = function(x) paste0(x, "%"),
    na.value = "grey85",
    name = "% honey bee\nobservations"
  ) +
  coord_sf(
    xlim = st_bbox(broward_urban)[c("xmin", "xmax")],
    ylim = st_bbox(broward_urban)[c("ymin", "ymax")],
    expand = FALSE
  ) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

apis_map

ggsave(
  "../Output/broward_urban_apis_relative_abundance_map.png",
  apis_map,
  width = 5.5,
  height = 5.5,
  dpi = 300,
  bg = "transparent"
)


####extra maps not used ####
# 5.0 Quick map of all bee observations##

bee_points_map <- ggplot() +
  geom_sf(
    data = parks,
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_sf(
    data = bees_sf,
    alpha = 0.5,
    size = 0.8
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    title = "Broward County bee observations",
    subtitle = "All Anthophila observations before filtering to park interiors"
  )

bee_points_map

ggsave(
  "../Output/broward_bee_observations_map.png",
  bee_points_map,
  width = 8,
  height = 6,
  dpi = 300
)

# 6.0 Optional density map of bee observations #

bees_coords <- bees_sf %>%
  st_transform(3857) %>%
  cbind(st_coordinates(.))

bee_density_map <- ggplot() +
  geom_sf(
    data = st_transform(parks, 3857),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  stat_density_2d(
    data = bees_coords,
    aes(
      X,
      Y,
      fill = after_stat(level),
      alpha = after_stat(level)
    ),
    geom = "polygon",
    contour = TRUE
  ) +
  geom_point(
    data = bees_coords,
    aes(X, Y),
    size = 0.4,
    alpha = 0.4
  ) +
  scale_fill_gradientn(
    colors = hb_ramp,
    guide = "none"
  ) +
  scale_alpha(
    range = c(0.05, 0.35),
    guide = "none"
  ) +
  coord_sf(crs = st_crs(3857)) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    title = "Spatial density of Broward County bee observations"
  )

bee_density_map

ggsave(
  "../Output/broward_bee_observation_density_map.png",
  bee_density_map,
  width = 8,
  height = 6,
  dpi = 300
)

