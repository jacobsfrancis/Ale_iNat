# File for exploring iNat data that Ale Downloaded
# Authors: JSF + AA
# 8 April 2026

#Clear Workspace if you want

rm(list=ls())

#0.1 Packages ####

library(sf)
library(dplyr)
library(readr)
library(ggplot2)

#0.2 Read in data ####

## apid obs for broward co ####
bees <- read.delim("../Data/Broward_apid_obs_unfiltered.csv") 

## shapefile for broward parks ####
parks <- st_read("../Data/urban_parks/urban_parks.shp")

#1.0 Geospatial Join ####
## make bee points in ST format ####
bees_sf <- bees %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

## transform bees to same CRS as parks ####
bees_sf <- st_transform(bees_sf, st_crs(parks))

## spatial join: attach park attributes to each bee point ####
bees_by_park <- st_join(bees_sf, parks, join = st_within)

## look at joined data ####
head(bees_by_park)


#2.0 Remove obs not in parks####
bees_in_parks <- bees_by_park %>% filter(!is.na(poly_id))

#3.0 Calculate percentage of apidae obs that are HB

Apid_Summary <- bees_in_parks %>%
  group_by(poly_id) %>%
  summarise(
    n_apidae = n(),
    n_hb = sum(genus == "Apis", na.rm = TRUE),
    percentageHB = 100 * n_hb / n_apidae
  )

Apid_Summary_20 <-  Apid_Summary %>% filter(n_apidae>5)

# Plot it ####

explore_Plot1 <- ggplot(
  data = Apid_Summary_20,
  aes(x = reorder(poly_id, percentageHB), y = percentageHB, color = percentageHB)
) +
  geom_point(aes(size=n_apidae)) +
  #scale_y_log10() +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(
    x = "Park",
    y = "Percent honey bee observations",
    color = "% honey bee"
  )

explore_Plot1

as.data.frame(Apid_Summary_20)

write.csv(
  as.data.frame(Apid_Summary_20)[,1:4],
  "../Output/ApisRelAbund.csv",
  row.names = FALSE
)
# Gutter ####

library(sf)
library(ggplot2)
library(dplyr)

ggplot() +
  geom_sf(data = parks, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data = bees_sf, alpha = 0.5, size = 0.8) +
  theme_minimal()

library(sf)
library(ggplot2)
library(dplyr)

bees_coords <- bees_sf %>%
  st_transform(3857) %>%   # projected CRS is better for spatial density
  cbind(st_coordinates(.))

ggplot() +
  geom_sf(data = st_transform(parks, 3857), fill = "grey95", color = "grey70", linewidth = 0.2) +
  stat_density_2d(
    data = bees_coords,
    aes(X, Y, fill = after_stat(level), alpha = after_stat(level)),
    geom = "polygon",
    contour = TRUE
  ) +
  geom_point(
    data = bees_coords,
    aes(X, Y),
    size = 0.4,
    alpha = 0.4
  ) +
  scale_alpha(range = c(0.05, 0.35), guide = "none") +
  coord_sf(crs = st_crs(3857)) +
  theme_minimal()


library(sf)
library(dplyr)
library(ggplot2)
library(dbscan)

bees_proj <- bees_sf %>% st_transform(3857)
coords <- st_coordinates(bees_proj)


ggplot() +
  geom_sf(data = st_transform(parks, 3857), fill = "grey95", color = "grey75", linewidth = 0.2) +
  geom_sf(data = bees_proj, aes(color = cluster), size = 1, alpha = 0.7) +
  theme_minimal()
