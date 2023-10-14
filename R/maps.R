## Map Figures
# June 28, 2023

library(tidyverse)
library(ggplot2)
library(ggmap)
library(measurements)
library(maptools)
library(rmapshaper)
library(mapdata)


sets <- readRDS(here::here("data", "cleanSetData.RDS")) %>% 
  mutate(month = lubridate::month(date_time_local)) %>%
  # remove sets not on a troller
  filter(!grepl("rec", event))

# catch data pooled by tow
catch_size <- readRDS(here::here("data", "catch_size_pre.rds"))
catch_stock <- readRDS(here::here("data", "catch_stock_pre.rds"))
catch_origin <- readRDS(here::here("data", "catch_origin_pre.rds"))

# coastline file
crop_coast <- readRDS(here::here("data", "crop_coast_sf.RDS")) %>% 
  sf::st_crop(., xmin = -126.3, ymin = 48.45, xmax = -125.1, ymax = 49)

# base map used for study area and catches
base_map <- ggplot() + 
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank(),
        legend.position = "top") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))


## STUDY AREA ------------------------------------------------------------------

# map data
coast <- rbind(rnaturalearth::ne_states( "United States of America",
                                         returnclass = "sf"),
               rnaturalearth::ne_states( "Canada", returnclass = "sf"))
coast_trim <- sf::st_crop(
  coast,
  xmin = -128.5, ymin = 48, xmax = -122.5, ymax = 51.5
) %>% 
  sf::st_transform(., crs = sf::st_crs("+proj=utm +zone=10 +units=m"))

bath_grid <- readRDS(here::here("data", "pred_bathy_grid_1000m.RDS")) %>% 
  mutate(id = row_number(),
         shore_dist_km = shore_dist / 1000) %>% 
  filter(X < max(sets$xUTM - 500),
         X > min(sets$xUTM + 500),
         Y < max(sets$yUTM + 500),
         Y > min(sets$yUTM - 500))

crit_hab <- sf::st_read(
  here::here("data", "critical_habitat_trim", 
             "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016_shrunk.shp")) %>% 
  sf::st_transform(., crs = sf::st_crs("+proj=utm +zone=10 +units=m")) %>% 
  sf::st_crop(., 
              xmin = min(sets$xUTM + 500), 
              ymin = min(sets$yUTM - 500), 
              xmax = max(sets$xUTM - 500), 
              ymax = max(sets$yUTM + 500))


main_map <- base_map +
  geom_sf(data = coast_trim, fill = "darkgrey") +
  geom_sf(data = crit_hab, colour = "red", fill = "transparent", lty = 2) #+
  # theme(axis.text = element_blank())
depth_map <- base_map +
  geom_raster(data = bath_grid, aes(x = X, y = Y, fill = depth)) +
  scale_fill_viridis_c(name = "Bottom\nDepth (m)", direction = -1)  +
  theme(legend.position = "right",
        axis.text = element_blank())
slope_map <- base_map +
  geom_raster(data = bath_grid, aes(x = X, y = Y, fill = slope)) +
  scale_fill_viridis_c(name = "Bottom\nSlope\n(degrees)", direction = -1)  +
  theme(legend.position = "right",
        axis.text = element_blank())

crit_hab_lat <- crit_hab %>% 
  sf::st_transform(., crs = sf::st_crs(crop_coast))

sets_only <- ggplot() +
  geom_sf(data = crop_coast, color = "black", fill = "white", size = 1.25) +
  geom_point(data = sets, aes(x = lon, y = lat, fill = as.factor(year)),
             inherit.aes = FALSE, alpha = 0.4, shape = 21) +
  geom_sf(data = crit_hab, colour = "red", fill = "transparent", lty = 2) +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(guide = "none") +
  theme(axis.title = element_blank(),
        legend.title = element_blank(),
        axis.text = element_blank()) 

p1 <- cowplot::plot_grid(
  depth_map, slope_map,
  ncol = 1
)
p2 <- cowplot::plot_grid(
  main_map, sets_only,
  ncol = 1
)
p3 <- cowplot::plot_grid(
  plotlist = list(p2, p1),
  ncol = 2,
  rel_widths = c(1, 1.2)
)

png(here::here("figs", "ms_figs", "study_area.png"), res = 250, units = "in", 
    height = 5.5, width = 7.5)
p3
dev.off()


## CATCH MAPS ------------------------------------------------------------------

# size
blank_size <- catch_size %>% 
  filter(catch == "0")
pos_size <- catch_size %>% 
  filter(!catch == "0")

size_map <- base_map + 
  geom_sf(data = crop_coast, color = "black", fill = "white", size = 1.25) +
  geom_point(data = blank_size, 
             aes(x = lon, y = lat),
             shape = 3, alpha = 0.1,
             inherit.aes = FALSE) +
  geom_point(data = pos_size,
             aes(x = lon, y = lat, size = catch),
             shape = 21, fill = "#377eb8", alpha = 0.75,
             inherit.aes = FALSE) +
  facet_grid(size_bin~month_f) +
  theme(legend.position = "top",
        axis.text = element_blank())

png(here::here("figs", "ms_figs", "size_map.png"), res = 250, units = "in", 
    height = 5.5, width = 7.5)
size_map
dev.off()


# stock
blank_stock <- catch_stock %>% 
  filter(catch == "0")
pos_stock <- catch_stock %>% 
  filter(!catch == "0")

stock_map <- base_map + 
  geom_sf(data = crop_coast, color = "black", fill = "white", size = 1.25) +
  geom_point(data = blank_stock, 
             aes(x = lon, y = lat),
             shape = 3, alpha = 0.1,
             inherit.aes = FALSE) +
  geom_point(data = pos_stock,
             aes(x = lon, y = lat, size = catch),
             shape = 21, fill = "#4daf4a", alpha = 0.75,
             inherit.aes = FALSE) +
  facet_grid(agg~month_f) +
  theme(legend.position = "top",
        axis.text = element_blank())

png(here::here("figs", "ms_figs", "stock_map.png"), res = 250, units = "in", 
    height = 8.5, width = 7.5)
stock_map
dev.off()


# origin
blank_origin <- catch_origin %>% 
  filter(catch == "0")
pos_origin <- catch_origin %>% 
  filter(!catch == "0")

origin_map <- base_map + 
  geom_sf(data = crop_coast, color = "black", fill = "white", size = 1.25) +
  geom_point(data = blank_origin, 
             aes(x = lon, y = lat),
             shape = 3, alpha = 0.1,
             inherit.aes = FALSE) +
  geom_point(data = pos_origin,
             aes(x = lon, y = lat, size = catch),
             shape = 21, fill = "#984ea3", alpha = 0.75,
             inherit.aes = FALSE) +
  facet_grid(origin~month_f) +
  theme(legend.position = "top",
        axis.text = element_blank())

png(here::here("figs", "ms_figs", "origin_map.png"), res = 250, units = "in", 
    height = 4.25, width = 7.5)
origin_map
dev.off()
