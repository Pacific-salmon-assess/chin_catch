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
  sf::st_crop(., xmin = -126.3, ymin = 48.45, xmax = -125.1, ymax = 49) %>% 
  sf::st_transform(., crs = sf::st_crs("+proj=utm +zone=10 +units=m"))


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

# pull coordinates to fix projection
coast_coords <- sf::st_coordinates(coast_trim)
crop_coast_coords <- sf::st_coordinates(crop_coast)

bath_grid <- readRDS(here::here("data", "pred_bathy_grid_sa_1000m.RDS")) %>%
  mutate(id = row_number()) #%>%
  # filter(X < max(sets$xUTM - 500),
  #        X > min(sets$xUTM + 500),
  #        Y < max(sets$yUTM + 500),
  #        Y > min(sets$yUTM - 500))

crit_hab <- sf::st_read(
  here::here("data", "critical_habitat_trim2", 
             "Proposed_RKW_CriticalHabitat_SWVI_extended.shp")) %>% 
  sf::st_transform(., crs = sf::st_crs("+proj=utm +zone=10 +units=m")) %>% 
  sf::st_crop(., 
              xmin = min(sets$xUTM + 500), 
              ymin = min(sets$yUTM - 500), 
              xmax = max(sets$xUTM - 500), 
              ymax = max(sets$yUTM + 500))


main_map <- base_map +
  geom_sf(data = coast_trim, fill = "darkgrey") +
  geom_sf(data = crit_hab, colour = "red", fill = "transparent", lty = 2) +
  theme(axis.text = element_blank()) +
  scale_x_continuous(limits = c(min(coast_coords[ , "X"]), 
                                max(coast_coords[ , "X"]) - 10000), 
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(coast_coords[ , "Y"]) + 5000, 
                                max(coast_coords[ , "Y"]) - 15000), 
                     expand = c(0, 0))

sets_only <- ggplot() +
  geom_sf(data = crop_coast, color = "black", fill = "white", size = 1.25) +
  geom_point(data = sets, aes(x = xUTM, y = yUTM, fill = as.factor(year)),
             inherit.aes = FALSE, alpha = 0.4, shape = 21) +
  geom_contour(data = bath_grid %>% 
                 filter(depth < 401), 
               aes(x = X, y = Y, z = depth), 
               breaks = seq(50, 400, by = 50),
               colour = "black") +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(guide = "none") +
  theme(axis.title = element_blank(),
        legend.title = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "grey60")) +
  scale_x_continuous(limits = c(min(crop_coast_coords[ , "X"]) + 2500, 
                                max(crop_coast_coords[ , "X"]) - 5000), 
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(crop_coast_coords[ , "Y"]) + 4500, 
                                max(crop_coast_coords[ , "Y"]) - 2500), 
                     expand = c(0, 0))


png(here::here("figs", "ms_figs", "study_area1.png"), res = 250, units = "in", 
    height = 5.5, width = 7)
cowplot::ggdraw() +
  cowplot::draw_plot(sets_only) +
  cowplot::draw_plot(main_map,
            height = 0.33, x = -0.33, y = 0.08)
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
             aes(x = xUTM, y = yUTM),
             shape = 3, alpha = 0.1,
             inherit.aes = FALSE) +
  geom_point(data = pos_size,
             aes(x = xUTM, y = yUTM, size = catch),
             shape = 21, fill = "#377eb8", alpha = 0.75,
             inherit.aes = FALSE) +
  facet_grid(size_bin~month_f) +
  theme(legend.position = "top",
        axis.text = element_blank(),
        panel.background = element_rect(fill = "grey60")) +
  scale_x_continuous(limits = c(min(crop_coast_coords[ , "X"]) + 5000, 
                                max(crop_coast_coords[ , "X"]) - 5000), 
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(crop_coast_coords[ , "Y"]) + 4500, 
                                max(crop_coast_coords[ , "Y"]) - 2500), 
                     expand = c(0, 0))

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
             aes(x = xUTM, y = yUTM),
             shape = 3, alpha = 0.1,
             inherit.aes = FALSE) +
  geom_point(data = pos_stock,
             aes(x = xUTM, y = yUTM, size = catch),
             shape = 21, fill = "#4daf4a", alpha = 0.75,
             inherit.aes = FALSE) +
  facet_grid(agg~month_f) +
  theme(legend.position = "top",
        axis.text = element_blank(),
        panel.background = element_rect(fill = "grey60")) +
  scale_x_continuous(limits = c(min(crop_coast_coords[ , "X"]) + 5000, 
                                max(crop_coast_coords[ , "X"]) - 5000), 
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(crop_coast_coords[ , "Y"]) + 4500, 
                                max(crop_coast_coords[ , "Y"]) - 2500), 
                     expand = c(0, 0))

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
             aes(x = xUTM, y = yUTM),
             shape = 3, alpha = 0.1,
             inherit.aes = FALSE) +
  geom_point(data = pos_origin,
             aes(x = xUTM, y = yUTM, size = catch),
             shape = 21, fill = "#984ea3", alpha = 0.75,
             inherit.aes = FALSE) +
  facet_grid(origin~month_f) +
  theme(legend.position = "top",
        axis.text = element_blank(),
        panel.background = element_rect(fill = "grey60")) +
  scale_x_continuous(limits = c(min(crop_coast_coords[ , "X"]) + 5000, 
                                max(crop_coast_coords[ , "X"]) - 5000), 
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(crop_coast_coords[ , "Y"]) + 4500, 
                                max(crop_coast_coords[ , "Y"]) - 2500), 
                     expand = c(0, 0))

png(here::here("figs", "ms_figs", "origin_map.png"), res = 250, units = "in", 
    height = 4.25, width = 7.5)
origin_map
dev.off()
