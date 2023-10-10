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

# chin <- readRDS(here::here("data", "cleanTagData_GSI.RDS")) %>% 
#   mutate(
#     month = lubridate::month(deployment_time),
#     month = case_when(
#       month %in% c(4, 5) ~ 5,
#       month %in% c(8, 9) ~ 8,
#       TRUE ~ month
#     ),
#     month_f = as.factor(month),
#     size_class = case_when(
#       fl < 65 ~ "small",
#       fl >= 65 & fl < 75 ~ "medium",
#       fl >= 75 ~ "large"
#     ),
#     size_class = as.factor(size_class),
#     year = as.factor(year),
#     agg_name = fct_relevel(agg_name, "WCVI", "ECVI", "Fraser", "PugetSo",
#                            "Col", "WA_OR", "Cali")
#   ) %>% 
#   filter(event %in% sets$event)


## STUDY AREA ------------------------------------------------------------------

set_dat <- readRDS(here::here("data", "cleanSetData.RDS")) 

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
  filter(X < max(set_dat$xUTM - 500),
         X > min(set_dat$xUTM + 1000),
         Y < max(set_dat$yUTM - 500),
         Y > min(set_dat$yUTM + 500))

crit_hab <- sf::st_read(
  here::here("data", "critical_habitat_trim", 
             "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016_shrunk.shp")) %>% 
  sf::st_transform(., crs = sf::st_crs("+proj=utm +zone=10 +units=m")) %>% 
  sf::st_crop(., 
              xmin = min(set_dat$xUTM + 1000), 
              ymin = max(set_dat$yUTM - 500), 
              xmax = max(set_dat$xUTM - 500), 
              ymax = min(set_dat$yUTM + 500))


base_map <- ggplot() + 
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank(),
        legend.position = "top") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
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

crop_coast <- readRDS(here::here("data", "crop_coast_sf.RDS")) %>% 
# crop_coast <- coast %>%
  sf::st_crop(., xmin = -126.3, ymin = 48.45, xmax = -125.1, ymax = 49)

sets_only <- ggplot() +
  geom_sf(data = crop_coast, color = "black", fill = "white", size = 1.25) +
  geom_point(data = sets, aes(x = lon, y = lat, fill = as.factor(year)),
             inherit.aes = FALSE, alpha = 0.4, shape = 21) +
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

# stock_catch <- ggplot() +
#   geom_sf(data = crop_coast, color = "black", fill = "white", size = 1.25) +
#   geom_point(data = chin,
#              aes(x = lon, y = lat, fill = size_class),
#              inherit.aes = FALSE, shape = 21, alpha = 0.5) +
#   ggsidekick::theme_sleek() +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_fill_viridis_d() +
#   theme(axis.title = element_blank(),
#         legend.title = element_blank(),
#         axis.text = element_blank()) +
#   facet_grid(agg_name~month)
# 
# 
# png(here::here("figs", "ms_figs", "stock_catch_map.png"), 
#     res = 250, units = "in", height = 5.5, width = 5.5)
# stock_catch
# dev.off()