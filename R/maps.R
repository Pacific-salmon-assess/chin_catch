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

chin <- readRDS(here::here("data", "cleanTagData_GSI.RDS")) %>% 
  mutate(
    month = lubridate::month(deployment_time),
    month = case_when(
      month %in% c(4, 5) ~ 5,
      month %in% c(8, 9) ~ 8,
      TRUE ~ month
    ),
    month_f = as.factor(month),
    size_class = case_when(
      fl < 65 ~ "small",
      fl >= 65 & fl < 75 ~ "medium",
      fl >= 75 ~ "large"
    ),
    size_class = as.factor(size_class),
    year = as.factor(year),
    agg_name = fct_relevel(agg_name, "WCVI", "ECVI", "Fraser", "PugetSo",
                           "Col", "WA_OR", "Cali")
  ) %>% 
  filter(event %in% sets$event)


## STUDY AREA ------------------------------------------------------------------

wCan <- map_data("worldHires", region = c("usa", "canada")) %>%
  filter(long < -110) %>%
  fortify(.)

crop_coast <- readRDS(here::here("data", "crop_coast_sf.rds")) %>% 
  st_crop(., xmin = -126.3, ymin = 48.45, xmax = -125.15, ymax = 49) 


sets_only <- ggplot() +
  geom_sf(data = crop_coast, color = "black", fill = "white", size = 1.25) +
  geom_point(data = sets, aes(x = lon, y = lat, colour = as.factor(year)),
             inherit.aes = FALSE, alpha = 0.2) +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title = element_blank(),
        legend.title = element_blank()) 

stock_catch <- ggplot() +
  geom_sf(data = crop_coast, color = "black", fill = "white", size = 1.25) +
  geom_point(data = chin,
             aes(x = lon, y = lat, fill = size_class),
             inherit.aes = FALSE, shape = 21, alpha = 0.5) +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_d() +
  theme(axis.title = element_blank(),
        legend.title = element_blank(),
        axis.text = element_blank()) +
  facet_grid(agg_name~month)


png(here::here("figs", "ms_figs", "set_map.png"), res = 250, units = "in", 
    height = 5.5, width = 5.5)
sets_only
dev.off()

png(here::here("figs", "ms_figs", "stock_catch_map.png"), 
    res = 250, units = "in", height = 5.5, width = 5.5)
stock_catch
dev.off()