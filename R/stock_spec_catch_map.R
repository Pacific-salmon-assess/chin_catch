## Initial exploration of tagging data
# July 22, 2019

library(tidyverse)
library(ggplot2)
library(ggmap)
library(measurements)
library(maptools)
library(rmapshaper)
library(mapdata)

wCan <- map_data("worldHires", region = c("usa", "canada")) %>%
  filter(long < -110) %>%
  fortify(.)

sets <- readRDS(here::here("data", "cleanSetData.RDS")) %>% 
  mutate(month = lubridate::month(date_time_local))
chin <- readRDS(here::here("data", "cleanTagData_GSI.RDS")) %>% 
  mutate(month = lubridate::month(deployment_time),
         size_class = case_when(
           fl > 80 ~ "large",
           fl >= 68 ~ "med",
           fl < 65 ~ "small"
         ),
         size_class = as.factor(size_class),
         year = as.factor(year)) 

p <- ggplot(data = wCan, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-126.4, -125.25), ylim = c(48.4, 49), ratio = 1.3) + 
  geom_polygon(color = "black", fill = "gray80") +
  geom_point(data = sets, aes(x = lon, y = lat),
             inherit.aes = FALSE, shape = 3, alpha = 0.2) +
  labs(x = "Longitude", y = "Latitude") +
  ggsidekick::theme_sleek() +
  theme(axis.text = element_blank(),
        strip.background = element_rect(colour="white", fill="white")#,
        # legend.position=c(0.1, 0.2)
        ) +
  geom_point(data = chin %>% filter(!is.na(agg_name)),
             aes(x = lon, y = lat, fill = size_class),
             inherit.aes = FALSE, shape = 21) +
  scale_fill_viridis_d() +
  facet_grid(agg_name~as.factor(year))
p_no_smalls <- ggplot(data = wCan, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-126.4, -125.25), ylim = c(48.4, 49), ratio = 1.3) + 
  geom_polygon(color = "black", fill = "gray80") +
  geom_point(data = sets, aes(x = lon, y = lat),
             inherit.aes = FALSE, shape = 3, alpha = 0.2) +
  labs(x = "Longitude", y = "Latitude") +
  ggsidekick::theme_sleek() +
  theme(axis.text = element_blank(),
        strip.background = element_rect(colour="white", fill="white")#,
        # legend.position=c(0.1, 0.2)
  ) +
  geom_point(data = chin %>% filter(!is.na(agg_name), !size_class == "small"),
             aes(x = lon, y = lat, fill = size_class),
             inherit.aes = FALSE, shape = 21) +
  scale_fill_viridis_d() +
  facet_grid(agg_name~as.factor(year))

pdf(here::here("figs", "maps", "stock_specific_catch.pdf"), height = 7,
    width = 8)
p
p_no_smalls
dev.off()


# set locations only
sets_trim <- sets %>% 
  filter(charter == "troller")
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., xmin = -126.3, ymin = 48.45, xmax = -125.25, ymax = 49) 

bathy <- marmap::getNOAA.bathy(lon1 = -126.3, lon2 = -125.25, lat1 = 48.45, 
                               lat2 = 49, 
                       resolution = 0.1)  
trim_wcvi <- fortify(bathy) %>% 
  mutate(z = ifelse(z > 0 , 0, z),
         z2 = ifelse(z < -501, 501, -1 * z))

sets_out <- ggplot() + 
  geom_raster(data = trim_wcvi, aes(x = x, y = y, fill = z2)) +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_point(data = sets_trim, aes(x = lon, y = lat),
             inherit.aes = FALSE, shape = 21, alpha = 0.2, fill = "red") +
  ggsidekick::theme_sleek() +
  scale_fill_continuous(trans = "sqrt") +
  theme(strip.background = element_rect(colour="white", fill="white"),
        legend.position = "none",
        axis.title = element_blank()
  ) 

pdf(here::here("figs", "set_map.pdf"))
sets_out
dev.off()