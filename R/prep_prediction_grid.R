
library(sf)
library(raster)
library(tidyverse)
library(ncdf4)
library(rmapshaper)
library(mapdata)
library(rnaturalearth)
library(rnaturalearthdata)


chin <- readRDS(here::here("data", "clean_catch.RDS"))

# shapefiles for IPES survey grid and combined WCVI/IPES grid
# ipes_grid_raw <- raster::shapefile(
#   here::here("data", "ipes_shapefiles", "IPES_Grid_UTM9.shp"))
# 
# ipes_grid_trim <- subset(ipes_grid_raw, LAT > 48.25 & LAT < 49.1) 

crit_hab <- raster::shapefile(
  here::here("data", "critical_habitat_trim2", 
             "Proposed_RKW_CriticalHabitat_SWVI_extended.shp")) %>% 
  # here::here("data", "critical_habitat_trim", 
  #            "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016_shrunk.shp")) %>% 
  sp::spTransform(., CRS("+proj=utm +zone=10 +units=m"))


# parallelize based on operating system (should speed up some of the spatial
# processing calculations)
# library("parallel")
# ncores <- detectCores() - 2
if (Sys.info()['sysname'] == "Windows") {
  # library("doParallel")
  # cl <- makeCluster(ncores)
  # registerDoParallel(cl)
  big_bathy_path <- "C:/Users/FRESHWATERC/Documents/drive docs/spatial/BC NetCDF"
  creel_path <- "C:/Users/FRESHWATERC/Documents/drive docs/spatial/creel_areas/"
} else {
  # doMC::registerDoMC(ncores)
  big_bathy_path <- "/Users/cam/Google Drive/spatial/BC NetCDF"
  creel_path <- "/Users/cam/Google Drive/spatial/creel_areas/"
}

ncin <- nc_open(
  paste(
    big_bathy_path, "gebco_2022_n58.9197_s47.033_w-137.9622_e-121.9747_ipes.nc",
    sep = "/"
  )
)

#specify lat/long
dep_list <- list(
  lon = ncvar_get(ncin, "lon"),
  lat = ncvar_get(ncin, "lat"),
  dep = ncvar_get(ncin, "elevation")
)

# function create dataframe
dep_dat_f <- function(x) {
  expand.grid(lon = x$lon, lat = x$lat) %>%
    as.matrix(.) %>%
    cbind(., depth = as.vector(x$dep)) %>%
    data.frame()
}

dep_dat_full <- dep_dat_f(dep_list) %>% 
  mutate(depth = -1 * depth) %>% 
  # remove land data
  filter(depth > 0)


# convert each to raster, downscale, and add terrain data to both
# then convert back to dataframes
bc_raster <- rasterFromXYZ(dep_dat_full, 
                           crs = "+proj=longlat +datum=WGS84")

bc_raster_utm <- projectRaster(
  bc_raster,
  crs = "+proj=utm +zone=10 +units=m",
  # convert to 1000 m resolution
  res = 1000
)

plot(bc_raster_utm)
plot(crit_hab, 
     add = T,
     border = "blue")


# crop to survey grid
# dum <- crop(bc_raster_utm, extent(ipes_grid_trim))
# ipes_raster_utm <- mask(dum, ipes_grid_trim)
study_area_utm <- crop(bc_raster_utm, extent(crit_hab))
ch_utm <- mask(dum, crit_hab)


# merge and add aspect/slope
study_area_slope <- terrain(
  study_area_utm, opt = 'slope', unit = 'degrees',
  neighbors = 8
)
ch_slope <- terrain(
  ch_utm, opt = 'slope', unit = 'degrees',
  neighbors = 8
)

# convert to SF
study_area_sf_list <- purrr::map(
  list(
    depth = study_area_utm,
    slope = study_area_slope
  ),
  function (x) {
    # leave at 1km x 1km res for inlets
    # aggregate(x, fact = 2) %>% 
    as(x, 'SpatialPixelsDataFrame') %>%
      as.data.frame() %>%
      st_as_sf(., coords = c("x", "y"),
               crs = sf::st_crs("+proj=utm +zone=10 +units=m"))
  }
)
ch_sf_list <- purrr::map(
  list(
    depth = ch_utm,
    slope = ch_slope
  ),
  function (x) {
    # leave at 1km x 1km res for inlets
    # aggregate(x, fact = 2) %>% 
    as(x, 'SpatialPixelsDataFrame') %>%
      as.data.frame() %>%
      st_as_sf(., coords = c("x", "y"),
               crs = sf::st_crs("+proj=utm +zone=10 +units=m"))
  }
)


# join depth and slope data
study_area_sf <- st_join(study_area_sf_list$depth, study_area_sf_list$slope) 
ch_sf <- st_join(ch_sf_list$depth, ch_sf_list$slope) 


# coast sf for plotting and calculating distance to coastline 
# (has to be lat/lon for dist2Line)
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -126.3, ymin = 48, 
              xmax = -125, ymax = 49.5) 

## convert to lat/lon for coast distance function 
ch_sf_deg <- ch_sf %>%
  st_transform(., crs = st_crs(coast))


# calculate distance to coastline
coast_dist <- geosphere::dist2Line(p = sf::st_coordinates(ch_sf_deg),
                                   line = as(coast, 'Spatial'))

coast_utm <- coast %>% 
  sf::st_transform(., crs = sf::st_crs("+proj=utm +zone=10 +units=m"))


# combine all data
ch_grid <- data.frame(
  st_coordinates(ch_sf[ , 1]),
  depth = ch_sf$depth,
  slope = ch_sf$slope,
  shore_dist = coast_dist[, "distance"]
)
study_grid <- data.frame(
  st_coordinates(study_area_sf[ , 1]),
  depth = study_area_sf$depth,
  slope = study_area_sf$slope
)

# interpolate missing data 
study_grid_interp <- VIM::kNN(study_grid, k = 5) %>% 
  dplyr::select(-ends_with("imp"))
ch_grid_interp <- VIM::kNN(ch_grid, k = 5) %>% 
  dplyr::select(-ends_with("imp"))


ggplot() + 
  geom_sf(data = coast_utm) +
  geom_raster(data = ch_grid, aes(x = X, y = Y, fill = depth)) +
  scale_fill_viridis_c() +
  ggsidekick::theme_sleek()


# export grid
saveRDS(study_grid_interp,
        here::here("data", "pred_bathy_grid_sa_1000m.RDS"))
saveRDS(ch_grid_interp,
        here::here("data", "pred_bathy_grid_1000m.RDS"))
saveRDS(coast_utm, here::here("data", "coast_trim_utm.RDS"))
