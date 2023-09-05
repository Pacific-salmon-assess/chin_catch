
library(sf)
library(raster)
# library(rgdal)
library(tidyverse)
library(ncdf4)
# library(maptools)
library(rmapshaper)
library(mapdata)
library(rnaturalearth)
library(rnaturalearthdata)


chin <- readRDS(here::here("data", "clean_catch.RDS")) %>% 
  rename(agg = agg_name)

# shapefiles for IPES survey grid and combined WCVI/IPES grid
# ipes_grid_raw <- raster::shapefile(
#   here::here("data", "ipes_shapefiles", "IPES_Grid_UTM9.shp"))
# 
# ipes_grid_trim <- subset(ipes_grid_raw, LAT > 48.25 & LAT < 49.1) 

crit_hab <- raster::shapefile(
  here::here("data", "critical_habitat", 
             "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016.shp")) %>% 
  sp::spTransform(., CRS("+proj=utm +zone=9 +units=m"))



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
                           crs = sp::CRS("+proj=longlat +datum=WGS84"))
bc_raster_utm <- projectRaster(bc_raster,
                               crs = sp::CRS("+proj=utm +zone=9 +units=m"),
                               # convert to 1000 m resolution
                               res = 1000)


plot(bc_raster_utm)
plot(crit_hab, 
     add = T,
     border = "blue")
plot(ipes_grid_trim, 
     add = T,
     border = "blue")


# crop to survey grid
# dum <- crop(bc_raster_utm, extent(ipes_grid_trim))
# ipes_raster_utm <- mask(dum, ipes_grid_trim)
dum <- crop(bc_raster_utm, extent(crit_hab))
ch_utm <- mask(dum, crit_hab)


# merge and add aspect/slope
ch_slope <- terrain(
  ch_utm, opt = 'slope', unit = 'degrees',
  neighbors = 8
)
ch_aspect <- terrain(
  ch_utm, opt = 'aspect', unit = 'degrees',
  neighbors = 8
)
ch_list <- list(
  depth = ch_utm,
  slope = ch_slope,
  aspect = ch_aspect
)



#TODO:

ipes_sf_list <- purrr::map(
  ipes_raster_list,
  function (x) {
    # leave at 1km x 1km res for inlets
    # aggregate(x, fact = 2) %>% 
    as(x, 'SpatialPixelsDataFrame') %>%
      as.data.frame() %>%
      st_as_sf(., coords = c("x", "y"),
               crs = sp::CRS("+proj=utm +zone=9 +units=m"))
  }
)


# join depth and slope data
ipes_sf <- st_join(ipes_sf_list$depth, ipes_sf_list$slope) 

# coast sf for plotting and calculating distance to coastline 
# (has to be lat/lon for dist2Line)
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = min(dat_trim$lon), ymin = 48, 
              xmax = max(dat_trim$lon), ymax = max(dat_trim$lat)) 

## convert to lat/lon for coast distance function 
ipes_sf_deg <- ipes_sf %>%
  st_transform(., crs = st_crs(coast))


# calculate distance to coastline
coast_dist <- geosphere::dist2Line(p = sf::st_coordinates(ipes_sf_deg),
                                   line = as(coast, 'Spatial'))

coast_utm <- coast %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))


# combine all data
ipes_grid <- data.frame(
  st_coordinates(ipes_sf[ , 1]),
  depth = ipes_sf$depth,
  slope = ipes_sf$slope,
  shore_dist = coast_dist[, "distance"]
)
# ipes_grid_trim <- ipes_grid %>% filter(!depth > 500)

# interpolate missing data 
ipes_grid_interp <- VIM::kNN(ipes_grid, k = 5)


ggplot() + 
  geom_sf(data = coast_utm) +
  geom_raster(data = ipes_grid, aes(x = X, y = Y, fill = depth)) +
  scale_fill_viridis_c() +
  # geom_point(data = dat_trim, aes(x = utm_x, y = utm_y),
  #            fill = "white",
  #            shape = 21) +
  ggsidekick::theme_sleek()


# subset WCVI predictive grid to make IPES only
ipes_grid_raw_sf <- st_read(
  here::here("data", "spatial", "ipes_shapefiles", "IPES_Grid_UTM9.shp"))

ipes_wcvi_sf <- st_as_sf(ipes_grid_interp %>% 
                           select(-ends_with("imp")),
                         coords = c("X", "Y"),
                         crs = st_crs(ipes_grid_raw_sf))

dd <- st_intersection(ipes_grid_raw_sf, ipes_wcvi_sf)

ipes_grid_only <- data.frame(
  st_coordinates(dd[ , 1]),
  depth = dd$depth,
  slope = dd$slope,
  shore_dist = dd$shore_dist
)

grid_list <- list(
  ipes_grid = ipes_grid_only,
  wcvi_grid = ipes_grid_interp %>% 
    select(-ends_with("imp"))
)


# export grid
saveRDS(grid_list,
        here::here("data", "spatial", "pred_ipes_grid.RDS"))
saveRDS(coast_utm, here::here("data", "spatial", "coast_trim_utm.RDS"))
