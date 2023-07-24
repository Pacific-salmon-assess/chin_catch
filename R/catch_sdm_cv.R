## Fit preliminary spatial distribution models
# August 7, 2020
# updated Feb 22, 2023
# updated July 11, 2023
# Ultimately estimate stock- or size-specific spatial distribution maps
# Here develop basic structure using cross validation to determine whether 
# temporal variability better characterized by week/month or year as FEs vs.
# spatiotemporal field fit only to total legal individuals
# Also include static spatial variables (e.g. depth/slope) and dynamic 
# fine-scale (e.g. time of day) and annual indices (e.g. PDO)


library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(raster)
library(sf)
library(sdmTMBextra)


chin <- readRDS(here::here("data", "clean_catch.RDS")) %>% 
  rename(agg = agg_name)


chin %>% 
  group_by(size_bin) %>% 
  tally()

# calculate total catch across size bins
catch_size1 <- chin %>% 
  group_by(event, size_bin) %>% 
  summarize(catch = n(), .groups = "drop") %>% 
  ungroup()


# clean and bind to set data
set_dat <- readRDS(here::here("data", "cleanSetData.RDS")) %>% 
  mutate(
    week = lubridate::week(date_time_local)
  )


catch_size <- expand.grid(
  event = set_dat$event,
  size_bin = unique(catch_size1$size_bin)
) %>%
  arrange(event) %>%
  left_join(., catch_size1, by = c("event", "size_bin")) %>%
  replace_na(., replace = list(catch = 0)) %>%
  left_join(., set_dat, by = "event") %>%
  # remove sets not on a troller and tacking
  filter(!grepl("rec", event),
         !tack == "yes") %>% 
  mutate(
    size_bin = as.factor(size_bin),
    # pool undersampled months
    month_f = case_when(
      month %in% c(4, 5) ~ 5,
      month %in% c(8, 9) ~ 8,
      TRUE ~ month
    ) %>% 
      as.factor(),
    year_f = as.factor(year),
    yUTM_ds = yUTM / 1000,
    xUTM_ds = xUTM / 1000,
    offset = log((set_dist / 1000) * lines),
    depth_z = scale(mean_depth)[ , 1],
    slope_z = scale(mean_slope)[ , 1],
    slack_z = scale(hours_from_slack)[ , 1],
    week_z = scale(week)[ , 1],
    moon_z = scale(moon_illuminated)[ , 1],
    dist_z = scale(coast_dist_km)
  ) 



## TOTAL CATCH SPATIOTEMPORAL MODELS -------------------------------------------

# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)


# test model
cv_6 <- sdmTMB_cv(
  catch ~ 0 + (1 | year_f) + size_bin + poly(slack_z, 2) +
    depth_z + poly(moon_z, 2) +
    slope_z + week_z:size_bin,
  offset = "offset",
  data = catch_size,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom1(),
  spatial = "off",
  spatial_varying = ~ 0 + size_bin + month_f,
  anisotropy = TRUE,
  control = sdmTMBcontrol(
    newton_loops = 1,
    map = list(
      ln_tau_Z = factor(
        c(1, 2, 3, rep(4, times = length(unique(catch_size$month_f)) - 1))
      )
    )
  ),
  silent = FALSE,
  use_initial_fit = TRUE,
  k_folds = 5
)

cv_8 <- sdmTMB_cv(
  catch ~ 0 + (1 | year_f) + size_bin + poly(slack_z, 2) +
    depth_z + poly(moon_z, 2) +
    slope_z + week_z:size_bin,
  offset = "offset",
  data = catch_size,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom1(),
  spatial = "off",
  time = "month",
  spatial_varying = ~ 0 + size_bin,
  spatiotemporal = "ar1",
  anisotropy = TRUE,
  control = sdmTMBcontrol(
    newton_loops = 1,
    map = list(
      ln_tau_Z = factor(
        c(1, 2, 3)
      )
    )
  ),
  silent = FALSE,
  use_initial_fit = TRUE,
  k_folds = 5
)

cv_9 <- sdmTMB_cv(
  catch ~ 0 + (1 | year_f) + size_bin + poly(slack_z, 2) +
    depth_z + poly(moon_z, 2) +
    slope_z + week_z:size_bin,
  offset = "offset",
  data = catch_size,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom1(),
  spatial = "off",
  spatial_varying = ~ 0 + size_bin + month_f + year_f,
  anisotropy = TRUE,
  control = sdmTMBcontrol(
    newton_loops = 1,
    map = list(
      ln_tau_Z = factor(
        c(1, 2, 3, rep(4, times = length(unique(catch_size$month_f)) - 1),
          rep(5, times = length(unique(catch_size$year_f)) - 1))
      )
    )
  ),
  silent = FALSE,
  use_initial_fit = TRUE,
  k_folds = 5
)

cv_list <- list(cv_6, cv_8, cv_9)
saveRDS(cv_list, here::here("data", "model_fits", "cv_list.rds"))

cv_list <- readRDS(here::here("data", "model_fits", "cv_list.rds"))

purrr::map(
  cv_list,
  ~ .x$elpd
)

purrr::map(
  cv_list,
  ~ .x$sum_loglik
)
