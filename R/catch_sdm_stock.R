## Fit stock-based spatial distribution models
# July 5, 2023
# Uses parameter estimates generated in catch_sdm as priors for stock-specific
# models

library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(raster)
library(sf)
library(sdmTMBextra)


chin <- readRDS(here::here("data", "clean_catch.RDS")) %>% 
  rename(agg = agg_name)


chin %>% 
  group_by(agg) %>% 
  tally()


# calculate total catch across size bins
catch_stock1 <- chin %>%
  filter(!is.na(agg),
         fl > 65,
         # remove rare stocks
         !agg %in% c("ECVI", "WCVI", "WA_OR", "Cali")) %>%
  group_by(event, agg) %>%
  summarize(catch = n(), .groups = "drop") %>%
  ungroup()


# clean and bind to set data
set_dat <- readRDS(here::here("data", "cleanSetData.RDS")) %>% 
  mutate(
    week = lubridate::week(date_time_local)
  )


catch_stock <- expand.grid(
  event = set_dat$event,
  agg = unique(catch_stock1$agg)
) %>%
  arrange(event) %>%
  left_join(., catch_stock1, by = c("event", "agg")) %>%
  replace_na(., replace = list(catch = 0)) %>%
  left_join(., set_dat, by = "event") %>%
  # remove sets not on a troller and tacking
  filter(!grepl("rec", event),
         !tack == "yes") %>%
  mutate(
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
    moon_z = scale(moon_illuminated)[ , 1]
  )


# SPATIAL PREDICTION GRID ------------------------------------------------------

# import 1.8 X 1.8 km grid generated in prep_bathymetry.R in chinTagging repo
pred_bathy_grid <- readRDS(
  here::here("data", "pred_bathy_grid_utm.RDS"))

# generate combinations for month/year
pred_dat <- expand.grid(
  week = seq(min(catch_stock$week), max(catch_stock$week), length = 5),
  year_f = unique(catch_stock$year_f)
) %>% 
  filter(year_f == "2021") 

pred_grid_list <- vector(mode = "list", length = nrow(pred_dat))
for (i in seq_along(pred_grid_list)) {
  pred_grid_list[[i]] <- pred_bathy_grid %>% 
    filter(X < max(catch_stock$xUTM + 1000) & X > min(catch_stock$xUTM - 1000),
           Y < max(catch_stock$yUTM + 1000) & Y > min(catch_stock$yUTM - 1000)) %>% 
    mutate(week = pred_dat$week[i],
           year_f = pred_dat$year_f[i])
}

pred_grid <- pred_grid_list %>% 
  bind_rows() %>% 
  mutate(
    xUTM_ds = X / 1000,
    yUTM_ds = Y / 1000,
    slack_z = 0,
    moon_z = 0,
    month_f = case_when(
      week <= 22 ~ "5",
      week > 22 & week <= 26 ~ "6",
      week > 26 & week <= 30 ~ "7",
      week > 30 ~ "8"
    ),
    week_z = (week - mean(catch_stock$week)) / sd(catch_stock$week),
    depth_z = (depth - mean(catch_stock$mean_depth)) / 
      sd(catch_stock$mean_depth),
    slope_z = (slope - mean(catch_stock$mean_slope)) / 
      sd(catch_stock$mean_slope)
  ) 


# helper function for simple predictive maps
plot_map <- function(dat, column) {
  ggplot() +
    geom_raster(data = dat, aes(X, Y, fill = {{ column }})) +
    coord_fixed() +
    ggsidekick::theme_sleek() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank()
    )
}


## IMPACT OF PRIORS ------------------------------------------------------------

fr_sub <- catch_stock %>% filter(agg == "Fraser Sub.")

sdm_mesh1 <- make_mesh(fr_sub,
                       c("xUTM_ds", "yUTM_ds"),
                       n_knots = 250)

f1 <-  sdmTMB(
  catch ~ week_z,
  offset = "offset",
  data = fr_sub,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  anisotropy = FALSE,
  priors = sdmTMBpriors(
    phi = halfnormal(0, 5),
    matern_s = pc_matern(range_gt = 40, sigma_lt = 10)
  ),
  silent = FALSE
)

f2 <- sdmTMB(
  catch ~ week_z,
  offset = "offset",
  data = fr_sub,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  anisotropy = FALSE,
  priors = sdmTMBpriors(
    phi = halfnormal(0, 5),
    b = normal(
      #mean
      c(-2.5, 0.25),
      #sd
      c(2, 1)
    ),
    matern_s = pc_matern(range_gt = 40, sigma_lt = 10)
  ),
  silent = FALSE
)

f3 <- sdmTMB(
  catch ~ week_z + (1 | year_f),
  offset = "offset",
  data = fr_sub,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  anisotropy = FALSE,
  priors = sdmTMBpriors(
    phi = halfnormal(0, 5),
    b = normal(
      #mean
      c(-2.5, 0.25),
      #sd
      c(2, 1)
    ),
    matern_s = pc_matern(range_gt = 40, sigma_lt = 10)
  ),
  silent = FALSE
)

f4 <- sdmTMB(
  catch ~ week_z + (1 | year_f),
  offset = "offset",
  data = fr_sub,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom1(),
  spatial = "off",
  spatial_varying = ~ 0 + month_f,
  anisotropy = FALSE,
  priors = sdmTMBpriors(
    phi = halfnormal(0, 5),
    b = normal(
      #mean
      c(-2.5, 0.25),
      #sd
      c(2, 1)
    ),
    matern_s = pc_matern(range_gt = 40, sigma_lt = 10)
  ),
  control = sdmTMBcontrol(
    newton_loops = 1,
    map = list(
      ln_tau_Z = factor(
        rep(1, times = length(unique(catch_stock$month_f)))
      )
    )
  ),
  silent = FALSE
)

f5 <- sdmTMB(
  catch ~ week_z + (1 | year_f),
  offset = "offset",
  data = fr_sub,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  time = "month",
  spatiotemporal = "iid",
  anisotropy = FALSE,
  priors = sdmTMBpriors(
    phi = halfnormal(0, 5),
    b = normal(
      #mean
      c(-2.5, 0.25),
      #sd
      c(2, 1)
    ),
    matern_s = pc_matern(range_gt = 40, sigma_lt = 10),
    matern_st = pc_matern(range_gt = 40, sigma_lt = 10)
  ),
  silent = FALSE
)
