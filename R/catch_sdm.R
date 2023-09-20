## Fit size-based spatial distribution models
# August 7, 2020
# updated Feb 22, 2023
# Given cross-validation (catch_sdm_cv.R) and goals of analysis, utilize
# month focused spatio temporal model
# Also include static spatial variables (e.g. depth/slope) and dynamic 
# fine-scale (e.g. time of day) 
# Fit separate models to each size class, unlike _mvrw which fits multivariate
# random field random walk to all sizes simultaneously

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


# import bycatch data representing sublegal catch
chin_juv <- read.csv(
  here::here("data", "bycatchData.csv"),
  stringsAsFactors = F) %>% 
  filter(species == "chinook",
         !grepl("2023", event)) %>% 
  mutate(size_bin = "sublegal") %>% 
  dplyr::select(
    event, size_bin, catch = count
  )


# calculate total catch across size bins
catch_size1 <- chin %>% 
  group_by(event, size_bin) %>% 
  summarize(catch = n(), .groups = "drop") %>% 
  ungroup() %>% 
  rbind(., chin_juv)


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
    month = case_when(
      month %in% c(4, 5) ~ 5,
      # month %in% c(8, 9) ~ 8,
      TRUE ~ month
    ),
    month_f = month %>% 
      as.factor(),
    year_f = as.factor(year),
    yUTM_ds = yUTM / 1000,
    xUTM_ds = xUTM / 1000,
    offset = log((set_dist / 1000) * lines),
    # effort of one corresponds to ~2 lines towed for 200 m
    depth_z = scale(mean_depth)[ , 1],
    slope_z = scale(mean_slope)[ , 1],
    slack_z = scale(hours_from_slack)[ , 1],
    week_z = scale(week)[ , 1],
    moon_z = scale(moon_illuminated)[ , 1],
    dist_z = scale(coast_dist_km)[ , 1]
  ) 

catch_tbl <- catch_size %>% 
  group_by(size_bin) %>% 
  group_nest()
 

# BUILD MESHES -----------------------------------------------------------------

# construct a common mesh
sdm_mesh1 <- make_mesh(catch_tbl$data[[1]],
                      c("xUTM_ds", "yUTM_ds"),
                      n_knots = 125)



# SPATIAL PREDICTION GRID ------------------------------------------------------

# import 1 x 1 km grid based on shrunk CH shape files and generated in 
# prep_prediction_grid.R 
pred_bathy_grid <- readRDS(
  here::here("data", "pred_bathy_grid_1000m.RDS"))

# generate combinations for month/year
pred_dat <- expand.grid(
  # week = seq(min(catch_size$week), max(catch_size$week), length = 5),
  week = c(18.5, 23, 27.5, 32, 36.5),
  year_f = unique(catch_size$year_f)
) %>% 
  filter(year_f == "2021") 

pred_grid_list <- vector(mode = "list", length = nrow(pred_dat))
for (i in seq_along(pred_grid_list)) {
  pred_grid_list[[i]] <- pred_bathy_grid %>% 
    filter(X < max(catch_size$xUTM + 1000) & X > min(catch_size$xUTM - 1000),
           Y < max(catch_size$yUTM + 1000) & Y > min(catch_size$yUTM - 1000)) %>% 
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
    month = case_when(
      week < 23 ~ 5,
      week >= 23 & week < 27 ~ 6,
      week >= 27 & week < 32 ~ 7,
      week >= 32 & week < 35 ~ 8,
      week >= 35 ~ 9
    ),
    month_f = as.factor(month),
    week_z = (week - mean(catch_size$week)) / sd(catch_size$week),
    depth_z = (depth - mean(catch_size$mean_depth)) / 
      sd(catch_size$mean_depth),
    slope_z = (slope - mean(catch_size$mean_slope)) / 
      sd(catch_size$mean_slope)
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


# SPATIAL MODELS ---------------------------------------------------------------

fit_list <- purrr::map(
  catch_tbl$data,
  ~ sdmTMB(
    catch ~ (1 | year_f) + poly(week_z, 2),
    offset = "offset",
    data = .x,
    mesh = sdm_mesh1,
    family = sdmTMB::nbinom1(),
    spatial = "on",
    # spatial_varying = ~ 0 + size_bin,
    time = "month",
    spatiotemporal = "rw",
    anisotropy = FALSE,
    # share_range = FALSE,
    silent = FALSE
  )
)
fit_m <- sdmTMB(
  catch ~ (1 | year_f) + poly(week_z, 2),
  offset = "offset",
  data = catch_tbl$data[[2]],
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  # spatial_varying = ~ 0 + size_bin,
  time = "month",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  # share_range = FALSE,
  silent = FALSE
)

## fails to converge for full so focus on medium size class for comp
saveRDS(fit_list[[2]], here::here("data", "model_fits", "fit_medium.rds"))


fit_m <- readRDS(here::here("data", "model_fits", "fit_medium.rds"))



## SPATIAL PREDICTIONS ---------------------------------------------------------

upsilon_pred <- predict(fit_m,
                        newdata = pred_grid[1:3, ])

pp <- predict(fit, newdata = pred_grid, se_fit = FALSE, re_form = NULL)
# for now remove extra rows but check w/ Sean
pp$upsilon_stc <- as.numeric(upsilon_pred$proj_upsilon_st_A_vec)[1:nrow(pp)]

week_key <- data.frame(
  week = unique(pred_grid$week)
) %>% 
  mutate(
    date = c("May 1", "June 5", "July 10", "Aug 15", "Sep 10") %>% 
      as.factor() %>% 
      fct_reorder(., week)
  )
  

size_omega <- pp %>%
  dplyr::select(-c(size_bin, week, month)) %>%
  distinct() %>% 
  plot_map(., omega_s) +
  scale_fill_gradient2(name = "Spatial\nRF Effect") +
  theme(legend.position = "top")

month_omega <- pp %>%
  plot_map(., upsilon_stc) +
  scale_fill_gradient2(name = "Spatial\nRF Effect") +
  facet_grid(size_bin~month) +
  theme(legend.position = "top")

full_preds <- pp %>% 
  group_by(size_bin) %>% 
  mutate(scale_est = exp(est) / max(exp(est))) %>%
  left_join(., week_key, by = "week") %>% 
  plot_map(., scale_est) +
  scale_fill_viridis_c(trans = "sqrt", name = "Scaled\nAbundance") +
  facet_grid(date~size_bin)


png(here::here("figs", "ms_figs", "month_omega.png"), res = 250, units = "in", 
    height = 4.5, width = 4.5)
month_omega
dev.off()

png(here::here("figs", "ms_figs", "size_omega.png"), res = 250, units = "in", 
    height = 7.5, width = 7.5)
size_omega
dev.off()

png(here::here("figs", "ms_figs", "spatial_preds.png"), res = 250, units = "in", 
    height = 7.5, width = 7.5)
full_preds
dev.off()

