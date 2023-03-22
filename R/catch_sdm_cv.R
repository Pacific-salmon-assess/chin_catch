## Fit preliminary spatial distribution models
# August 7, 2020
# updated Feb 22, 2023
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


# Import stage data and fitted model generated in gen_detection_histories.R
stage_dat <- readRDS(here::here("data", "generated_data", 
                                "agg_lifestage_df.RDS")) %>% 
  dplyr::select(vemco_code, stage)


chin_raw <- readRDS(here::here("data", "tagging_data", 
                               "cleanTagData_GSI.RDS")) %>%
  mutate(year_day = lubridate::yday(date)) %>% 
  rename(vemco_code = acoustic_year, agg = agg_name) %>% 
  left_join(., stage_dat, by = "vemco_code") 

## predict life stage based on fitted model 
stage_mod <- readRDS(here::here("data", "model_outputs", "stage_fl_hierA.RDS"))
# include RIs when stock ID is known
pred_dat <- tidybayes::add_epred_draws(
  object = stage_mod,
  newdata = chin_raw %>% 
    filter(is.na(stage) | stage == "unknown",
           !is.na(agg)),
  re_formula = NULL, #scale = "response", 
  ndraws = 100
)
# otherwise marginalize RIs
pred_dat_na <- tidybayes::add_epred_draws(
  object = stage_mod,
  newdata = chin_raw %>% 
    filter(is.na(stage) | stage == "unknown",
           is.na(agg)),
  re_formula = NA, #scale = "response", 
  ndraws = 100
)

fl_preds_mean <- rbind(pred_dat, pred_dat_na) %>% 
  group_by(fish) %>% 
  summarize(med = median(.epred), 
            lo = quantile(.epred, prob = 0.05),
            up = quantile(.epred, prob = 0.95),
            .groups = "drop") %>% 
  mutate(stage_predicted = ifelse(med > 0.50 & lo > 0.5, "mature", "immature"))

chin <- left_join(chin_raw, fl_preds_mean, by = "fish") %>% 
  mutate(
    stage = ifelse(is.na(stage) | stage == "unknown", stage_predicted, stage),
    month = lubridate::month(date)
  ) %>% 
  dplyr::select(fish, vemco_code, event, date, year, month, year_day,
                fl, clip, stage, stock:agg_prob) %>% 
  filter(stage == "mature")
    
# calculate total catch 
catch1 <- chin %>% 
  group_by(event) %>% 
  tally(name = "abund") %>% 
  mutate(agg = "total") %>% 
  ungroup()


# clean and bind to set data
set_dat <- readRDS(here::here("data", "tagging_data", "cleanSetData.RDS")) %>% 
  mutate(
    week = lubridate::week(date_time_local)
  )


# import lighthouse data
sst <- readRDS(here::here("data", "spatial_data", "trim_amphitrite.rds"))


catch <- expand.grid(event = set_dat$event) %>% 
  arrange(event) %>% 
  left_join(., catch1, by = c("event")) %>%
  replace_na(., replace = list(abund = 0)) %>% 
  left_join(., set_dat, by = "event") %>%
  left_join(., sst, by = c("month", "year")) %>% 
  mutate(
    year_f = as.factor(year),
    yUTM_ds = yUTM / 1000,
    xUTM_ds = xUTM / 1000,
    offset = log(set_dist),
    sst_z = scale(sst)[ , 1],
    week_z = scale(week)[ , 1]
  ) %>% 
  # remove sets not on a troller
  filter(!grepl("rec", event))



# FIT SIMPLE GAMS --------------------------------------------------------------

## check correlations among explanatory variables
# corr <- catch_trim %>% 
#   dplyr::select(year_day_z,
#                 mean_depth_z,
#                 mean_slope_z,
#                 tide_z,
#                 hr_slack_z,
#                 moon_z,
#                 coast_dist_z) %>% 
#   cor()
# ggcorrplot::ggcorrplot(corr)
# no red flags


## explore 
ggplot(catch) +
  geom_point(aes(x = year_day, y = abund))
ggplot(catch) +
  geom_point(aes(x = mean_depth, y = abund))
ggplot(catch) +
  geom_point(aes(x = mean_slope, y = abund))
ggplot(catch) +
  geom_point(aes(x = hours_from_slack, y = abund))
ggplot(catch) +
  geom_point(aes(x = moon_illuminated, y = abund))
ggplot(catch) +
  geom_point(aes(x = coast_dist, y = abund))


## explore relationships between depth, sd depth, yday, hour, and time to slack
# to determine whether GAMs are appropriate and which FEs to include
# ambiguous support for including both hr to slack and tide, remove latter for
# now for simplicity
mod1 <- mgcv::gam(abund ~ s(year_f, bs = "re") + s(year_day, bs = "tp") + 
                   s(mean_depth, bs = "tp") + 
                   s(mean_slope, bs = "tp") + #s(tide, bs = "tp") + 
                   s(hours_from_slack, bs = "tp") + 
                   s(moon_illuminated, bs = "tp") + 
                   offset(offset),
                 data = catch,
                 family = mgcv::nb(link = "log"),
                 method = "ML")

# alternative model uses lower resolution week instead of year_day
mod2 <- mgcv::gam(abund ~ s(year_f, bs = "re") + s(week, bs = "tp") + 
                    s(mean_depth, bs = "tp") + 
                    s(mean_slope, bs = "tp") + #s(tide, bs = "tp") + 
                    s(hours_from_slack, bs = "tp") + 
                    s(moon_illuminated, bs = "tp") + 
                    offset(offset),
                  data = catch,
                  family = mgcv::nb(link = "log"),
                  method = "ML")
# functionally equivalent

summary(mod1)
hist(resid(mod1))
mgcv::gam.check(mod1)



# all covariates significant but some worrisome diagnostics; try fitting spatial
# model and check fit


# BUILD MESHES -----------------------------------------------------------------

# construct a few alternative meshes
sdm_mesh1 <- make_mesh(catch,
                      c("xUTM_ds", "yUTM_ds"),
                      n_knots = 250)
sdm_mesh2 <- make_mesh(catch,
                       c("xUTM_ds", "yUTM_ds"),
                       cutoff = 1)


# TOTAL CATCH SPATIAL MODEL ----------------------------------------------------

fit_full <- sdmTMB(
  abund ~ 1 + (1 | year_f) +
    s(year_day, bs = "tp", k = 4) + 
    poly(mean_depth, 2) + 
    mean_slope + 
    poly(moon_illuminated, 2) +
    hours_from_slack , 
  offset = catch$offset,
  data = catch,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  anisotropy = TRUE,
  control = sdmTMBcontrol(
    newton_loops = 1
    # nlminb_loops = 2
  ),
  silent = FALSE
)
sanity(fit_full)
# converges well


## TOTAL CATCH SPATIOTEMPORAL MODELS -------------------------------------------

# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)


# test model
cv_year_ar1 <- sdmTMB_cv(
  abund ~ s(year_day, bs = "tp", k = 4) +
    mean_depth +
    mean_slope +
    poly(moon_illuminated, 2) +
    hours_from_slack, 
  offset = "offset",
  data = catch,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "ar1",
  anisotropy = TRUE,
  share_range = TRUE,
  time = "year",
  control = sdmTMBcontrol(
    newton_loops = 1
  ),
  silent = FALSE,
  use_initial_fit = TRUE,
  k_folds = 5
)

## DOES NTO CONVERGE
# cv_month_iid <- sdmTMB_cv(
#   abund ~ 1 + 
#     (1 | year_f) + 
#     poly(mean_depth, 2) + 
#     mean_slope + 
#     poly(moon_illuminated, 2) +
#     hours_from_slack, 
#   offset = "offset",
#   data = catch,
#   mesh = sdm_mesh1,
#   family = sdmTMB::nbinom2(),
#   spatial = "on",
#   spatiotemporal = "iid",
#   anisotropy = TRUE,
#   share_range = TRUE,
#   time = "month",
#   control = sdmTMBcontrol(
#     newton_loops = 1
#     # nlminb_loops = 2
#   ),
#   silent = FALSE,
#   use_initial_fit = TRUE,
#   k_folds = 5
# )

cv_month_ar1 <- sdmTMB_cv(
  abund ~ (1 | year_f) + 
    mean_depth + 
    mean_slope + 
    poly(moon_illuminated, 2) +
    hours_from_slack, 
  offset = "offset",
  data = catch,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "ar1",
  anisotropy = FALSE,
  share_range = TRUE,
  time = "month",
  control = sdmTMBcontrol(
    # newton_loops = 1
    nlminb_loops = 2
  ),
  silent = FALSE,
  use_initial_fit = TRUE,
  k_folds = 5
)

# does not converge
# cv_week_ar1 <- sdmTMB_cv(
#   abund ~ 1 + 
#     (1 | year_f) + 
#     mean_depth + 
#     mean_slope + 
#     poly(moon_illuminated, 2) +
#     hours_from_slack, 
#   offset = "offset",
#   data = catch,
#   mesh = sdm_mesh1,
#   family = sdmTMB::nbinom2(),
#   spatial = "on",
#   spatiotemporal = "ar1",
#   anisotropy = FALSE,
#   share_range = TRUE,
#   time = "week",
#   control = sdmTMBcontrol(
#     newton_loops = 1
#     # nlminb_loops = 2
#   ),
#   silent = FALSE,
#   use_initial_fit = TRUE,
#   k_folds = 5
# )

# allow yearly distributions to vary spatially (DOES NOT CONVERGE)
# cv_month_ar1B <- sdmTMB_cv(
#   abund ~ 0 + 
#     year_f + 
#     poly(mean_depth, 2) + 
#     hours_from_slack, 
#   offset = "offset",
#   data = catch,
#   mesh = sdm_mesh2,
#   family = sdmTMB::nbinom2(),
#   spatial = "on",
#   spatiotemporal = "ar1",
#   spatial_varying = ~ 0 + year_f,
#   anisotropy = TRUE,
#   share_range = TRUE,
#   time = "month",
#   control = sdmTMBcontrol(
#     # newton_loops = 1
#     nlminb_loops = 2
#   ),
#   silent = FALSE,
#   use_initial_fit = TRUE,
#   k_folds = 5
# )


cv_year_ar1$elpd
cv_month_ar1$elpd


st_year_ar1 <- sdmTMB(
  abund ~ 1 + 
    s(year_day, bs = "tp", k = 4) +
    mean_depth +
    mean_slope +
    poly(moon_illuminated, 2) +
    hours_from_slack,
  offset = "offset",
  data = catch,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "ar1",
  anisotropy = TRUE,
  share_range = TRUE,
  time = "year",
  control = sdmTMBcontrol(
    newton_loops = 1
  ),
  silent = FALSE
)


st_month_ar1 <- sdmTMB(
  abund ~ 0 + 
    year_f + 
    mean_depth +
    mean_slope + 
    poly(moon_illuminated, 2) +
    hours_from_slack, 
  offset = "offset",
  data = catch,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "ar1",
  anisotropy = TRUE,
  share_range = TRUE,
  time = "month",
  control = sdmTMBcontrol(
    # newton_loops = 1
    nlminb_loops = 2
  ),
  silent = FALSE
)


AIC(st_month_ar1, st_year_ar1)
