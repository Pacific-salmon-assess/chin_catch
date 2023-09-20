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
  ) %>% 
  droplevels()


# BUILD MESHES -----------------------------------------------------------------

# construct a few alternative meshes
sdm_mesh1 <- make_mesh(catch_stock,
                       c("xUTM_ds", "yUTM_ds"),
                       n_knots = 125)


# SPATIAL PREDICTION GRID ------------------------------------------------------

# import 1 x 1 km grid based on shrunk CH shape files and generated in 
# prep_prediction_grid.R 
pred_bathy_grid <- readRDS(
  here::here("data", "pred_bathy_grid_1000m.RDS"))

# generate combinations for month/year
pred_dat <- expand.grid(
  week = seq(min(catch_size$week), max(catch_size$week), length = 5),
  year_f = unique(catch_size$year_f),
  size_bin = unique(catch_size$size_bin)
) %>% 
  filter(year_f == "2021") 

pred_grid_list <- vector(mode = "list", length = nrow(pred_dat))
for (i in seq_along(pred_grid_list)) {
  pred_grid_list[[i]] <- pred_bathy_grid %>% 
    filter(X < max(catch_size$xUTM + 1000) & X > min(catch_size$xUTM - 1000),
           Y < max(catch_size$yUTM + 1000) & Y > min(catch_size$yUTM - 1000)) %>% 
    mutate(week = pred_dat$week[i],
           year_f = pred_dat$year_f[i],
           size_bin = pred_dat$size_bin[i])
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


## FIT PRELIM MODELS -----------------------------------------------------------


# fit MVRW version w/ top four stocks
fit_full <- sdmTMB(
  catch ~ 0 + (1 | year_f) + agg + depth_z + slope_z + week_z,
  offset = "offset",
  data = catch_stock,
  mesh = sdm_mesh1,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  # spatial_varying = ~ 0 + size_bin,
  time = "month",
  spatiotemporal = "rw",
  groups = "agg",
  anisotropy = TRUE,
  share_range = TRUE,
  silent = FALSE
)
fit_full2 <- update(fit_full, share_range = FALSE)


# fit individual stock models

stock_tbl <- catch_stock %>% 
  group_by(agg) %>% 
  group_nest()

fit_list <- purrr::map(
  stock_tbl$data,
  
  ~ sdmTMB(
    catch ~ week_z + (1 | year_f),
    offset = "offset",
    data = .x,
    mesh = sdm_mesh1,
    family = sdmTMB::nbinom1(),
    spatial = "on",
    time = "month",
    spatiotemporal = "rw",
    anisotropy = TRUE,
    silent = FALSE
  )
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
f3b <- sdmTMB(
  catch ~ week_z + depth_z + slope_z + (1 | year_f),
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
      c(-2.5, 0.25, -0.25, 0.25),
      #sd
      c(2, 1, 1, 1)
    ),
    matern_s = pc_matern(range_gt = 40, sigma_lt = 10)
  ),
  silent = FALSE
)

# fails to converge
# f4 <- sdmTMB(
#   catch ~ week_z + (1 | year_f),
#   offset = "offset",
#   data = fr_sub,
#   mesh = sdm_mesh1,
#   family = sdmTMB::nbinom1(),
#   spatial = "off",
#   spatial_varying = ~ 0 + month_f,
#   anisotropy = FALSE,
#   priors = sdmTMBpriors(
#     b = normal(
#       #mean
#       c(-2.5, 0.25),
#       #sd
#       c(2, 1)
#     ),
#     matern_s = pc_matern(range_gt = 40, sigma_lt = 10)
#   ),
#   control = sdmTMBcontrol(
#     newton_loops = 1,
#     start = list(
#       ln_phi = log(0.16) #approximate value estimated in catch_sdm
#     ),
#     map = list(
#       ln_phi = factor(NA),
#       ln_tau_Z = factor(
#         seq(1, length(unique(catch_stock$month_f)), by = 1)
#       )
#     )
#   ),
#   silent = FALSE
# )

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
    b = normal(
      #mean
      c(-2.5, 0.25),
      #sd
      c(2, 1)
    ),
    matern_s = pc_matern(range_gt = 40, sigma_lt = 10),
    matern_st = pc_matern(range_gt = 40, sigma_lt = 10)
  ),
  control = sdmTMBcontrol(
    newton_loops = 1,
    start = list(
      ln_phi = log(0.16) #approximate value estimated in catch_sdm
    ),
    map = list(
      ln_phi = factor(NA)
      )
  ),
  silent = FALSE
)
# fails to converge
# f5b <- sdmTMB(
#   catch ~ week_z + depth_z + slope_z + (1 | year_f),
#   offset = "offset",
#   data = fr_sub,
#   mesh = sdm_mesh1,
#   family = sdmTMB::nbinom1(),
#   spatial = "on",
#   time = "month",
#   spatiotemporal = "iid",
#   anisotropy = FALSE,
#   priors = sdmTMBpriors(
#     b = normal(
#       #mean
#       c(-2.5, 0.25, -0.25, 0.25),
#       #sd
#       c(2, 1, 1, 1)
#     ),
#     matern_s = pc_matern(range_gt = 40, sigma_lt = 10),
#     matern_st = pc_matern(range_gt = 40, sigma_lt = 10)
#   ),
#   control = sdmTMBcontrol(
#     newton_loops = 1,
#     start = list(
#       ln_phi = log(0.16) #approximate value estimated in catch_sdm
#     ),
#     map = list(
#       ln_phi = factor(NA)
#     )
#   ),
#   silent = FALSE
# )


# check most complicated model
mod_in <- f5
samp <- sample_mle_mcmc(mod_in, mcmc_iter = 130, mcmc_warmup = 100)

obj <- mod_in$tmb_obj
random <- unique(names(obj$env$par[obj$env$random]))
pl <- as.list(mod_in$sd_report, "Estimate")
fixed <- !(names(pl) %in% random)
map <- lapply(pl[fixed], function(x) factor(rep(NA, length(x))))
obj <- TMB::MakeADFun(obj$env$data, pl, map = map, DLL = "sdmTMB")
obj_mle <- mod_in
obj_mle$tmb_obj <- obj
obj_mle$tmb_map <- map
ss <- simulate(obj_mle, mcmc_samples = extract_mcmc(samp), nsim = 30)


pred_fixed <- mod_in$family$linkinv(
  predict(mod_in, newdata = fr_sub)$est_non_rf
)
dharma_res_out <- DHARMa::createDHARMa(
  simulatedResponse = ss,
  observedResponse = mod_in$data$catch,
  fittedPredictedResponse = pred_fixed
)
DHARMa::testResiduals(dharma_res_out)



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
    month = case_when(
      week <= 22 ~ 5,
      week > 22 & week <= 26 ~ 6,
      week > 26 & week <= 30 ~ 7,
      week > 30 ~ 8
    ),
    month_f = as.factor(month),
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


## SPATIAL PREDICTIONS ---------------------------------------------------------

p3 <- predict(f3, newdata = pred_grid, se_fit = FALSE, re_form = NULL) %>% 
  mutate(model = "f3",
         epsilon_st = NA)
p5 <- predict(f5, newdata = pred_grid, se_fit = FALSE, re_form = NULL) %>% 
  mutate(model = "f5")
pp <- rbind(p3, p5) %>% 
  mutate(exp_est = exp(est))


plot_map(pp, omega_s) +
  scale_fill_gradient2(name = "Spatial\nRF Effect") +
  facet_wrap(~model)

plot_map(p5, epsilon_st) +
  scale_fill_gradient2(name = "Spatiotemporal\nRF Effect") +
  facet_wrap(~month)

plot_map(pp, exp_est) +
  facet_grid(model~month) +
  scale_fill_viridis_c()

month_omega <- pp %>%
  dplyr::select(-year_f) %>% 
  distinct() %>% 
  pivot_longer(cols = starts_with("zeta_s_year"),
               names_prefix = "zeta_s_year_f",
               names_to = "year") %>% 
  plot_map(., value) +
  scale_fill_gradient2(name = "Spatial\nRF Effect") +
  facet_wrap(~year) +
  theme(legend.position = "top")

full_preds <- pp %>% 
  group_by(size_bin) %>% 
  mutate(scale_est = exp(est) / max(exp(est))) %>%
  left_join(., week_key, by = "week") %>% 
  plot_map(., scale_est) +
  scale_fill_viridis_c(trans = "sqrt", name = "Scaled\nAbundance") +
  facet_grid(date~size_bin)


png(here::here("figs", "ms_figs", "month_omega.png"), res = 250, units = "in", 
    height = 7.5, width = 7.5)
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