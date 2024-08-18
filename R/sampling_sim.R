## Unbalanced Sampling Simulation
# Use parameter estimates from catch_sdm_mvrw to generate simple simulation
# evaluating how well biased vs. unbiased samples recover true parameters
# NOTE: excludes temporal FEs and size bins, and spatiotemporal REs to speed up 
# fitting
# Simulate from two datesets, both excluding FEs except size bin, week, and year
# 1) Original dataset 
# 2) Identical but locations sampled at random from pred_grid


library(tidyverse)
library(sf)
library(sdmTMB)

fit_size <- readRDS(here::here("data", "model_fits", "fit_mvrfrw_size.rds"))

summary(fit_size)

total_catch <- fit_size$data %>% 
  filter(!bin == "sublegal") %>% 
  group_by(event) %>% 
  summarize(sum_catch = sum(catch))
dd <- fit_size$data %>% 
  dplyr::select(-c(size_bin, bin, catch)) %>% 
  distinct() %>% 
  left_join(., total_catch, by = "event")
mesh <- make_mesh(dd, xy_cols = c("xUTM_ds", "yUTM_ds"), 
                       n_knots = 175)

fit_size <- sdmTMB(
  sum_catch ~ 1 + (1 | year_f) + 
    (poly(slack_z, 2) * tide_z) +
    depth_z + sunrise_z + slope_z + poly(week_z, 2),
  offset = "offset",
  data = dd,
  mesh = mesh,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  # spatial_varying = ~ 0 + size_bin,
  # time = "month",
  # spatiotemporal = "rw",
  # groups = "bin",
  anisotropy = FALSE,
  # share_range = FALSE,
  silent = FALSE
)


# use predictive spatial grid to define simulation 
pred_grid <- readRDS(here::here("data", "pred_grid_size_model.RDS")) %>% 
  filter(bin == "medium", week == 18.5, depth < 225)
# mesh <- make_mesh(pred_grid, xy_cols = c("xUTM_ds", "yUTM_ds"), 
#                   n_knots = 250)


# original sampling design
orig_locs_sf <- fit_size$data %>% 
  dplyr::select(xUTM, yUTM) %>% 
  distinct() %>% 
  st_as_sf(., coords = c("xUTM", "yUTM"), 
           crs = st_crs("+proj=utm +zone=10 +units=m")) 
# use sf to find nearest neighbors between each
pred_grid_sf <- pred_grid %>% 
  st_as_sf(., coords = c("X", "Y"), 
           crs = st_crs("+proj=utm +zone=10 +units=m"))
nearest_indices <- st_nearest_feature(orig_locs_sf, sim_dat_sf)
orig_grid <- pred_grid_sf[unique(nearest_indices), ] %>% 
  sample_n(size = 500) %>% 
  st_drop_geometry()
orig_mesh <- make_mesh(orig_grid, xy_cols = c("xUTM_ds", "yUTM_ds"),
                       n_knots = 100)

# random design
random_grid <- pred_grid %>% 
  sample_n(500) %>% 
  st_drop_geometry
random_mesh <- make_mesh(random_grid, xy_cols = c("xUTM_ds", "yUTM_ds"),
                          n_knots = 100)


# make empty lists to store par values
niter <- 20
orig_pars <- sim_pars <- vector(mode = "list", length = niter)
orig_pass <- sim_pass <- rep(NA, length = niter)


for(i in seq_len(niter)) {
  # simulate using only spatial relationships
  sim_dat_orig <- sdmTMB_simulate(
    formula = ~ 1 + depth_z + slope_z,
    data = orig_grid,
    mesh = orig_mesh,
    family = nbinom2(),
    range = 25,
    phi = 2,
    sigma_O = 0.7,
    seed = 42 + i,
    B = c(-1, -0.36, 0.1) # intercept, depth and slope effects
  )
  
  ggplot() +
    geom_point(data = sim_dat_orig,
               aes(x = xUTM_ds, y = yUTM_ds, size = observed),
               shape = 21, fill = "#377eb8", alpha = 0.75,
               inherit.aes = FALSE)
  
  fit_orig <- sdmTMB(
    observed ~ 1 + depth_z + slope_z,
    data = sim_dat_orig,
    mesh = orig_mesh,
    family = nbinom2(link = "log")
  )
  
  orig_pass[i] <- ifelse(fit_orig$pos_def_hessian == FALSE, 0, 1) 
  orig_pars[[i]] <- rbind(
    tidy(fit_orig), tidy(fit_orig, effects = "ran_par")
  ) %>% 
    mutate(iteration = i)
  
  
  # random sampling
  sim_dat_random <- sdmTMB_simulate(
    formula = ~ 1 + depth_z + slope_z,
    data = random_grid,
    mesh = random_mesh,
    family = nbinom2(),
    range = 25,
    phi = 2,
    sigma_O = 0.7,
    seed = 42 + i,
    B = c(-1, -0.36, 0.1) # intercept, depth and slope effects
  )
  
  ggplot() +
    geom_point(data = sim_dat_random,
               aes(x = xUTM_ds, y = yUTM_ds, size = observed),
               shape = 21, fill = "red", alpha = 0.75,
               inherit.aes = FALSE)

  fit_sim <- sdmTMB(
    observed ~ 1 + depth_z + slope_z,
    data = sim_dat_random,
    mesh = random_mesh,
    family = nbinom2(link = "log")
  )
  
  sim_pass[i] <- ifelse(fit_sim$pos_def_hessian == FALSE, 0, 1) 
  sim_pars[[i]] <- rbind(
    tidy(fit_sim), tidy(fit_sim, effects = "ran_par")
  ) %>% 
    mutate(iteration = i)
}

sum(sim_pass)
sum(orig_pass)

sim_dat <- sim_pars[sim_pass == 1] %>% 
  bind_rows() %>% 
  mutate(
    sample = "original"
  )

orig_dat <- orig_pars[orig_pass == 1] %>% 
  bind_rows() %>% 
  mutate(
    sample = "random"
  )

ggplot() +
  geom_boxplot(
    data = rbind(sim_dat, orig_dat),
    aes(x = sample, y = estimate)
  ) +
  facet_wrap(~term, scales = "free_y")
