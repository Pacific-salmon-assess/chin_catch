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

dat <- fit_size$data

set_dat <- dat %>% 
  mutate(month = ifelse(month == "9", 8, month)) %>% 
  group_by(event, xUTM_ds, yUTM_ds, depth_z, slope_z, week_z, month) %>% 
  summarize(catch = sum(catch))


mesh_set <- make_mesh(set_dat, xy_cols = c("xUTM_ds", "yUTM_ds"), 
                       n_knots = 100)
mesh_orig <- make_mesh(dat, xy_cols = c("xUTM_ds", "yUTM_ds"), 
                      n_knots = 100)


fit_simp <- sdmTMB(
  catch ~ depth_z + slope_z + poly(week_z, 2),
  offset = "offset",
  data = set_dat,
  mesh = mesh_set,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  anisotropy = FALSE,
  silent = FALSE
)
  
fit_size <- sdmTMB(
  catch ~  bin + depth_z + slope_z + poly(week_z, 2):bin,
  offset = "offset",
  data = dat,
  mesh = mesh_orig,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  time = "month",
  spatiotemporal = "rw",
  groups = "bin",
  anisotropy = TRUE,
  # share_range = FALSE,
  silent = FALSE
)
  

# sample new locations for each event 
pred_bathy_grid <- readRDS(
  here::here("data", "pred_bathy_grid_1000m.RDS")) 
trim_grid <- pred_bathy_grid %>% 
  filter(X < max(dat$xUTM + 1500) & X > min(dat$xUTM - 1500),
         Y < max(dat$yUTM + 1500) &   Y > min(dat$yUTM - 1500))


# subset full grid, by year and week, to represent a) original sampling events
# and b) random sampling events
set.seed(123)
pred_dat <- expand.grid(
  week = unique(dat$week),
  year_f = unique(dat$year_f)
) 
pred_list_orig <- pred_list_rand <- vector(
  mode = "list", length = nrow(pred_dat)
  )
for (i in seq_along(pred_list_orig)) {
  dum_sf <- trim_grid %>% 
    mutate(week = pred_dat$week[i],
           year_f = pred_dat$year_f[i]) %>%
    # add additional variables
    left_join(
      ., 
      dat %>% 
        select(week, week_z, month) %>% 
        distinct(),
      by = "week"
    ) %>% 
    mutate(
      xUTM_ds = X / 1000,
      yUTM_ds = Y / 1000,
      slack_z = 0,
      sunrise_z = 0,
      tide_z = 0,
      month_f = as.factor(month),
      depth_z = (depth - mean(dat$mean_depth)) / sd(dat$mean_depth),
      slope_z = (slope - mean(dat$mean_slope)) / sd(dat$mean_slope)
    ) %>% 
    st_as_sf(., coords = c("X", "Y"), 
             crs = st_crs("+proj=utm +zone=10 +units=m")) 
  orig_dum <- dat %>% 
    filter(week == pred_dat$week[i] & year_f == pred_dat$year_f[i]) %>% 
    select(yUTM, xUTM) %>% 
    distinct() 
  
  
  if (nrow(orig_dum) > 1) {
    orig_dum_sf <- st_as_sf(orig_dum, coords = c("xUTM", "yUTM"), 
                            crs = st_crs("+proj=utm +zone=10 +units=m"))
    nearest_indices <- st_nearest_feature(orig_dum_sf, dum_sf)
    pred_list_orig[[i]]  <- dum_sf[nearest_indices, ]
    pred_list_rand[[i]] <- dum_sf %>% 
      sample_n(nrow(orig_dum_sf))
  }
}

# condense to dataframe and replicate for each size bin
size_bins <- unique(dat$bin)
rand_dd <- orig_dd <- vector(
  mode = "list", length = length(size_bins)
)
for (i in seq_along(size_bins)) {
  rand_dd[[i]] <- pred_list_rand %>% bind_rows() %>% st_drop_geometry() %>%
    mutate(bin = size_bins[i])
  orig_dd[[i]] <- pred_list_orig %>% bind_rows() %>% st_drop_geometry() %>%
    mutate(bin = size_bins[i])
}

# rand_dd <- pred_list_rand %>% bind_rows() %>% st_drop_geometry()
# orig_dd <- pred_list_orig %>% bind_rows() %>% st_drop_geometry()

pred_dat_rand <- rand_dd %>% bind_rows()
pred_dat_orig <- orig_dd %>% bind_rows()

mesh_new <- make_mesh(pred_dat_rand, xy_cols = c("xUTM_ds", "yUTM_ds"), 
                      n_knots = 100)
new_sim <- simulate(fit_simp, nsim = 10, newdata = pred_dat_rand)
pred_dat_rand$sim_catch <- new_sim[,2]

ggplot() +
  geom_point(data = pred_dat_rand %>% filter(!sim_catch == 0),
             aes(x = xUTM_ds, y = yUTM_ds, size = sim_catch),
             shape = 21, fill = "#377eb8", alpha = 0.75,
             inherit.aes = FALSE) +
  facet_wrap(~as.factor(month))
# facet_grid(bin~as.factor(month))

fit_rand <- sdmTMB(
  sim_catch ~ depth_z + slope_z + poly(week_z, 2),
  data = pred_dat_rand,
  mesh = mesh_new,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  anisotropy = FALSE,
  silent = FALSE
)



mesh_new <- make_mesh(pred_dat_orig, xy_cols = c("xUTM_ds", "yUTM_ds"), 
                      n_knots = 100)
new_sim <- simulate(fit_simp, nsim = 10, newdata = pred_dat_orig)
pred_dat_orig$sim_catch <- new_sim[,1]

ggplot() +
  geom_point(data = pred_dat_orig %>% filter(!sim_catch == 0),
             aes(x = xUTM_ds, y = yUTM_ds, size = sim_catch),
             shape = 21, fill = "#377eb8", alpha = 0.75,
             inherit.aes = FALSE) +
  facet_wrap(~as.factor(month))
  # facet_grid(bin~as.factor(month))

fit_new <- sdmTMB(
  sim_catch ~ depth_z + slope_z + poly(week_z, 2),
  data = pred_dat_orig,
  mesh = mesh_new,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  anisotropy = FALSE,
  silent = FALSE
)





new_sim <- simulate(fit_size, nsim = 10, newdata = pred_dat_orig)
pred_dat_orig$sim_catch <- new_sim[,1]
mesh_orig <- make_mesh(pred_dat_orig, xy_cols = c("xUTM_ds", "yUTM_ds"), 
                      n_knots = 100)

fit_size_sim <- sdmTMB(
  sim_catch ~  bin + depth_z + slope_z + poly(week_z, 2):bin,
  # offset = "offset",
  data = pred_dat_orig,
  mesh = mesh_orig,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  time = "month",
  spatiotemporal = "rw",
  groups = "bin",
  anisotropy = TRUE,
  # share_range = FALSE,
  silent = FALSE
)


# take original data, set FEs to zero
pred_dat_orig <- dat %>% 
  select(colnames(pred_dat_new)) #%>% 
  # mutate(sunrise_z = 0, tide_z = 0, slack_z = 0)
mesh_orig <- make_mesh(pred_dat_orig, xy_cols = c("xUTM_ds", "yUTM_ds"), 
                      n_knots = 175)

orig_sim <- simulate(fit_size, nsim = 10, newdata = pred_dat_orig)
pred_dat_orig$sim_catch <- orig_sim[,1]

fit_obs <- sdmTMB(
  sim_catch ~ 1 + (1 | year_f) + bin + 
    depth_z + slope_z + poly(week_z, 2):bin
    ,
  offset = "offset",
  data = pred_dat_orig,
  mesh = mesh_orig,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  # time = "month",
  # spatiotemporal = "rw",
  groups = "bin",
  anisotropy = TRUE,
  # share_range = FALSE,
  silent = FALSE
)


mesh <- make_mesh(dat, xy_cols = c("xUTM_ds", "yUTM_ds"), 
                      n_knots = 175)
fit_new <- sdmTMB(
  catch ~ 1,
  data = dat,
  mesh = mesh,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  # time = "month",
  # spatiotemporal = "rw",
  # groups = "bin",
  anisotropy = FALSE,
  share_range = TRUE,
  silent = FALSE
)





# use predictive spatial grid to define simulation 
pred_grid <- readRDS(here::here("data", "pred_grid_size_model.RDS")) %>% 
  filter(bin == "medium", week == 18.5)
mesh <- make_mesh(pred_grid, xy_cols = c("xUTM_ds", "yUTM_ds"), n_knots = 125)

# use coefs from original model fit to simulate
fit_size <- readRDS(here::here("data", "model_fits", "fit_mvrfrw_size.rds"))


sum_dat <- dat %>% 
  group_by(event) %>% 
  summarize(catch = sum(catch)) %>% 
  left_join(., dat %>% select(event, year_f:sunrise_z) %>% distinct(), by  = "event")
sum_mesh <- 
fit2 <- sdmTMB(
  catch ~ 1,
  data = dat,
  mesh = mesh,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  # time = "month",
  # spatiotemporal = "rw",
  # groups = "bin",
  anisotropy = FALSE,
  share_range = TRUE,
  silent = FALSE
)


fes <- tidy(fit_size, effects = "fixed") %>% 
  filter(
    term %in% c("binmedium", "depth_z", "slope_z")
  )
res <- tidy(fit_size, effects = "ran_pars")
  
# simulate using only spatial relationships
sim_dat <- sdmTMB_simulate(
  formula = ~ 1 + depth_z + slope_z,
  data = pred_grid,
  mesh = mesh,
  family = nbinom1(),
  range = res$estimate[[1]],
  # sigma_E = 0.1,
  phi = res$estimate[[3]],
  sigma_O = res$estimate[[4]],
  seed = 42,
  B = fes$estimate # intercept, depth and slope effects
)



# define original sampling frame
orig_locs_sf <- fit_size$data %>% 
  dplyr::select(xUTM, yUTM) %>% 
  distinct() %>% 
  st_as_sf(., coords = c("xUTM", "yUTM"), 
           crs = st_crs("+proj=utm +zone=10 +units=m")) 

ggplot() + 
  geom_point(data = sim_dat %>% filter(!observed == 0),
             aes(x = xUTM_ds, y = yUTM_ds, size = observed),
             shape = 21, fill = "#377eb8", alpha = 0.75,
             inherit.aes = FALSE) 


# use sf to find nearest neighbors between each
sim_dat_sf <- sim_dat %>% 
  mutate(xUTM = xUTM_ds * 1000, 
         yUTM = yUTM_ds * 1000) %>% 
  st_as_sf(., coords = c("xUTM", "yUTM"), 
           crs = st_crs("+proj=utm +zone=10 +units=m"))
nearest_indices <- st_nearest_feature(orig_locs_sf, sim_dat_sf)

sim_dat_orig <- sim_dat_sf[nearest_indices, ] %>% 
  st_drop_geometry()
mesh_orig <- make_mesh(sim_dat_orig, xy_cols = c("xUTM_ds", "yUTM_ds"),
                       n_knots = 175)
fit_orig <- sdmTMB(
  observed ~ 1 + depth_z + slope_z,
  data = sim_dat_orig,
  mesh = mesh_orig,
  family = nbinom1(link = "log")
)

sim_dat_new <- sim_dat_sf %>% 
  sample_n(size = nrow(sim_dat_orig), replace = TRUE) %>% 
  st_drop_geometry()
mesh_sim <- make_mesh(sim_dat_new, xy_cols = c("xUTM_ds", "yUTM_ds"),
                       n_knots = 175)

fit_orig <- sdmTMB(
  observed ~ 1 + depth_z + slope_z,
  data = sim_dat_new,
  mesh = mesh_orig,
  family = nbinom1()
)
