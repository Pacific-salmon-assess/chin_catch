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

# original dataset, filtered to include only medium size class
set_dat <- dat %>% 
  select(event, year_f, week_z, month) %>% 
  distinct()


# sample new locations for each event 
pred_grid1 <- readRDS(here::here("data", "pred_grid_size_model.RDS"))
pred_grid_sub <- pred_grid1 %>% 
  filter(week == 18.5, bin == "medium") %>%
  sample_n(., size = nrow(set_dat), replace = TRUE) %>% 
  select(xUTM_ds, yUTM_ds, depth_z, slope_z) 


pred_dat_new <- purrr::map(
  unique(dat$bin), function (x) {
    cbind(set_dat, pred_grid_sub) %>% 
      mutate(bin = x,
             sunrise_z = 0,
             tide_z = 0,
             slack_z = 0)
  }
) %>% 
  bind_rows()
mesh_new <- make_mesh(pred_dat_new, xy_cols = c("xUTM_ds", "yUTM_ds"), 
                      n_knots = 175)

new_sim <- simulate(fit_size, nsim = 10, newdata = pred_dat_new, 
                    params = "mvn")
pred_dat_new$sim_catch <- new_sim[,1]


library(sdmTMBextra)
samp <- sample_mle_mcmc(
  fit_size, mcmc_iter = 110, mcmc_warmup = 100, 
  mcmc_chains = 4,
  stan_args = list(thin = 5, cores = 4)
)
obj <- fit_size$tmb_obj
random <- unique(names(obj$env$par[obj$env$random]))
pl <- as.list(fit_size$sd_report, "Estimate")
fixed <- !(names(pl) %in% random)
map <- lapply(pl[fixed], function(x) factor(rep(NA, length(x))))
obj <- TMB::MakeADFun(obj$env$data, pl, map = map, DLL = "sdmTMB")
obj_mle <- fit_size
obj_mle$tmb_obj <- obj
obj_mle$tmb_map <- map
sim_out <- simulate(
  obj_mle, mcmc_samples = sdmTMBextra::extract_mcmc(samp), nsim = 10
)



ggplot() + 
  geom_point(data = pred_dat_orig %>% filter(!sim_catch == 0),
             aes(x = xUTM_ds, y = yUTM_ds, size = sim_catch),
             shape = 21, fill = "#377eb8", alpha = 0.75,
             inherit.aes = FALSE) +
  facet_grid(bin~as.factor(month)) 

fit_new <- sdmTMB(
  sim_catch ~ 1,
  data = pred_dat_new,
  mesh = mesh_new,
  family = sdmTMB::nbinom1(),
  spatial = "on",
  # time = "month",
  # spatiotemporal = "rw",
  # groups = "bin",
  anisotropy = FALSE,
  share_range = TRUE,
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
  filter(bin == "medium")
mesh <- make_mesh(pred_grid, xy_cols = c("xUTM_ds", "yUTM_ds"), n_knots = 175)

# use coefs from original model fit to simulate
fit_size <- readRDS(here::here("data", "model_fits", "fit_mvrfrw_size.rds"))

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
  family = nbinom1(link = "log"),
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
