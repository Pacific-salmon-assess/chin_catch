## Fit preliminary spatial distribution models
# August 7, 2020
# updated Feb 22, 2023
# Ultimately estimate stock- or size-specific spatial distribution maps
# Given cross-validation (catch_sdm_cv.R) and goals of analysis, utilize
# month focused spatio temporal model
# Also include static spatial variables (e.g. depth/slope) and dynamic 
# fine-scale (e.g. time of day) 

library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(raster)
library(sf)
library(sdmTMBextra)

# Import stage data and fitted model generated in gen_detection_histories.R
stage_dat <- readRDS(here::here("data", "agg_lifestage_df.RDS")) %>% 
  dplyr::select(vemco_code, stage)


chin_raw <- readRDS(here::here("data", "cleanTagData_GSI.RDS")) %>%
  mutate(year_day = lubridate::yday(date)) %>% 
  rename(vemco_code = acoustic_year, agg = agg_name) %>% 
  left_join(., stage_dat, by = "vemco_code") 


## predict life stage based on fitted model 
stage_mod <- readRDS(here::here("data", "stage_fl_hierA.RDS"))
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
    month = lubridate::month(date),
    size_bin = case_when(
      fl < 65 ~ "small",
      fl >= 65 & fl < 74 ~ "medium",
      fl >= 74 ~ "large"
    )
  ) %>% 
  dplyr::select(fish, vemco_code, event, date, year, month, year_day,
                fl, size_bin, clip, stage, stock:agg_prob) 

chin %>% 
  group_by(size_bin) %>% 
  tally()
catch_stock %>% 
  group_by(agg) %>% 
  tally()
    
# calculate total catch across size bins
catch_size <- chin %>% 
  group_by(event, size_bin) %>% 
  summarize(catch = n(), .groups = "drop") %>% 
  ungroup()
catch_stock <- chin %>%
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


# import lighthouse data
# sst <- readRDS(here::here("data", "spatial_data", "trim_amphitrite.rds"))


catch <- expand.grid(event = set_dat$event,
                     size_bin = unique(catch_size$size_bin)) %>% 
  arrange(event) %>% 
  left_join(., catch_size, by = c("event", "size_bin")) %>%
  replace_na(., replace = list(catch = 0)) %>% 
  left_join(., set_dat, by = "event") %>%
  mutate(
    year_f = as.factor(year),
    yUTM_ds = yUTM / 1000,
    xUTM_ds = xUTM / 1000,
    offset = log(set_dist),
    depth_z = scale(mean_depth)[ , 1],
    slack_z = scale(hours_from_slack)[ , 1],
    week_z = scale(week)[ , 1]
  ) %>% 
  # remove sets not on a troller
  filter(!grepl("rec", event))

catch_stock <- expand.grid(
  event = set_dat$event,
  agg = unique(catch_stock$agg)
) %>% 
  arrange(event) %>% 
  left_join(., catch_stock, by = c("event", "agg")) %>%
  replace_na(., replace = list(catch = 0)) %>% 
  left_join(., set_dat, by = "event") %>%
  mutate(
    year_f = as.factor(year),
    yUTM_ds = yUTM / 1000,
    xUTM_ds = xUTM / 1000,
    offset = log(set_dist),
    depth_z = scale(mean_depth)[ , 1],
    slack_z = scale(hours_from_slack)[ , 1],
    week_z = scale(week)[ , 1]
  ) %>% 
  # remove sets not on a troller
  filter(!grepl("rec", event))


catch_tbl_size <- as_tibble(catch) %>% 
  group_by(size_bin) %>% 
  group_nest() %>% 
  rename(group_var = size_bin) %>% 
  mutate(grouping = "size")
catch_tbl_stock <- as_tibble(catch_stock) %>% 
  group_by(agg) %>% 
  group_nest() %>%
  rename(group_var = agg) %>% 
  mutate(grouping = "stock")
catch_tbl <- rbind(catch_tbl_size, catch_tbl_stock)


# BUILD MESHES -----------------------------------------------------------------

# construct a few alternative meshes ( can be applied to all size bins so only 
# one needed)
sdm_mesh1 <- make_mesh(catch %>% filter(size_bin == "large"),
                      c("xUTM_ds", "yUTM_ds"),
                      n_knots = 250)
sdm_mesh2 <- make_mesh(catch %>% filter(size_bin == "large"),
                       c("xUTM_ds", "yUTM_ds"),
                       n_knots = 150)


# SPATIAL MODELS ---------------------------------------------------------------

ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)

sub_tbl <- catch_tbl %>% filter(grouping == "stock")

sp_fits <- furrr::future_map(
  sub_tbl$data,
  ~ sdmTMB(
    catch ~ (1 | year_f) +
      mean_depth +
      hours_from_slack
      , 
    offset = "offset",
    data = .x,
    mesh = sdm_mesh1,
    family = sdmTMB::nbinom2(),
    spatial = "on",
    anisotropy = FALSE,
    control = sdmTMBcontrol(
      newton_loops = 1
      # nlminb_loops = 2
    ),
    silent = FALSE
  )
)

purrr::map(sp_fits, sanity)
# converges well

# residual plots
purrr::map(
  sp_fits, 
  ~ {
   sims <- simulate(.x, nsim = 30)
   dharma_residuals(sims, .x)
  }
)
# look good


# spatial residuals 
sp_plots <- purrr::map2(
  catch_tbl$data,
  catch_tbl$sp_fits, 
  ~ {
    .x$resids <- residuals(.y)
    ggplot(data = .x) +
      geom_point(aes(x = xUTM, y = yUTM, colour = resids)) +
      scale_color_gradient2() +
      facet_wrap(~year) +
      ggsidekick::theme_sleek()
  }
)
# also look good


## TOTAL CATCH SPATIOTEMPORAL MODELS -------------------------------------------

catch_tbl$st_fits[4:6] <- furrr::future_map(
  catch_tbl$data[4:6],
  ~ sdmTMB(
    catch ~ #0 + year_f 
      (1 | year_f)# +
      # mean_depth #+
      # hours_from_slack
    , 
    offset = "offset",
    data = .x,
    mesh = sdm_mesh2,
    family = sdmTMB::nbinom2(),
    spatial = "on",
    spatiotemporal = "iid",
    anisotropy = FALSE,
    time = "month",
    control = sdmTMBcontrol(
      newton_loops = 1
      # nlminb_loops = 2
    ),
    silent = FALSE
  )
)

purrr::map(catch_tbl$st_fits[4:6], sanity)
# converges well


## RESIDUAL CHECKS -------------------------------------------------------------

# residual plots
purrr::map(
  catch_tbl$st_fits, 
  ~ {
    sims <- simulate(.x, nsim = 30)
    dharma_residuals(sims, .x)
  }
)
# look good


# spatial residuals 
sp_plots <- purrr::map2(
  catch_tbl$data,
  catch_tbl$st_fits, 
  ~ {
    .x$resids <- residuals(.y)
    ggplot(data = .x) +
      geom_point(aes(x = xUTM, y = yUTM, colour = resids)) +
      scale_color_gradient2() +
      facet_wrap(~year) +
      ggsidekick::theme_sleek()
  }
)
# also look good

purrr::map(
  catch_tbl$st_fits, 
  ~ {
    sims <- simulate(.x, nsim = 30)
    breaks_vec <- c(seq(0, 16, by = 1))
    hist(sims[ , 1], col="green", pch=20, cex=4, 
         breaks = breaks_vec)
    hist(catch$catch, pch=20, cex=4, col=rgb(1,0,0,0.5), add=TRUE, 
         breaks = breaks_vec)
    }
)



## FIXED EFFECTS ---------------------------------------------------------------

for (i in 1:nrow(catch_tbl)) {
  file_name <- paste(catch_tbl$size_bin[[i]], "st_fixed_effects.pdf", sep = "_")
  pdf(here::here("figs", file_name))
  visreg::visreg(catch_tbl$st_fits[[i]], 
                 xvar = "hours_from_slack", scale = "response")
  visreg::visreg(catch_tbl$st_fits[[i]], 
                 xvar = "mean_depth", scale = "response")
  dev.off()
}
  

# PREDICTION GRID --------------------------------------------------------------

# import 1.8 X 1.8 km grid generated in prep_bathymetry.R in chinTagging repo
pred_bathy_grid <- readRDS(
  here::here("data", "pred_bathy_grid_utm.RDS"))

# generate combinations for month/year
month_year <- expand.grid(
  month = c(5, 6, 7, 8, 9),
  year_day = c(166, 196, 227, 258),
  year = unique(catch$year)
)

pred_grid_list <- vector(mode = "list", length = nrow(month_year))
for (i in seq_along(pred_grid_list)) {
  pred_grid_list[[i]] <- pred_bathy_grid %>% 
    filter(X < max(catch$xUTM + 1000) & X > min(catch$xUTM - 1000),
           Y < max(catch$yUTM + 1000) & Y > min(catch$yUTM - 1000)) %>% 
    mutate(year = month_year$year[i],
           month = month_year$month[i],
           year_day = month_year$year_day[i])
}

pred_grid <- pred_grid_list %>% 
  bind_rows() %>% 
  mutate(
    xUTM_ds = X / 1000,
    yUTM_ds = Y / 1000,
    coast_dist_km = shore_dist / 1000,
    hours_from_slack = median(catch$hours_from_slack),
    moon_illuminated = 0.5,
    year_f = as.factor(year)
  ) %>% 
  dplyr::select(
    X, Y, xUTM_ds, yUTM_ds, month, year, year_f, year_day, hours_from_slack, 
    moon_illuminated, mean_depth = depth, mean_slope = slope
  )


## SPATIAL PREDICTIONS ---------------------------------------------------------

catch_tbl$sim_preds <- purrr::map(
  catch_tbl$st_fits,
  ~ predict(.x, newdata = pred_grid, nsim = 20, type = "response")
)
catch_tbl$preds <- purrr::map(
  catch_tbl$st_fits,
  ~ predict(.x, newdata = pred_grid)
)

catch_preds <- catch_tbl %>% 
  dplyr::select(size_bin, preds) %>% 
  unnest(cols = c(preds))


# helper function for simple predictive maps
plot_map <- function(dat, column) {
  ggplot() +
    geom_raster(data = dat, aes(X, Y, fill = {{ column }})) +
    coord_fixed() +
    ggsidekick::theme_sleek()
}


#Predictions incorporating all fixed and random effects
# Generally the impacts of different meshes considered here were negligible
month_preds <- plot_map(catch_preds %>% filter(year == "2019"), exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)") +
  facet_grid(size_bin~month) 
year_preds <- plot_map(catch_preds %>% filter(month == "7"), exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)") +
  facet_grid(size_bin~year) #+


# Spatial random effects (i.e. independent of time)
omega <- plot_map(catch_preds, omega_s) +
  scale_fill_gradient2() +
  ggtitle("Spatial Random Effects") +
  facet_wrap(~size_bin) #+

# Spatiotemporal random effects
eps <- plot_map(catch_preds, epsilon_st) +
  scale_fill_gradient2() +
  ggtitle("Spatiotemporal Random Effects") +
  facet_grid(size_bin~month) 


pdf(here::here("figs", "spatial_preds_sizebins.pdf"))
month_preds
year_preds
omega
eps
dev.off()