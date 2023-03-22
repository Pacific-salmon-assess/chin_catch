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
    month = lubridate::month(date),
    size_bin = case_when(
      fl < 65 ~ "small",
      fl >= 65 & fl < 75 ~ "medium",
      fl >= 75 ~ "large"
    )
  ) %>% 
  dplyr::select(fish, vemco_code, event, date, year, month, year_day,
                fl, size_bin, clip, stage, stock:agg_prob) 

chin %>% 
  group_by(size_bin) %>% 
  tally()
    
# calculate total catch across size bins
catch_size <- chin %>% 
  group_by(event, size_bin) %>% 
  summarize(catch = n(), .groups = "drop") %>% 
  ungroup()
catch_stock <- chin %>%
  filter(!is.na(agg)) %>% 
  group_by(event, agg) %>% 
  summarize(catch = n(), .groups = "drop") %>% 
  ungroup()


# clean and bind to set data
set_dat <- readRDS(here::here("data", "tagging_data", "cleanSetData.RDS")) %>% 
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
  # left_join(., sst, by = c("month", "year")) %>% 
  mutate(
    year_f = as.factor(year),
    yUTM_ds = yUTM / 1000,
    xUTM_ds = xUTM / 1000,
    offset = log(set_dist),
    # sst_z = scale(sst)[ , 1],
    week_z = scale(week)[ , 1]
  ) %>% 
  # remove sets not on a troller
  filter(!grepl("rec", event))


catch_tbl <- as_tibble(catch) %>% 
  group_by(size_bin) %>% 
  group_nest()


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

catch_tbl$sp_fits <- furrr::future_map(
  catch_tbl$data,
  ~ sdmTMB(
    catch ~ #0 + year_f 
      (1 | year_f) +
      mean_depth +
      mean_slope +
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

purrr::map(catch_tbl$sp_fits, sanity)
# converges well

# residual plots
purrr::map(
  catch_tbl$sp_fits, 
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

catch_tbl$st_fits <- furrr::future_map(
  catch_tbl$data,
  ~ sdmTMB(
    catch ~ #0 + year_f 
      (1 | year_f) +
      mean_depth +
      # mean_slope +
      hours_from_slack
    , 
    offset = "offset",
    data = .x,
    mesh = sdm_mesh2,
    family = sdmTMB::nbinom2(),
    spatial = "on",
    spatiotemporal = "iid",
    anisotropy = FALSE,
    share_range = TRUE,
    time = "month",
    control = sdmTMBcontrol(
      newton_loops = 1
      # nlminb_loops = 2
    ),
    silent = FALSE
  )
)

purrr::map(catch_tbl$st_fits, sanity)
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


# residuals look good 
catch$resids <- residuals(st_month_ar1)
ggplot(data = catch) +
  geom_point(aes(x = xUTM, y = yUTM, colour = resids)) +
  scale_color_gradient2() +
  facet_wrap(~year) +
  ggsidekick::theme_sleek()


## FIXED EFFECTS ---------------------------------------------------------------

for (i in 1:nrow(catch_tbl)) {
  file_name <- paste(catch_tbl$size_bin[[i]], "st_fixed_effects.pdf", sep = "_")
  pdf(here::here("figs", "sdm_catch", file_name))
  visreg::visreg(catch_tbl$st_fits[[i]], 
                 xvar = "hours_from_slack", scale = "response")
  visreg::visreg(catch_tbl$st_fits[[i]], 
                 xvar = "mean_depth", scale = "response")
  dev.off()
}
  

# PREDICTION GRID --------------------------------------------------------------

# import 1.8 X 1.8 km grid generated in prep_bathymetry.R
pred_bathy_grid <- readRDS(
  here::here("data", "spatial_data", "pred_bathy_grid_utm.RDS"))

# generate combinations for month/year
month_year <- expand.grid(
  month = c(6, 7, 8, 9),
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

preds_year <- predict(st_year_ar1, newdata = pred_grid)
preds_month <- predict(st_month_ar1, newdata = pred_grid)

preds_both <- rbind(
  preds_year %>% 
    mutate(
      model = "year_st"
    ),
  preds_month %>% 
    mutate(
      model = "month_st"
    )
  )

# helper function for simple predictive maps
plot_map <- function(dat, column) {
  ggplot() +
    geom_raster(data = dat, aes(X, Y, fill = {{ column }})) +
    coord_fixed() +
    ggsidekick::theme_sleek()
}

#Predictions incorporating all fixed and random effects
# Generally the impacts of different meshes considered here were negligible
plot_map(preds_both %>% filter(year == "2020"), exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)") +
  facet_grid(model~month) 
plot_map(preds_both %>% filter(month == "7"), exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)") +
  facet_grid(model~year) #+



# Spatial random effects (i.e. independent of time)
plot_map(preds_both, omega_s) +
  scale_fill_gradient2() +
  ggtitle("Prediction (fixed effects + all random effects)") +
  facet_grid(model~year) #+

# Spatiotemporal random effects
plot_map(preds_both, epsilon_st) +
  scale_fill_gradient2() +
  ggtitle("Prediction (fixed effects + all random effects)") +
  facet_grid(model~month) 
plot_map(preds_both, epsilon_st) +
  scale_fill_gradient2() +
  ggtitle("Prediction (fixed effects + all random effects)") +
  facet_grid(model~year) 



# plot against map
coast <- rbind(rnaturalearth::ne_states( "United States of America",
                                         returnclass = "sf"),
               rnaturalearth::ne_states( "Canada", returnclass = "sf"))
box <- c(xmin = -128.5, ymin = 48, xmax = -122, ymax = 51)
coast_trim <- st_crop(coast, box)
coast_trimUTM <- st_transform(coast_trim, 
                              crs = "+proj=utm +zone=10 +datum=WGS84")

plot_map(preds3$data, "exp(est)") +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)") +
  facet_wrap(~year_f) +
  geom_sf(data = coast_trimUTM)

plot_map(preds3$data, "omega_s") +
  scale_fill_viridis_c() +
  ggtitle("Spatial random effects only")+
  geom_point(data = catch_trim, aes(x = xUTM, y = yUTM), colour = "red")


## PRESENTATION MAP ------------------------------------------------------------

#inset map
minLat3 <- min(catch$lat)
maxLat3 <- max(catch$lat)
minLong3 <- min(catch$long)
maxLong3 <- max(catch$long)
inset_map <- readRDS(here::here("data", "generated_data", "base_map.RDS")) +
  geom_rect(aes(xmin = minLong3, xmax = maxLong3, ymin = minLat3, 
                ymax = maxLat3),
            fill = NA, lty = 2, col = "black") +
  theme(strip.background = element_rect(colour="white", fill="white"),
        legend.position = "top",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  coord_fixed(xlim = c(-128.5, -122), ylim = c(48, 51), ratio = 1.3)


#primary map
min_utm_y <- min(preds3$data$Y) 
max_utm_y <- max(preds3$data$Y) 
min_utm_x <- min(preds3$data$X) 
max_utm_x <- max(preds3$data$X) 
spatial_effects <- ggplot() +
  geom_raster(data = preds3$data, aes_string("X", "Y", fill = "omega_s")) +
  geom_sf(data = coast_trimUTM) +
  scale_fill_viridis_c(name = "Chinook\nDensity\nAnomaly") +
  coord_sf(ylim = c(min_utm_y, max_utm_y), xlim = c(min_utm_x, max_utm_x)) +
  ggsidekick::theme_sleek() + 
  theme(axis.title = element_blank())

png(here::here("figs", "sdm", "rand_effs_inset.png"), height = 5, width = 5,
    units = "in", res = 250)
cowplot::ggdraw() +
  cowplot::draw_plot(spatial_effects) +
  cowplot::draw_plot(inset_map, x = 0.1, y = 0.16, width = 0.3, height = 0.3)
dev.off()


## NON-SPATIAL PREDICTIONS -----------------------------------------------------

# starting dataset
blank_new_dat <- expand.grid(
  year = unique(catch_trim$year),
  year_day_z = 0,
  hour_z = 0,
  hr_slack_z = 0,
  mean_depth_z = 0, 
  sd_depth_z = 0,
  year_f = levels(catch_trim$year_f),
  offset = mean(catch_trim$offset)
)

marg_plot <- function(model, exp_var) {
  dum <- blank_new_dat %>% 
    dplyr::select(!{{ exp_var }}) %>% 
    expand_grid(., 
                xx = seq(min(catch_trim[ , exp_var]), 
                         max(catch_trim[ , exp_var]),
                         length.out = 100)) 
  dum[ , exp_var] <- dum$xx
  
  p <- predict(model, newdata = dum, se_fit = TRUE, re_form = NA)
  
  ggplot(p,
         aes(xx, exp(est),
                    ymin = exp(est - 1.96 * est_se),
                    ymax = exp(est + 1.96 * est_se),
                    colour = year_f)) +
    coord_cartesian(ylim = c(0, 3.5)) +
    geom_line(aes(colour = year_f)) +
    geom_ribbon(aes(fill = year_f), alpha = 0.4)
}

depth_marg <- marg_plot(model = mod3, exp_var = "mean_depth_z")
sd_depth_marg <- marg_plot(model = mod3, exp_var = "sd_depth_z")
year_day_marg <- marg_plot(model = mod3, exp_var = "year_day_z")
hour_marg <- marg_plot(model = mod3, exp_var = "hour_z")
hr_slack_marg <- marg_plot(model = mod3, exp_var = "hr_slack_z")

