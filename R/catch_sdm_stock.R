## Fit stock-based spatial distribution models
# July 5, 2023
# As in catch_sdm_mvrw, fits different groups using MVRW

library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(raster)
library(sf)
library(sdmTMBextra)


chin <- readRDS(here::here("data", "clean_catch.RDS")) %>% 
  rename(agg = agg_name) %>% 
  filter(
    !is.na(agg),
    fl > 65
  )


chin %>% 
  group_by(agg) %>% 
  tally()


# calculate total catch across size bins
catch_stock1 <- chin %>%
  filter(# remove rare stocks
         !agg %in% c("ECVI", "Fraser Year.")) %>%
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
    effort = (set_dist / 1000) * lines,
    offset = log(effort),
    # effort of one corresponds to ~2 lines towed for 200 m
    depth_z = scale(mean_depth)[ , 1],
    slope_z = scale(mean_slope)[ , 1],
    slack_z = scale(hours_from_slack)[ , 1],
    week_z = scale(week)[ , 1],
    moon_z = scale(moon_illuminated)[ , 1],
    dist_z = scale(coast_dist_km)[ , 1]
  ) %>% 
  droplevels()


# check average cpue by week/size bin
catch_stock %>% 
  group_by(week, agg, year) %>% 
  summarize(
    sum_catch = sum(catch),
    sum_effort = sum(effort)
  ) %>% 
  mutate(
    cpue = sum_catch / sum_effort
  ) %>% 
  ggplot(.) +
  geom_point(aes(x = week, y = log(cpue))) +
  facet_wrap(~agg)



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
  week = c(18.5, 23, 27.5, 32, 36.5),
  year_f = unique(catch_stock$year_f),
  agg = unique(catch_stock$agg)
) %>% 
  filter(year_f == "2021") 

pred_grid_list <- vector(mode = "list", length = nrow(pred_dat))
for (i in seq_along(pred_grid_list)) {
  pred_grid_list[[i]] <- pred_bathy_grid %>% 
    filter(X < max(catch_stock$xUTM + 1000) & X > min(catch_stock$xUTM - 1000),
           Y < max(catch_stock$yUTM + 1000) & Y > min(catch_stock$yUTM - 1000)) %>% 
    mutate(week = pred_dat$week[i],
           year_f = pred_dat$year_f[i],
           agg = pred_dat$agg[i])
}

pred_grid <- pred_grid_list %>% 
  bind_rows() %>% 
  mutate(
    xUTM_ds = X / 1000,
    yUTM_ds = Y / 1000,
    # slack_z = 0,
    # moon_z = 0,
    month = case_when(
      week < 23 ~ 5,
      week >= 23 & week < 27 ~ 6,
      week >= 27 & week < 32 ~ 7,
      week >= 32 & week < 35 ~ 8,
      week >= 35 ~ 9
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


## FIT PRELIM MODELS -----------------------------------------------------------


# fit MVRW version w/ top four stocks
# fit_full <- sdmTMB(
#   catch ~ 0 + (1 | year_f) + agg + depth_z + slope_z + poly(week_z, 2):agg,
#   offset = "offset",
#   data = catch_stock,
#   mesh = sdm_mesh1,
#   family = sdmTMB::nbinom1(),
#   spatial = "on",
#   # spatial_varying = ~ 0 + size_bin,
#   time = "month",
#   spatiotemporal = "rw",
#   groups = "agg",
#   anisotropy = TRUE,
#   share_range = TRUE, 
#   silent = FALSE
# )
# 
# fit_full2 <- update(fit_full, share_range = FALSE)
# 
# saveRDS(fit_full2, here::here("data", "fit_mvrfrw_stock.rds"))

fit_full <- readRDS(here::here("data", "fit_mvrfrw_stock.rds"))


## MODEL CHECKS ----------------------------------------------------------------

sims_nb1 <- simulate(fit_full, nsim = 100, newdata = catch_stock)
sims_nb1 %>% 
  dharma_residuals(fit_full)


mod_in <- fit_full
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



## FIXED EFFECT PREDS ----------------------------------------------------------

# NOTE month and year values don't matter since random variables and integrated 
# out
pred_foo <- function(x = "week", nd, fit) {
  p <- predict(fit, newdata = nd, se_fit = FALSE, re_form = NA, 
               re_form_iid = NA, nsim = 100)
  nd %>% 
    mutate(
      est = apply(p, 1, mean),
      est_se = apply(p, 1, sd),
      variable = x
    )
}

nd_week <- expand.grid(
  year_f = unique(catch_stock$year_f), 
  week_z = seq(-2, 2, length = 30),
  month = unique(catch_stock$month), 
  slack_z = 0,
  moon_z = 0,
  agg = unique(catch_stock$agg),
  depth_z = 0,
  slope_z = 0
) %>% 
  filter(year_f == "2020", 
         month == "7")

p_week <- pred_foo(x = "week", nd = nd_week, fit = fit_full)


# fixed depth effects
nd_depth <- expand.grid(
  year_f = unique(catch_stock$year_f),
  week_z = 0,
  month = unique(catch_stock$month),
  slack_z = 0,
  moon_z = 0,
  agg = unique(catch_stock$agg),
  depth_z = seq(-2, 2, length = 30),
  slope_z = 0
) %>% 
  filter(agg == "PugetSo", year_f == "2020", month == "7")

p_depth <- pred_foo(x = "depth", nd = nd_depth, fit = fit_full)


# fixed slope effects
nd_slope <- nd_depth %>% 
  mutate(
    slope_z = seq(-2, 2, length = 30),
    depth_z = 0
  )
p_slope <- pred_foo(x = "slope", nd = nd_slope, fit = fit_full)



full_p <- list(p_week, p_depth, p_slope) %>% 
  do.call(rbind, .) %>%
  group_by(agg, variable) %>% 
  mutate(
    week = (week_z * sd(catch_stock$week)) + mean(catch_stock$week),
    depth = (depth_z * sd(catch_stock$mean_depth)) + mean(catch_stock$mean_depth),
    slope = (slope_z * sd(catch_stock$mean_slope)) + mean(catch_stock$mean_slope),
    exp_est = exp(est),
    max_est = max(exp_est),
    scale_est = exp_est / max_est,
    up = (est + (stats::qnorm(0.975) * est_se)),
    lo = (est + (stats::qnorm(0.025) * est_se)),
    exp_up = exp(up),
    exp_lo = exp(lo),
    scale_up = exp(up) / max_est,
    scale_lo = exp(lo) / max_est
  ) 

# plots
plot_foo <- function(var_in = "week", x_lab = "Week") {
  var_in2 <- rlang::enquo(var_in)
  dum <- full_p %>% 
    filter(variable == !!var_in2)
  ggplot(
    dum, 
    aes(x = dum[[var_in]], y = scale_est, ymin = scale_up, 
        ymax = scale_lo)
  ) +
    geom_line() +
    labs(x = x_lab) +
    geom_ribbon(alpha = 0.4) +
    ggsidekick::theme_sleek() +
    coord_cartesian(y = c(0.03, 3)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.title.y = element_blank())
}


week_plot <- plot_foo(var_in = "week", x_lab = "Week") +
  facet_wrap(~agg) +
  theme(axis.title.y = element_text(angle = 90)) +
  ylab("Scaled Abundance")
depth_plot <- plot_foo(var_in = "depth", x_lab = "Bottom Depth (m)")
slope_plot <- plot_foo(var_in = "slope", x_lab = "Bottom Slope (degrees)")

p1 <- cowplot::plot_grid(
  depth_plot, slope_plot,
  ncol = 2
)

png(here::here("figs", "ms_figs_stock", "week_fe.png"), res = 250, units = "in", 
    height = 5.5, width = 5.5)
week_plot
dev.off()

png(here::here("figs", "ms_figs_stock", "other_fes.png"), res = 250, units = "in", 
    height = 5.5, width = 5.5)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    p1, 
    left = grid::textGrob("Scaled Abundance", rot = 90))
)
dev.off()


## SPATIAL PREDICTIONS ---------------------------------------------------------

upsilon_pred <- predict(fit_full,
                        newdata = pred_grid,
                        return_tmb_report = TRUE)

pp <- predict(fit_full, newdata = pred_grid, se_fit = FALSE, re_form = NULL)
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


omega <- pp %>%
  dplyr::select(-c(agg, week, month)) %>%
  distinct() %>% 
  plot_map(., omega_s) +
  scale_fill_gradient2(name = "Spatial\nRF Effect") +
  theme(legend.position = "top")

epsilon <- pp %>%
  plot_map(., upsilon_stc) +
  scale_fill_gradient2(name = "Spatial\nRF Effect") +
  facet_grid(agg~month) +
  theme(legend.position = "top")

full_preds <- pp %>% 
  group_by(agg) %>% 
  mutate(scale_est = exp(est) / max(exp(est))) %>%
  left_join(., week_key, by = "week") %>% 
  plot_map(., scale_est) +
  scale_fill_viridis_c(trans = "sqrt", name = "Scaled\nAbundance") +
  facet_grid(date~agg)


png(here::here("figs", "ms_figs_stock", "omega.png"), res = 250, units = "in", 
    height = 4.5, width = 4.5)
omega
dev.off()

png(here::here("figs", "ms_figs_stock", "epsilon.png"), res = 250, units = "in", 
    height = 7.5, width = 7.5)
epsilon
dev.off()

png(here::here("figs", "ms_figs_stock", "spatial_preds.png"), res = 250, units = "in", 
    height = 7.5, width = 7.5)
full_preds
dev.off()


## subsetted version of above for presentation
full_preds_trim <- pp %>% 
  filter(agg %in% c("Fraser Sub.", "PugetSo")) %>% 
  group_by(agg) %>% 
  mutate(scale_est = exp(est) / max(exp(est))) %>%
  left_join(., week_key, by = "week") %>% 
  plot_map(., scale_est) +
  scale_fill_viridis_c(trans = "sqrt", name = "Scaled\nAbundance") +
  facet_grid(agg~date)  +
  theme(legend.position = "top",
        legend.key.size = unit(0.9, 'cm'))

png(here::here("figs", "spatial_preds_pres.png"), res = 250, 
    units = "in", 
    height = 4.5, width = 7.5)
full_preds_trim
dev.off()