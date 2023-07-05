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


chin <- readRDS(here::here("data", "clean_catch.RDS")) %>% 
  rename(agg = agg_name)


chin %>% 
  group_by(size_bin) %>% 
  tally()
chin %>% 
  group_by(agg) %>% 
  tally()
    

# calculate total catch across size bins
catch_size1 <- chin %>% 
  group_by(event, size_bin) %>% 
  summarize(catch = n(), .groups = "drop") %>% 
  ungroup()
# catch_stock1 <- chin %>%
#   filter(!is.na(agg),
#          fl > 65,
#          # remove rare stocks
#          !agg %in% c("ECVI", "WCVI", "WA_OR", "Cali")) %>% 
#   group_by(event, agg) %>% 
#   summarize(catch = n(), .groups = "drop") %>% 
#   ungroup()


# clean and bind to set data
set_dat <- readRDS(here::here("data", "cleanSetData.RDS")) %>% 
  mutate(
    week = lubridate::week(date_time_local)
  )


catch_size <- expand.grid(event = set_dat$event,
                     size_bin = unique(catch_size1$size_bin)) %>%
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
    month_f = case_when(
      month %in% c(4, 5) ~ 5,
      month %in% c(8, 9) ~ 8,
      TRUE ~ month
    ) %>% 
      as.factor(),
    year_f = as.factor(year),
    yUTM_ds = yUTM / 1000,
    xUTM_ds = xUTM / 1000,
    offset = log(set_dist),
    depth_z = scale(mean_depth)[ , 1],
    slope_z = scale(mean_slope)[ , 1],
    slack_z = scale(hours_from_slack)[ , 1],
    week_z = scale(week)[ , 1],
    moon_z = scale(moon_illuminated)[ , 1]
  ) 

# catch_stock <- expand.grid(
#   event = set_dat$event,
#   agg = unique(catch_stock1$agg)
# ) %>%
#   arrange(event) %>%
#   left_join(., catch_stock1, by = c("event", "agg")) %>%
#   replace_na(., replace = list(catch = 0)) %>%
#   left_join(., set_dat, by = "event") %>%
#   # remove sets not on a troller and tacking
#   filter(!grepl("rec", event),
#          !tack == "yes") %>%
#   mutate(
#     # pool undersampled months
#     month_f = case_when(
#       month %in% c(4, 5) ~ 5,
#       month %in% c(8, 9) ~ 8,
#       TRUE ~ month
#     ) %>% 
#       as.factor(),
#     year_f = as.factor(year),
#     yUTM_ds = yUTM / 1000,
#     xUTM_ds = xUTM / 1000,
#     offset = log(set_dist),
#     depth_z = scale(mean_depth)[ , 1],
#     slope_z = scale(mean_slope)[ , 1],
#     slack_z = scale(hours_from_slack)[ , 1],
#     week_z = scale(week)[ , 1],
#     moon_z = scale(moon_illuminated)[ , 1]
#   ) 
 

# EXP FIGS ---------------------------------------------------------------------

ggplot(catch_size) +
  geom_point(aes(x = week, y = log(catch))) +
  facet_wrap(~size_bin)

ggplot(catch_size) +
  geom_point(aes(x = hours_from_slack, y = log(catch))) +
  facet_wrap(~size_bin)

ggplot(catch_size) +
  geom_point(aes(x = mean_slope, y = log(catch))) +
  facet_wrap(~size_bin)

ggplot(catch_size) +
  geom_point(aes(x = mean_depth, y = log(catch))) +
  facet_wrap(~size_bin)

ggplot(catch_size) +
  geom_point(aes(x = moon_illuminated, y = log(catch))) +
  facet_wrap(~size_bin)



# BUILD MESHES -----------------------------------------------------------------

# construct a few alternative meshes
sdm_mesh1 <- make_mesh(catch_size,
                      c("xUTM_ds", "yUTM_ds"),
                      n_knots = 250)
# sdm_mesh2 <- make_mesh(catch_size,
#                        c("xUTM_ds", "yUTM_ds"),
#                        n_knots = 150)


# SPATIAL PREDICTION GRID ------------------------------------------------------

# import 1.8 X 1.8 km grid generated in prep_bathymetry.R in chinTagging repo
pred_bathy_grid <- readRDS(
  here::here("data", "pred_bathy_grid_utm.RDS"))

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
    hours_from_slack = median(catch_size$hours_from_slack),
    moon_illuminated = 0.5,
    month_f = case_when(
      week <= 22 ~ "5",
      week > 22 & week <= 26 ~ "6",
      week > 26 & week <= 30 ~ "7",
      week > 30 ~ "8"
    )) %>% 
  rename(mean_depth = depth, mean_slope = slope)


# helper function for simple predictive maps
plot_map <- function(dat, column) {
  ggplot() +
    geom_raster(data = dat, aes(X, Y, fill = {{ column }})) +
    coord_fixed() +
    ggsidekick::theme_sleek()
}

ggplot() +
  geom_raster(data = pred_grid, aes(X, Y, fill = mean_depth)) +
  geom_point(data = set_dat %>%
               filter(!grepl("rec", event)), 
             aes(xUTM, yUTM, colour = as.factor(year)),
             alpha = 0.4) +
  coord_fixed() +
  ggsidekick::theme_sleek()


# SPATIAL MODELS ---------------------------------------------------------------

# f6 <- sdmTMB(
#   catch ~ (1 | year_f) + size_bin + poly(hours_from_slack, 2) +
#     mean_depth + poly(moon_illuminated, 2) +
#     mean_slope + week:size_bin,
#   offset = "offset",
#   data = catch_size,
#   mesh = sdm_mesh1,
#   family = sdmTMB::nbinom2(),
#   spatial = "off",
#   spatial_varying = ~ 0 + size_bin + month_f,
#   anisotropy = FALSE,
#   control = sdmTMBcontrol(
#     newton_loops = 1,
#     map = list(
#       ln_tau_Z = factor(
#         c(1, 2, 3, rep(4, times = length(unique(catch_size$month_f)) - 1))
#       )
#     )
#   ),
#   silent = FALSE
# )
# 
# f6_nb1 <- update(f6, family = sdmTMB::nbinom1())

# f6_nb1 <- sdmTMB(
#   catch ~ 0 + (1 | year_f) + size_bin + poly(slack_z, 2) +
#     depth_z + poly(moon_z, 2) +
#     slope_z + week_z:size_bin,
#   offset = "offset",
#   data = catch_size,
#   mesh = sdm_mesh1,
#   family = sdmTMB::nbinom1(),
#   spatial = "off",
#   spatial_varying = ~ 0 + size_bin + month_f,
#   anisotropy = FALSE,
#   control = sdmTMBcontrol(
#     newton_loops = 1,
#     map = list(
#       ln_tau_Z = factor(
#         c(1, 2, 3, rep(4, times = length(unique(catch_size$month_f)) - 1))
#       )
#     )
#   ),
#   silent = FALSE
# )
# saveRDS(f6_nb1, here::here("data", "model_fits", "f6_nb1.rds"))

f6_nb1 <- readRDS(here::here("data", "model_fits", "f6_nb1.rds"))


## SIMULATION CHECKS -----------------------------------------------------------

# quick check
# sims_nb2 <- simulate(f6, nsim = 100)
# sims_nb2 %>% 
#   dharma_residuals(f6)

sims_nb1 <- simulate(f6_nb1, nsim = 100, newdata = catch_size)
sims_nb1 %>% 
  dharma_residuals(f6_nb1)


# sample from posterior using MCMC to avoid Laplace approximation
object <- f6_nb1
samp <- sample_mle_mcmc(object, mcmc_iter = 130, mcmc_warmup = 100)

obj <- object$tmb_obj
random <- unique(names(obj$env$par[obj$env$random]))
pl <- as.list(object$sd_report, "Estimate")
fixed <- !(names(pl) %in% random)
map <- lapply(pl[fixed], function(x) factor(rep(NA, length(x))))
obj <- TMB::MakeADFun(obj$env$data, pl, map = map, DLL = "sdmTMB")
obj_mle <- object
obj_mle$tmb_obj <- obj
obj_mle$tmb_map <- map
ss <- simulate(obj_mle, mcmc_samples = extract_mcmc(samp), nsim = 30)


pred_fixed <- f6_nb1$family$linkinv(
  predict(f6_nb1, newdata = catch_size)$est_non_rf
)
r_nb1_size <- DHARMa::createDHARMa(
  simulatedResponse = ss,
  observedResponse = f6_nb1$data$catch,
  fittedPredictedResponse = pred_fixed
)
DHARMa::testResiduals(r_nb1_size)


## FIXED EFFECT PREDICTIONS ----------------------------------------------------

# quick visualization of effect size estimates 
fes <- tidy(f6_nb1, effects = "fixed", conf.int = T)
fes %>% 
  filter(!term %in% c("size_binlarge", "size_binmedium", "size_binsmall")) %>%
  ggplot(., aes(y = term, x = estimate)) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high))


# quick visualization of effect size estimates 
res <- tidy(f6_nb1, effects = "ran_par", conf.int = T)

res %>% 
  filter(!term %in% c("range")) %>%
  ggplot(., aes(y = term, x = estimate)) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high))



# NOTE month and year values don't matter since random variables and integrated 
# out

# fixed week effects
# nd_week <- expand.grid(
#   year_f = unique(catch_size$year_f), #"2020",
#   week = seq(min(catch_size$week), max(catch_size$week), length = 30),
#   month_f = unique(catch_size$month_f), #"7",
#   hours_from_slack = 0,
#   moon_illuminated = 0.5,
#   size_bin = unique(catch_size$size_bin),
#   mean_depth = median(catch_size$mean_depth),
#   mean_slope = median(catch_size$mean_slope)
# ) %>% 
#   filter(year_f == "2020", 
#          month_f == "7")
nd_week <- expand.grid(
  year_f = unique(catch_size$year_f), #"2020",
  week_z = seq(-2, 2, length = 30),
  month_f = unique(catch_size$month_f), #"7",
  slack_z = 0,
  moon_z = 0,
  size_bin = unique(catch_size$size_bin),
  depth_z = 0,
  slope_z = 0
) %>% 
  filter(year_f == "2020", 
         month_f == "7")

p <- predict(f6_nb1, newdata = nd_week, se_fit = TRUE, 
             re_form = NA, re_form_iid = NA)
p$variable <- "week"


# fixed depth effects
nd_depth <- expand.grid(
  year_f = unique(catch_size$year_f),
  week_z = 0,
  month_f = unique(catch_size$month_f),
  slack_z = 0,
  moon_z = 0,
  size_bin = unique(catch_size$size_bin),
  depth_z = seq(-2, 2, length = 30),
  slope_z = 0
) %>% 
  filter(size_bin == "medium", year_f == "2020", month_f == "7")

p_depth <- predict(f6_nb1, newdata = nd_depth, se_fit = TRUE, 
                   re_form = NA, re_form_iid = NA)
p_depth$variable <- "depth"

# fixed slope effects
nd_slope <- nd_depth %>% 
  mutate(
    slope_z = seq(-2, 2, length = 30),
    depth_z = 0
  )
p_slope <- predict(f6_nb1, newdata = nd_slope, se_fit = TRUE, 
                   re_form = NA, re_form_iid = NA)
p_slope$variable <- "slope"

# fixed slack effects
nd_slack <- nd_depth %>% 
  mutate(
    slack_z = seq(-3, 3, length = 30),
    depth_z = 0
  )
p_slack <- predict(f6_nb1, newdata = nd_slack, se_fit = TRUE, 
                   re_form = NA, re_form_iid = NA)
p_slack$variable <- "slack"

# fixed lunar effects
nd_moon <- nd_depth %>% 
  mutate(
    moon_z = seq(-1.8, 1.22, length = 30),
    depth_z = 0
  )
p_moon <- predict(f6_nb1, newdata = nd_moon, se_fit = TRUE, 
                  re_form = NA, re_form_iid = NA)
p_moon$variable <- "moon"


full_p <- list(p, p_depth, p_slope, p_slack, p_moon) %>% 
  do.call(rbind, .) %>%
  # p %>% 
  group_by(size_bin, variable) %>% 
  mutate(
    week = (week_z * sd(catch_size$week)) + mean(catch_size$week),
    depth = (depth_z * sd(catch_size$mean_depth)) + mean(catch_size$mean_depth),
    slope = (slope_z * sd(catch_size$mean_slope)) + mean(catch_size$mean_slope),
    slack = (slack_z * sd(catch_size$hours_from_slack)) +
      mean(catch_size$hours_from_slack),
    moon = (moon_z * sd(catch_size$moon_illuminated)) +
      mean(catch_size$moon_illuminated),
    exp_est = exp(est),
    max_est = max(exp_est),
    scale_est = exp_est / max_est,
    up = (est + 1.96 * est_se),
    lo = (est - 1.96 * est_se),
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
    coord_cartesian(y = c(0.03, 2)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.title.y = element_blank())
}

week_plot <- plot_foo(var_in = "week", x_lab = "Week") +
  facet_wrap(~size_bin) +
  theme(axis.title.y = element_text(angle = 90)) +
  ylab("Scaled Abundance")
depth_plot <- plot_foo(var_in = "depth", x_lab = "Bottom Depth (m)")
slope_plot <- plot_foo(var_in = "slope", x_lab = "Bottom Slope (degrees)")
slack_plot <- plot_foo(var_in = "slack", x_lab = "Hours From Slack")
moon_plot <- plot_foo(var_in = "moon", 
                      x_lab = "Proportion Moon Illuminated")

p1 <- cowplot::plot_grid(
  depth_plot, slope_plot, slack_plot, moon_plot,
  ncol = 2
)

png(here::here("figs", "ms_figs", "week_fe.png"), res = 250, units = "in", 
    height = 2.75, width = 5.5)
week_plot
dev.off()

png(here::here("figs", "ms_figs", "other_fes.png"), res = 250, units = "in", 
    height = 5.5, width = 5.5)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(p1, 
                         left = grid::textGrob("Scaled Abundance", rot = 90)))
dev.off()



## SPATIAL PREDICTIONS ---------------------------------------------------------

pp <- predict(f6_nb1, newdata = pred_grid, se_fit = FALSE, re_form = NULL)

plot_map(pp, zeta_s_size_binlarge) +
  scale_fill_gradient2()
plot_map(pp, zeta_s_size_binmedium) +
  scale_fill_gradient2()
plot_map(pp, zeta_s_size_binsmall) +
  scale_fill_gradient2()

plot_map(pp, zeta_s_month_f6) +
  scale_fill_gradient2()
plot_map(pp, zeta_s_month_f7) +
  scale_fill_gradient2()
plot_map(pp, zeta_s_month_f8) +
  scale_fill_gradient2()

plot_map(pp, exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)") +
  facet_grid(size_bin~week)



# PREDICTION GRID --------------------------------------------------------------

# import 1.8 X 1.8 km grid generated in prep_bathymetry.R in chinTagging repo
pred_bathy_grid <- readRDS(
  here::here("data", "pred_bathy_grid_utm.RDS"))

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
    yUTM_ds = Y / 1000#,
    # coast_dist_km = shore_dist / 1000,
    # hours_from_slack = median(catch$hours_from_slack),
    # moon_illuminated = 0.5,
    # year_f = as.factor(year)
  ) %>% 
  rename(mean_depth = depth, mean_slope = slope)#%>% 
  # dplyr::select(
  #   X, Y, xUTM_ds, yUTM_ds, month, year, year_f, year_day, hours_from_slack, 
  #   moon_illuminated, mean_depth = depth, mean_slope = slope
  # )


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



#Predictions incorporating all fixed and random effects
# Generally the impacts of different meshes considered here were negligible
month_preds <- plot_map(pp, exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)") +
  facet_grid(size_bin~week) 
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