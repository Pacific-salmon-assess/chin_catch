## Fit size-based spatial distribution models
# August 7, 2020
# updated Feb 22, 2023
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


## prep size data
catch_size1 <- readRDS(here::here("data", "catch_size_pre.rds"))

# filter and adjust data made in data_figs.R
catch_size <- catch_size1 %>%
  # remove sets not on a troller and tacking
  filter(!grepl("rec", event),
         !tack == "yes") %>% 
  mutate(bin = as.factor(size_bin))


## prep stock data (small individuals already removed)
catch_stock1 <- readRDS(here::here("data", "catch_stock_pre.rds"))

catch_stock <- catch_stock1 %>%
  # remove sets not on a troller and tacking
  filter(!agg %in% c("ECVI", "Fraser Year."),
         !grepl("rec", event),
         !tack == "yes") %>% 
  mutate(
    bin = as.factor(agg)
  ) %>% 
  droplevels()


## prep origin data (small individuals already removed)
catch_origin1 <- readRDS(here::here("data", "catch_origin_pre.rds"))

catch_origin <- catch_origin1 %>%
  # remove sets not on a troller and tacking
  filter(!grepl("rec", event),
         !tack == "yes") %>% 
  mutate(bin = as.factor(origin))


# join and modify 
dat_tbl <- tibble(
  dataset = c("size", "stock", "origin")
) 
dat_tbl$data <- purrr::map(
  list(catch_size, catch_stock, catch_origin),
  ~ .x %>% 
    mutate(
      year_f = as.factor(year),
      yUTM_ds = yUTM / 1000,
      xUTM_ds = xUTM / 1000,
      effort = (set_dist / 1000) * lines,
      offset = log((set_dist / 1000) * lines),
      # effort of one corresponds to ~2 lines towed for 500 m
      depth_z = scale(mean_depth)[ , 1],
      slope_z = scale(mean_slope)[ , 1],
      slack_z = scale(hours_from_slack)[ , 1],
      week_z = scale(week)[ , 1],
      moon_z = scale(moon_illuminated)[ , 1],
      dist_z = scale(coast_dist_km)[ , 1]
    )
)


# figure colour palette
size_main <- "#377eb8"
origin_main <- "#4daf4a"
stock_main <- "#984ea3"


# BUILD MESHES -----------------------------------------------------------------

dat_tbl$mesh <- purrr::map(
  dat_tbl$data,
  ~ make_mesh(.x,
              c("xUTM_ds", "yUTM_ds"),
              n_knots = 125)
)


# SPATIAL PREDICTION GRID ------------------------------------------------------

# import 1 x 1 km grid based on shrunk CH shape files and generated in 
# prep_prediction_grid.R 
pred_bathy_grid <- readRDS(
  here::here("data", "pred_bathy_grid_1000m.RDS"))
trim_grid <- pred_bathy_grid %>% 
  filter(X < max(catch_size$xUTM + 1000) & X > min(catch_size$xUTM - 1000),
         Y < max(catch_size$yUTM + 1000) & Y > min(catch_size$yUTM - 1000))

pred_dat_size <- expand.grid(
  week = c(18.5, 23, 27.5, 32, 36.5),
  year_f = dat_tbl$data[[1]]$year_f[[1]], # year shouldn't matter since RI
  size_bin = unique(dat_tbl$data[[1]]$size_bin)
) 
pred_list_size <- vector(mode = "list", length = nrow(pred_dat_size))
for (i in seq_along(pred_list_size)) {
  pred_list_size[[i]] <- trim_grid %>% 
    mutate(week = pred_dat_size$week[i],
           year_f = pred_dat_size$year_f[i],
           size_bin = pred_dat_size$size_bin[i])
}

pred_dat_stock <- expand.grid(
  week = c(18.5, 23, 27.5, 32, 36.5),
  year_f = dat_tbl$data[[1]]$year_f[[1]], # year shouldn't matter since RI
  agg = unique(dat_tbl$data[[2]]$agg)
) 
pred_list_stock <- vector(mode = "list", length = nrow(pred_dat_stock))
for (i in seq_along(pred_list_stock)) {
  pred_list_stock[[i]] <- trim_grid %>% 
    mutate(week = pred_dat_stock$week[i],
           year_f = pred_dat_stock$year_f[i],
           agg = pred_dat_stock$agg[i])
}

pred_dat_origin <- expand.grid(
  week = c(18.5, 23, 27.5, 32, 36.5),
  year_f = dat_tbl$data[[1]]$year_f[[1]], # year shouldn't matter since RI
  origin = unique(dat_tbl$data[[3]]$origin)
) 
pred_list_origin <- vector(mode = "list", length = nrow(pred_dat_origin))
for (i in seq_along(pred_list_origin)) {
  pred_list_origin[[i]] <- trim_grid %>% 
    mutate(week = pred_dat_origin$week[i],
           year_f = pred_dat_origin$year_f[i],
           origin <- pred_dat_origin$origin[i])
}

# clean lists above
dat_tbl$pred_data <- purrr::map2(
  list(pred_list_size, pred_list_stock, pred_list_origin),
  dat_tbl$data,
  ~ .x %>% 
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
      week_z = (week - mean(.y$week)) / sd(.y$week),
      depth_z = (depth - mean(.y$mean_depth)) / sd(.y$mean_depth),
      slope_z = (slope - mean(.y$mean_slope)) / sd(.y$mean_slope)
    ) 
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



# FIT MODELS -------------------------------------------------------------------

fit_size <- sdmTMB(
  catch ~ 0 + (1 | year_f) + bin + poly(slack_z, 2) +
    depth_z + poly(moon_z, 2) +
    slope_z + poly(week_z, 2):bin,
  offset = "offset",
  data = dat_tbl$data[[1]],
  mesh = dat_tbl$mesh[[1]],
  family = sdmTMB::nbinom1(),
  spatial = "on",
  # spatial_varying = ~ 0 + size_bin,
  time = "month",
  spatiotemporal = "rw",
  groups = "bin",
  anisotropy = TRUE,
  share_range = FALSE,
  silent = FALSE
)

# CONSIDER MORE COMPLEX MODEL (1 share range off, 2 slack/moon included) w/ more
# stocks
fit_stock <- sdmTMB(
  catch ~ 0 + (1 | year_f) + bin +# poly(slack_z, 2) +
    depth_z +# poly(moon_z, 2) +
    slope_z + poly(week_z, 2):bin,
  offset = "offset",
  data = dat_tbl$data[[2]],
  mesh = dat_tbl$mesh[[2]],
  family = sdmTMB::nbinom1(),
  spatial = "on",
  time = "month",
  spatiotemporal = "rw",
  groups = "bin",
  anisotropy = TRUE,
  share_range = TRUE,
  silent = FALSE
)

fit_origin <- sdmTMB(
  catch ~ 0 + (1 | year_f) + bin +# poly(slack_z, 2) +
    depth_z + #poly(moon_z, 2) +
    slope_z + poly(week_z, 2):bin,
  offset = "offset",
  data = dat_tbl$data[[3]],
  mesh = dat_tbl$mesh[[3]],
  family = sdmTMB::nbinom1(),
  spatial = "on",
  time = "month",
  spatiotemporal = "rw",
  groups = "bin",
  anisotropy = TRUE,
  share_range = FALSE,
  silent = FALSE
)

saveRDS(fit_size, here::here("data", "model_fits", "fit_mvrfrw_size.rds"))
saveRDS(fit_stock, here::here("data", "model_fits", "fit_mvrfrw_stock.rds"))
saveRDS(fit_origin, here::here("data", "model_fits", "fit_mvrfrw_origin.rds"))


fit_size <- readRDS(here::here("data", "model_fits", "fit_mvrfrw_size.rds"))
fit_stock <- readRDS(here::here("data", "model_fits", "fit_mvrfrw_stock.rds"))
fit_origin <- readRDS(here::here("data", "model_fits", "fit_mvrfrw_origin.rds"))

fit_list <- list(fit_size, fit_stock, fit_origin)


## SIMULATION CHECKS -----------------------------------------------------------

## TODO: check all groups

sims_nb1 <- simulate(fit, nsim = 100, newdata = catch_size)
sims_nb1 %>% 
  dharma_residuals(fit)


# sample from posterior using MCMC to avoid Laplace approximation
object <- fit
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


pred_fixed <- fit$family$linkinv(
  predict(fit, newdata = catch_size)$est_non_rf
)
r_nb1_size <- DHARMa::createDHARMa(
  simulatedResponse = ss,
  observedResponse = fit$data$catch,
  fittedPredictedResponse = pred_fixed
)
DHARMa::testResiduals(r_nb1_size)
# all looks good


## FIXED EFFECTS ---------------------------------------------------------------

# pull parameter estimates
dat_tbl$fes <- purrr::map2(
  fit_list,
  dat_tbl$dataset,
  ~ tidy(.x, effects = "fixed", conf.int = T) %>% 
    filter(!grepl("poly", term),
           !(grepl("_z", term))) %>% 
    mutate(
      term = str_replace(term, "bin", "") %>% 
        as.factor(),
      exp_est = exp(estimate),
      exp_lo = exp(conf.low),
      exp_hi = exp(conf.high)
    ) 
)
dat_tbl$fes[[3]] <- dat_tbl$fes[[3]] %>% 
  mutate(term = fct_relevel(term, "unknown", after = Inf))
dat_tbl$res <- purrr::map2(
  fit_list,
  dat_tbl$dataset,
  ~ tidy(.x, effects = "ran_par", conf.int = T) 
)


# intercept plots
int_plots <- purrr::map2(
  dat_tbl$fes,
  list(size_main, stock_main, origin_main),
  ~ ggplot(.x) + 
    geom_pointrange(aes(x = term, y = exp_est, ymin = exp_lo, ymax = exp_hi),
                    shape = 21, fill = .y) +
    ggsidekick::theme_sleek() +
    labs(y = "Mean Catch Rate") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle=45, vjust=1, hjust=1)
    )
)

# slope estimates (supp table)
slope_dat <- dat_tbl %>% 
  dplyr::select(dataset, fes) %>% 
  unnest(cols = fes) %>% 
  filter(term %in% c("depth_z", "slope_z"))
  


# helper funcs
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

plot_foo <- function(dat_in, var_in = "week", x_lab = "Week", 
                     col_in = "darkgrey") {
  var_in2 <- rlang::enquo(var_in)
  dum <- dat_in %>% 
    filter(variable == !!var_in2)
  ggplot(
    dum, 
    aes(x = dum[[var_in]], y = scale_est, ymin = scale_up, ymax = scale_lo)
  ) +
    geom_line(colour = col_in) +
    labs(x = x_lab, y = "Scaled Abundance") +
    geom_ribbon(alpha = 0.3, fill = col_in) +
    ggsidekick::theme_sleek() +
    coord_cartesian(y = c(0.03, 3)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) 
}


dat_tbl$week_preds <- purrr::map2(
  dat_tbl$data,
  fit_list,
  function (x, y) {
    dum <- expand.grid(
      year_f = "2020",
      week_z = seq(-2, 2, length = 30),
      month = 7,
      slack_z = 0,
      moon_z = 0,
      bin = unique(x$bin),
      depth_z = 0,
      slope_z = 0
    ) 
    
    pred_foo(x = "week", nd = dum, fit = y) %>% 
      group_by(bin) %>% 
      mutate(
        week = (week_z * sd(x$week)) + mean(x$week),
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
  }
)


week_plots <- purrr::map2(
  dat_tbl$week_preds,
  list(size_main, stock_main, origin_main),
  ~ plot_foo(dat_in = .x, var_in = "week", x_lab = "Week", col_in = .y) +
    facet_wrap(~bin) +
    theme(
      legend.position = "none"
    )
) 


png(here::here("figs", "ms_figs", "size_week_fe.png"), res = 250, units = "in", 
    height = 4.5, width = 7.5)
cowplot::plot_grid(
  int_plots[[1]], week_plots[[1]], ncol = 2, rel_widths = c(0.55, 1)
)
dev.off()

png(here::here("figs", "ms_figs", "stock_week_fe.png"), res = 250, units = "in", 
    height = 4.5, width = 7.5)
cowplot::plot_grid(
  int_plots[[2]], week_plots[[2]], ncol = 2, rel_widths = c(0.55, 1)
)
dev.off()

png(here::here("figs", "ms_figs", "origin_week_fe.png"), res = 250, units = "in", 
    height = 4.5, width = 7.5)
cowplot::plot_grid(
  int_plots[[3]], week_plots[[3]], ncol = 2, rel_widths = c(0.55, 1)
)
dev.off()



# fixed depth effects
nd_depth <- expand.grid(
  year_f = unique(catch_size$year_f),
  week_z = 0,
  month = unique(catch_size$month),
  slack_z = 0,
  moon_z = 0,
  size_bin = unique(catch_size$size_bin),
  depth_z = seq(-2, 2, length = 30),
  slope_z = 0
) %>% 
  filter(size_bin == "medium", year_f == "2020", month == "7")

p_depth <- pred_foo(x = "depth", nd = nd_depth, fit = fit)


# fixed slope effects
nd_slope <- nd_depth %>% 
  mutate(
    slope_z = seq(-2, 2, length = 30),
    depth_z = 0
  )
p_slope <- pred_foo(x = "slope", nd = nd_slope, fit = fit)


# fixed slack effects
nd_slack <- nd_depth %>% 
  mutate(
    slack_z = seq(-3, 3, length = 30),
    depth_z = 0
  )
p_slack <- pred_foo(x = "slack", nd = nd_slack, fit = fit)


# fixed lunar effects
nd_moon <- nd_depth %>% 
  mutate(
    moon_z = seq(-1.8, 1.22, length = 30),
    depth_z = 0
  )
p_moon <- pred_foo(x = "moon", nd = nd_moon, fit = fit)


full_p <- list(p_week, p_depth, p_slope, p_slack, p_moon) %>% 
  do.call(rbind, .) %>%
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
    up = (est + (stats::qnorm(0.975) * est_se)),
    lo = (est + (stats::qnorm(0.025) * est_se)),
    exp_up = exp(up),
    exp_lo = exp(lo),
    scale_up = exp(up) / max_est,
    scale_lo = exp(lo) / max_est
  ) 



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
    height = 5.5, width = 5.5)
week_plot
dev.off()

png(here::here("figs", "ms_figs", "other_fes.png"), res = 250, units = "in", 
    height = 5.5, width = 5.5)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    p1, 
    left = grid::textGrob("Scaled Abundance", rot = 90))
  )
dev.off()



## SPATIAL PREDICTIONS ---------------------------------------------------------

upsilon_pred <- predict(fit,
                        newdata = pred_grid,
                        return_tmb_report = TRUE)

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

