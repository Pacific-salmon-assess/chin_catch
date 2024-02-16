## Fit size-based spatial distribution models
# August 7, 2020
# updated Feb 22, 2023
# Given cross-validation (catch_sdm_cv.R) and goals of analysis, utilize
# month focused spatio temporal model
# Also include static spatial variables (e.g. depth/slope) and dynamic 
# fine-scale (e.g. time of day) 


# ensure mvrfrw branch installed
devtools::install_github("https://github.com/pbs-assess/sdmTMB",
                         ref = "mvrfrw")

library(dplyr)
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
    bin = as.factor(agg),
    mig = ifelse(
      agg %in% c("Fraser Fall", "Puget Sound", "Low Col."), "res", "mig")
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
      dist_z = scale(coast_dist_km)[ , 1],
      sunrise_z = scale(time_since_sunrise)[ , 1]
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
              n_knots = 175)
)


# SPATIAL PREDICTION GRID ------------------------------------------------------

# import 1 x 1 km grid based on shrunk CH shape files and generated in 
# prep_prediction_grid.R 
pred_bathy_grid <- readRDS(
  here::here("data", "pred_bathy_grid_1000m.RDS"))
trim_grid <- pred_bathy_grid %>% 
  filter(X < max(catch_size$xUTM + 1500) & X > min(catch_size$xUTM - 1500),
         Y < max(catch_size$yUTM + 1500) & Y > min(catch_size$yUTM - 1500))


# clean lists above
dat_tbl$pred_data <- purrr::map(
  dat_tbl$data,
  function (x) {
    pred_dat <- expand.grid(
      week = c(18.5, 23, 27.5, 32, 36.5),
      year_f = x$year_f[[1]], # year shouldn't matter since RI
      bin = unique(x$bin)
    ) 
    pred_list <- vector(mode = "list", length = nrow(pred_dat))
    for (i in seq_along(pred_list)) {
      pred_list[[i]] <- trim_grid %>% 
        mutate(week = pred_dat$week[i],
               year_f = pred_dat$year_f[i],
               bin = pred_dat$bin[i])
    }
    
    pred_list %>% 
      bind_rows() %>% 
      mutate(
        xUTM_ds = X / 1000,
        yUTM_ds = Y / 1000,
        slack_z = 0,
        sunrise_z = 0,
        moon_z = 0,
        month = case_when(
          week < 23 ~ 5,
          week >= 23 & week < 27 ~ 6,
          week >= 27 & week < 32 ~ 7,
          week >= 32 & week < 35 ~ 8,
          week >= 35 ~ 9
        ),
        month_f = as.factor(month),
        week_z = (week - mean(x$week)) / sd(x$week),
        depth_z = (depth - mean(x$mean_depth)) / sd(x$mean_depth),
        slope_z = (slope - mean(x$mean_slope)) / sd(x$mean_slope)
        ) 
  }
)

# key to assign weeks dates
week_key <- data.frame(
  week = unique(dat_tbl$pred_data[[1]]$week)
) %>% 
  mutate(
    date = c("May 1", "June 5", "July 10", "Aug 15", "Sep 10") %>% 
      as.factor() %>% 
      fct_reorder(., week)
  )


# HELPER FUNCTIONS -------------------------------------------------------------

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


# FIT MODELS -------------------------------------------------------------------

## remove slack tide effects until 2023 data available

fit_size <- sdmTMB(
  catch ~ 0 + (1 | year_f) + bin + #poly(slack_z, 2) +
    depth_z + poly(moon_z, 2) + sunrise_z +
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


# CONSIDER MORE COMPLEX MODEL (2 slack/moon included)
# interaction with mig not supported by AIC (delta < 2)
fit_stock <- sdmTMB(
  catch ~ 0 + (1 | year_f) + bin +# poly(slack_z, 2) +
    depth_z + sunrise_z +# poly(moon_z, 2) +
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
  share_range = FALSE,
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
names(fit_list) <- c("size", "stock", "origin")


## SIMULATION CHECKS -----------------------------------------------------------

qq_list <- purrr::map2(
  fit_list,
  dat_tbl$data,
  function(x, y) {
    simulate(x, nsim = 100, newdata = y) %>% 
      dharma_residuals(x)
  }
)

# spatiotemporal resolution of fits
resid_plots <- purrr::map2(
  fit_list, names(fit_list), function (x, y) {
    dum <- x$data %>% 
      mutate(resid = resid(x))
    ggplot(dum) +
      geom_point(aes(x = xUTM_ds, y = yUTM_ds, colour = resid)) +
      facet_grid(month_f ~ year_f) +
      scale_colour_gradient2() +
      labs(title = y) +
      ggsidekick::theme_sleek()
  }
)

pdf(here::here("figs", "diagnostics", "resid_spatial.pdf"))
resid_plots
dev.off()


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
purrr::map2(
  fit_list,
  dat_tbl$dataset,
  ~ tidy(.x, effects = "fixed", conf.int = T) %>% 
    filter(term %in% c("depth_z", "slope_z")) %>% 
    mutate(dataset = .y)
) %>% 
  bind_rows()
  


dat_tbl$week_preds <- purrr::map2(
  dat_tbl$data,
  fit_list,
  function (x, y) {
    dum <- expand.grid(
      year_f = "2020",
      week_z = seq(-2, 2, length = 30),
      month = 7,
      sunrise_z = 0,
      # slack_z = 0,
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



# fixed depth effects (only for size model since similar with others)
nd_depth <- expand.grid(
  year_f = unique(dat_tbl$data[[1]]$year_f)[[1]],
  week_z = 0,
  month = 7, 
  slack_z = 0,
  moon_z = 0,
  bin = unique(dat_tbl$data[[1]]$bin),
  depth_z = seq(-2, 2, length = 30),
  slope_z = 0,
  sunrise_z = 0
) %>% 
  filter(bin == "medium")
p_depth <- pred_foo(x = "depth", nd = nd_depth, fit = fit_size)

# fixed slope effects
nd_slope <- nd_depth %>% 
  mutate(
    slope_z = seq(-2, 2, length = 30),
    depth_z = 0
  )
p_slope <- pred_foo(x = "slope", nd = nd_slope, fit = fit_size)

# fixed slack effects
# nd_slack <- nd_depth %>% 
#   mutate(
#     slack_z = seq(-3, 3, length = 30),
#     depth_z = 0
#   )
# p_slack <- pred_foo(x = "slack", nd = nd_slack, fit = fit_size)

# fixed lunar effects
nd_moon <- nd_depth %>% 
  mutate(
    moon_z = seq(-1.8, 1.22, length = 30),
    depth_z = 0
  )
p_moon <- pred_foo(x = "moon", nd = nd_moon, fit = fit_size)

# fixed sunrise effects
nd_sunrise <- nd_depth %>% 
  mutate(
    sunrise_z = seq(-2.1, 3.2, length = 30),
    moon_z = 0
  )
p_sunrise <- pred_foo(x = "sunrise", nd = nd_sunrise, fit = fit_size)


full_p <- list(p_depth, p_slope, #p_slack, 
               p_moon, p_sunrise) %>% 
  do.call(rbind, .) %>%
  group_by(variable) %>% 
  mutate(
    depth = (depth_z * sd(catch_size$mean_depth)) + mean(catch_size$mean_depth),
    slope = (slope_z * sd(catch_size$mean_slope)) + mean(catch_size$mean_slope),
    # slack = (slack_z * sd(catch_size$hours_from_slack)) +
    #   mean(catch_size$hours_from_slack),
    moon = (moon_z * sd(catch_size$moon_illuminated)) +
      mean(catch_size$moon_illuminated),
    sunrise = (sunrise_z * sd(catch_size$time_since_sunrise)) +
      mean(catch_size$time_since_sunrise),
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


depth_plot <- plot_foo(dat_in  = full_p, var_in = "depth", 
                       x_lab = "Bottom Depth (m)", col_in = size_main) +
  theme(axis.title.y = element_blank())
slope_plot <- plot_foo(dat_in  = full_p, var_in = "slope", 
                       x_lab = "Bottom Slope (degrees)", col_in = size_main) +
  theme(axis.title.y = element_blank())
# slack_plot <- plot_foo(dat_in  = full_p, var_in = "slack", 
#                        x_lab = "Hours From Slack", col_in = size_main) +
#   theme(axis.title.y = element_blank())
moon_plot <- plot_foo(dat_in  = full_p, var_in = "moon", 
                      x_lab = "Proportion Moon Illuminated",
                      col_in = size_main) +
  theme(axis.title.y = element_blank())
sunrise_plot <- plot_foo(dat_in  = full_p, var_in = "sunrise", 
                      x_lab = "Hours From Sunrise",
                      col_in = size_main) +
  theme(axis.title.y = element_blank())

p1 <- cowplot::plot_grid(
  depth_plot, slope_plot, #slack_plot, 
  moon_plot, sunrise_plot,
  ncol = 2
)


png(here::here("figs", "ms_figs", "size_other_fes.png"), res = 250, units = "in", 
    height = 5.5, width = 5.5)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    p1, 
    left = grid::textGrob("Scaled Abundance", rot = 90))
  )
dev.off()



## SPATIAL PREDICTIONS ---------------------------------------------------------


spatial_preds <- purrr::map2(
  fit_list,
  dat_tbl$pred_data,
  function (x, y) {
    upsilon_pred <- predict(x,
                            newdata = y,
                            return_tmb_report = TRUE)
    pp <- predict(x, newdata = y, se_fit = FALSE, re_form = NULL, 
                  re_form_iid = NA)
    pp$upsilon_stc <- as.numeric(upsilon_pred$proj_upsilon_st_A_vec)
    return(pp)
  }
)


omegas <- purrr::map(
  spatial_preds,
  ~ .x %>%
    dplyr::select(-c(bin, week, month)) %>%
    distinct() %>% 
    plot_map(., omega_s) +
    scale_fill_gradient2(name = "Spatial\nRF Effect") +
    # facet_grid(month~bin) +
    theme(legend.position = "top") 
)

epsilons <- purrr::map(
  spatial_preds,
  ~ .x %>% 
    plot_map(., upsilon_stc) +
    scale_fill_gradient2(name = "Spatial\nRF Effect") +
    facet_grid(bin ~ month) +
    theme(legend.position = "top") 
  )

full_preds <- purrr::map(
  spatial_preds,
  ~ .x %>%  
    group_by(bin) %>% 
    mutate(scale_est = exp(est) / max(exp(est))) %>%
    left_join(., week_key, by = "week") %>% 
    plot_map(., scale_est) +
    scale_fill_viridis_c(trans = "sqrt", name = "Scaled\nAbundance") +
    facet_grid(bin ~ date)  +
    theme(legend.position = "top",
          legend.key.size = unit(.85, 'cm')) 
)


# size figs
png(here::here("figs", "ms_figs", "size_omega.png"), res = 250, units = "in", 
    height = 4, width = 4)
omegas[[1]]
dev.off()

png(here::here("figs", "ms_figs", "size_epsilon.png"), res = 250, units = "in", 
    height = 7.5, width = 7.5)
epsilons[[1]]
dev.off()

png(here::here("figs", "ms_figs", "size_spatial_preds.png"), res = 250,
    units = "in", height = 7.5, width = 7.5)
full_preds[[1]]
dev.off()


# stock figs
png(here::here("figs", "ms_figs", "stock_omega.png"), res = 250, units = "in", 
    height = 4, width = 4)
omegas[[2]]
dev.off()

png(here::here("figs", "ms_figs", "stock_epsilon.png"), res = 250, units = "in", 
    height = 7.5, width = 7.5)
epsilons[[2]]
dev.off()

png(here::here("figs", "ms_figs", "stock_spatial_preds.png"), res = 250,
    units = "in", height = 7.5, width = 7.5)
full_preds[[2]]
dev.off()


# origin figs
png(here::here("figs", "ms_figs", "origin_omega.png"), res = 250, units = "in", 
    height = 4, width = 4)
omegas[[3]]
dev.off()

png(here::here("figs", "ms_figs", "origin_epsilon.png"), res = 250, units = "in", 
    height = 7.5, width = 7.5)
epsilons[[3]]
dev.off()

png(here::here("figs", "ms_figs", "origin_spatial_preds.png"), res = 250,
    units = "in", height = 7.5, width = 7.5)
full_preds[[3]]
dev.off()

