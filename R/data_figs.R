## Non-model figures
# 1) Fork length histograms
# 2) Box plots of size and energy density by stock
# 3) Stock composition by size bin and month


library(tidyverse)


## CLEAN -----------------------------------------------------------------------

# Import stage data and fitted model generated in gen_detection_histories.R
stage_dat <- readRDS(here::here("data", "agg_lifestage_df.RDS")) %>% 
  dplyr::select(vemco_code, stage)


chin_raw <- readRDS(here::here("data", "cleanTagData_GSI.RDS")) %>%
  mutate(year_day = lubridate::yday(date),
         # define hatchery
         origin = case_when(
           genetic_source == "PBT" ~ "hatchery",
           clip == "Y" ~ "hatchery",
           grepl("Fraser", agg_name) & genetic_source == "GSI" ~ 
             "wild",
           agg_name %in% c("Fraser", "WCVI", "ECVI") & genetic_source == "GSI" ~ 
             "wild",
           # cu %in% c("LFR-fall", "SWVI", "QP-fall", "CWCH-KOK", "STh-SHUR",
           #   "EVIGStr-sum") & genetic_source == "GSI" ~ "wild",
           agg_name %in% c("Puget Sound", "Up Col.","Low Col.") & clip == "N" ~ 
             "wild",
           # stocks in WA_OR coastal that are in WA and therefore def mass marked
           stock %in% c("FORKS_CREEK_HATCHERY", "SOL_DUC_RIVER", 
                        "TRASK_HATCHERY", "ELWHA_FALL") & clip == "N" ~ "wild",
           TRUE ~ "unknown"
         )) %>% 
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
    stage = ifelse(is.na(stage) | stage == "unknown", 
                   stage_predicted,
                   # paste("predicted", stage_predicted, sep = "_"),
                   stage),
    month = lubridate::month(date),
    year = as.factor(year),
    # mu = case_when(
    #   agg_name %in% c("WCVI", "ECVI") ~ agg_name,
    #   cu == "LFR-fall" ~ "Fraser\nFall",
    #   cu == "LTh" ~ "Fraser\nSpr. 1.2",
    #   cu %in% c("MFR-summer", "NTh-sum", "STh-1.3") ~ "Fraser\nSum. 1.3",
    #   cu %in% c("STh-0.3", "STh-SHUR") ~ "Fraser\nSum. 0.3",
    #   TRUE ~ cu
    # ),
    size_bin = case_when(
      fl < 65 ~ "small",
      fl >= 65 & fl < 75 ~ "medium",
      fl >= 75 ~ "large"
    ),
    agg_name = fct_recode(agg, NULL = "Fraser Sub."),
    agg = factor(
      agg, levels = c(
        "WCVI", "ECVI", "Fraser Year.", "Fraser 4.1", "Fraser Fall", 
        "Puget Sound", "Low Col.", "Up Col.", "WA_OR", "Cali"
      )
    ),
    origin = fct_relevel(
      origin, "hatchery", "wild", "unknown"
    ),
    stage = fct_relevel(
      stage, "mature", "immature"
    )
    ) %>% 
  dplyr::select(
    fish, vemco_code, event, date, year, month, year_day, deployment_time,
    fl, size_bin, lipid, origin, genetic_source, stage, stock, cu, agg
  ) 

# summary of sample sizes
chin %>% 
  group_by(year) %>% 
  summarize(n_chin = length(unique(fish)),
            n_tag = length(unique(vemco_code)),
            start_date = min(deployment_time, na.rm = T),
            end_date = max(deployment_time, na.rm = T))

saveRDS(chin, here::here("data", "clean_catch.RDS"))


# PREP INDIVIDUAL DATASETS -----------------------------------------------------

# clean and bind to set data
set_dat <- readRDS(here::here("data", "cleanSetData.RDS")) %>% 
  mutate(
    week = lubridate::week(date_time_local),
    # pool undersampled months
    month = case_when(
      month %in% c(4, 5) ~ 5,
      TRUE ~ month
    ),
    month_f = month %>% 
      as.factor()
  )

## size data (combine sublegal and adult)
# import bycatch data representing sublegal catch
chin_juv <- read.csv(
  here::here("data", "bycatchData.csv"),
  stringsAsFactors = F) %>% 
  filter(species == "chinook") %>% 
  mutate(size_bin = "sublegal") %>% 
  dplyr::select(
    event, size_bin, catch = count
  )

# calculate total catch across size bins
catch_size1 <- chin %>% 
  group_by(event, size_bin) %>% 
  summarize(catch = n(), .groups = "drop") %>% 
  ungroup() %>% 
  rbind(., chin_juv)

catch_size <- expand.grid(
  event = set_dat$event,
  size_bin = unique(catch_size1$size_bin)
) %>%
  arrange(event) %>%
  left_join(., catch_size1, by = c("event", "size_bin")) %>%
  replace_na(., replace = list(catch = 0)) %>%
  left_join(., set_dat, by = "event")

saveRDS(catch_size, here::here("data", "catch_size_pre.rds"))


## stock comp data
catch_stock1 <- chin  %>% 
  filter(
    !is.na(agg),
    !agg == "Fraser Sub.",
    size_bin %in% c("large", "medium")
  ) %>%
  group_by(event, agg) %>%
  summarize(catch = n(), .groups = "drop") %>%
  ungroup() %>% 
  droplevels()

catch_stock <- expand.grid(
  event = set_dat$event,
  agg = unique(catch_stock1$agg)
) %>%
  arrange(event) %>%
  left_join(., catch_stock1, by = c("event", "agg")) %>%
  replace_na(., replace = list(catch = 0)) %>%
  left_join(., set_dat, by = "event")

saveRDS(catch_stock, here::here("data", "catch_stock_pre.rds"))


## origin data
catch_origin1 <- chin  %>%
  filter(size_bin %in% c("large", "medium")) %>% 
  group_by(event, origin) %>%
  summarize(catch = n(), .groups = "drop") %>%
  ungroup()

catch_origin <- expand.grid(
  event = set_dat$event,
  origin = unique(catch_origin1$origin)
) %>%
  arrange(event) %>%
  left_join(., catch_origin1, by = c("event", "origin")) %>%
  replace_na(., replace = list(catch = 0)) %>%
  left_join(., set_dat, by = "event")

saveRDS(catch_origin, here::here("data", "catch_origin_pre.rds"))


# SUPP TABLE OF STOCK BREAKDOWN ------------------------------------------------

stock_table <- chin %>% 
  filter(!is.na(stock),
         !stock == "Fraser Sub.") %>% 
  arrange(agg, cu, stock) %>% 
  dplyr::select(stock_aggregate = agg, conservation_unit = cu, stock) %>% 
  mutate(
    stock_aggregate = stock_aggregate %>% 
      str_replace_all(., "\n", " "),
    stock = stock %>% 
      str_replace_all(., "_", " ") %>% 
      str_to_title()
  ) %>% 
  distinct()

write.csv(
  stock_table,
  here::here("data", "supp_stock_table.csv"),
  row.names = FALSE
)


# STOCK COMPOSITION ------------------------------------------------------------

# table of stock breakdown
total_comp_table <- chin %>% 
  group_by(agg, cu, stock) %>% 
  tally() %>% 
  print(n = Inf)

# table of origin
clip_table <- chin %>% 
  group_by(agg, cu, origin) %>% 
  tally() %>% 
  print(n = Inf)


# overall aggregate comp
full_comp <- chin %>% 
  group_by(size_bin) %>%
  mutate(size_n = n()) %>% 
  ungroup() %>% 
  group_by(agg, size_bin, size_n) %>%
  tally() %>%
  mutate(prop = n / size_n) 

full_comp_stacked <- ggplot(full_comp, 
                            aes(fill = agg, y = prop, x = size_bin)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d(name = "Stock Aggregate", na.value = "grey60") +
  labs(y = "Proportion Stock Composition", x = "Size Class") +
  ggsidekick::theme_sleek()


# stock composition by month
month_comp <- chin %>%
  group_by(month) %>%
  mutate(month_n = n()) %>%
  ungroup() %>%
  group_by(agg, month, month_n) %>%
  tally() %>%
  mutate(prop = n / month_n)

month_comp_stacked <- ggplot(month_comp,
                            aes(fill = agg, y = prop, x = month)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d(name = "Stock Aggregate", na.value = "grey60" ) +
  labs(y = "Proportion Stock Composition", x = "Month") +
  ggsidekick::theme_sleek()

# month_count_stacked <- ggplot(chin, 
#        aes(fill = agg, #y = prop, 
#            x = month)) + 
#   geom_bar(position="stack") +#, stat="identity") +
#   scale_fill_viridis_d(name = "Stock Aggregate", na.value = "grey60" ) +
#   labs(y = "Number of Individuals", x = "Origin") +
#   ggsidekick::theme_sleek()  +
#   theme(axis.title.y = element_blank())

# composition by origin (supplmentary)
origin_count_stacked <- ggplot(
  chin, 
  aes(fill = agg, #y = prop, 
      x = origin)
) + 
  geom_bar(position="stack") +#, stat="identity") +
  scale_fill_viridis_d(guide = "none", na.value = "grey60") +
  labs(x = "Origin", y = "Number of Individuals") +
  ggsidekick::theme_sleek() 

# origin by size class
# origin_size_comp <- chin %>% 
#   group_by(size_bin) %>%
#   mutate(size_n = n()) %>% 
#   ungroup() %>% 
#   group_by(origin, size_bin, size_n) %>%
#   tally() %>%
#   mutate(prop = n / size_n) 
# 
# ggplot(origin_size_comp, 
#        aes(fill = origin, y = prop, x = size_bin)) + 
#   geom_bar(position="stack", stat="identity") +
#   scale_fill_viridis_d(name = "Origin", na.value = "grey60" ) +
#   labs(y = "Proportion Origin Composition", x = "Size Bin") +
#   ggsidekick::theme_sleek() +
#   theme(axis.title.y = element_blank())


# composition by maturity stage
# stage_comp <- chin %>% 
#   group_by(stage) %>%
#   mutate(stage_n = n()) %>% 
#   ungroup() %>% 
#   group_by(agg_name, stage, stage_n) %>%
#   tally() %>%
#   mutate(prop = n / stage_n) 
# 
# stage_comp_stacked <- ggplot(stage_comp, 
#                             aes(fill = agg_name, y = prop, 
#                                 x = stage)) + 
#   geom_bar(position="stack", stat="identity") +
#   scale_fill_viridis_d(guide = "none", na.value = "grey60") +
#   labs(x = "Maturity Stage") +
#   ggsidekick::theme_sleek() +
#   theme(axis.title.y = element_blank())

stage_count_stacked <- ggplot(chin, 
                               aes(fill = agg, #y = prop, 
                                   x = stage)) + 
  geom_bar(position="stack") +#, stat="identity") +
  scale_fill_viridis_d(guide = "none", na.value = "grey60") +
  labs(x = "Maturity Stage", y = "Number of Individuals") +
  ggsidekick::theme_sleek() +
    theme(axis.title.y = element_blank())


p1 <- cowplot::plot_grid(
  origin_count_stacked, stage_count_stacked,
  ncol = 2
)
# y_grob <- grid::textGrob("Number of Individuals", 
#                          rot = 90)
# gridExtra::grid.arrange(
#   gridExtra::arrangeGrob(p1, left = y_grob))

# p2 <- cowplot::plot_grid(
#   month_comp_stacked, p1, nrow = 2
# )

png(here::here("figs", "ms_figs", "stock_comp.png"), res = 250, units = "in", 
    height = 5.5, width = 5.5)
# gridExtra::grid.arrange(
#   gridExtra::arrangeGrob(p2, left = y_grob))
cowplot::plot_grid(
  month_comp_stacked, p1, nrow = 2
)
dev.off()


# CONDITION --------------------------------------------------------------------

chin <- readRDS(here::here("data", "clean_catch.RDS"))
chin$agg2 <- str_replace(chin$agg, " ", "\n")

length_hist <- ggplot(chin, aes(x = fl, fill = stage)) +
  geom_histogram(col=I("black")) +
  geom_vline(aes(xintercept = 65), linetype = 1) +
  geom_vline(aes(xintercept = 75), linetype = 2) +
  scale_fill_viridis_d(option = "B") +
  labs(x = "Fork Length") +
  ggsidekick::theme_sleek() +
  theme(axis.title.y = element_blank())
 
png(here::here("figs", "ms_figs", "fl_histogram.png"), res = 250, units = "in", 
    height = 4, width = 5.5)
length_hist
dev.off()


## lipid content and fork length
lipid_fl <- ggplot(
  chin %>% 
    filter(!is.na(lipid), !is.na(agg), !agg == "Fraser Sub.")
) +
  geom_point(aes(x = fl, y = lipid, fill = stage), shape = 21, alpha = 0.6) +
  scale_fill_brewer(type = "div", palette = 5) +
  ggsidekick::theme_sleek() +
  facet_wrap(~agg) +
  labs(x = "Fork Length", y = "Lipid Content") +
  theme(legend.position = "top", legend.title = element_blank())


## size and lipid boxplots 
fl_bp <- chin %>%
  filter(!is.na(agg), !agg == "Fraser Sub.") %>%
  ggplot(., aes(x = fct_reorder(agg2, as.numeric(agg)),
                y = fl, fill = stage)) +
  geom_boxplot() +
  scale_fill_brewer(type = "div", palette = 5, guide = "none") +
  labs(y = "Fork Length") +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

lipid_bp <- chin %>%
  filter(!is.na(agg),
         !is.na(lipid),
         !agg == "Fraser Sub.") %>%
  ggplot(., aes(x =  fct_reorder(agg2, as.numeric(agg)), y = lipid, 
                fill = stage)) +
  geom_boxplot() +
  scale_fill_brewer(type = "div", palette = 5, guide = "none") +
  labs(y = "Lipid Content") +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank()
  )


p3 <- cowplot::plot_grid(
  fl_bp, lipid_bp, ncol = 1
)

png(here::here("figs", "ms_figs", "fl_lipid.png"), res = 250, units = "in", 
    height = 5, width = 8.5)
cowplot::plot_grid(
  lipid_fl, p3, ncol = 2
)
dev.off()


## FIT FL/CONDITION MODELs -----------------------------------------------------

library(mgcv)

## fit
chin_fl <- chin %>% 
  filter(!is.na(fl),
         !is.na(agg))

fit_fl <- gam(
  fl ~ s(year_day, bs = "tp", m = 2, k = 4) + 
    s(year_day, by = stage, bs = "tp", m = 1, k = 4) + stage + origin + 
    s(year_day, by = agg, bs = "tp", m = 1, k = 4) +  
    s(agg, bs = "re") + s(year, bs = "re"),
  data = chin_fl,
  method = "REML"
)

chin_lipid <- chin %>% 
  filter(!is.na(lipid),
         !is.na(agg))

fit_lipid <- gam(
  lipid ~ s(fl, bs = "tp", m = 2, k = 4) + 
    s(fl, by = agg, bs = "tp", m = 1, k = 4) +  
    s(year_day, bs = "tp", m = 2, k = 4) + 
    s(year_day, by = stage, bs = "tp", m = 1, k = 4) +
    stage + origin + 
    s(year_day, by = agg, bs = "tp", m = 1, k = 4) +  
    s(agg, bs = "re") + s(year, bs = "re"),
  data = chin_lipid,
  method = "REML",
  family = Gamma(link = "log")
)

## checks
sim_lipid <- simulate(fit_lipid, newdata = chin_lipid, nsim = 50) %>% 
  as.matrix()
fix_pred <- predict(fit_lipid, newdata = chin_lipid, type = "response") %>% 
  as.numeric() 
dharma_res <- DHARMa::createDHARMa(
      simulatedResponse = sim_lipid,
      observedResponse = chin_lipid$lipid,
      fittedPredictedResponse = fix_pred
    )
plot(dharma_res)

sim_fl <- simulate(fit_fl, newdata = chin_fl, nsim = 50) %>% 
  as.matrix()
fix_pred <- predict(fit_fl, newdata = chin_fl, type = "response") %>% 
  as.numeric() 
dharma_res <- DHARMa::createDHARMa(
  simulatedResponse = sim_fl,
  observedResponse = chin_fl$fl,
  fittedPredictedResponse = fix_pred
)
plot(dharma_res)


## posterior predictions
# lipid_post <- fl_post <- vector(length = 8, mode = "list")
# for (i in 1:8) {
#   lipid_post[[i]] <- chin_lipid %>%
#     mutate(
#       pred_lipid = sim_lipid[ , i]
#     ) %>%
#     pivot_longer(
#       cols = c(lipid, pred_lipid)
#     ) %>%
#     ggplot(.) +
#     geom_boxplot(aes(x = agg, y = value, fill = name)) +
#     facet_wrap(~stage, scales = "free_x")
#   fl_post[[i]] <- chin_fl %>%
#     mutate(
#       pred_fl = sim_fl[ , i]
#     ) %>%
#     pivot_longer(
#       cols = c(fl, pred_fl)
#     ) %>%
#     ggplot(.) +
#     geom_boxplot(aes(x = agg, y = value, fill = name)) +
#     facet_wrap(~stage, scales = "free_x")
# }
# 
# cowplot::plot_grid(lipid_post[[1]], lipid_post[[2]], lipid_post[[3]],
#                    lipid_post[[4]], lipid_post[[5]], lipid_post[[6]],
#                    ncol = 1)
# cowplot::plot_grid(fl_post[[1]], fl_post[[2]], fl_post[[3]],
#                    fl_post[[4]], fl_post[[5]], fl_post[[6]],
#                    ncol = 1)

## fixed effect sizes
fl_fe <- broom::tidy(fit_fl, conf.int = TRUE, parametric = TRUE) %>% 
  mutate(
    response = "fork length"
  )
lipid_fe <- broom::tidy(fit_lipid, conf.int = TRUE, parametric = TRUE) %>% 
  mutate(
    response = "lipid"
  )

fes <- rbind(fl_fe, lipid_fe) %>% 
  mutate(
    covariate = ifelse(grepl("origin", term), "origin", "maturity\nstage"),
    term = fct_recode(
      as.factor(term), "immature" = "stageimmature", "wild" = "originwild",
      "unknown" = "originunknown"
    )
  ) %>% 
  filter(
    !term == "(Intercept)"
  )

shape_pal <- c(21, 22)
names(shape_pal) <- unique(fes$covariate)

fe_plot <- ggplot(fes) +
  geom_pointrange(
    aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high, 
        shape = covariate)
  ) +
  facet_wrap(~response, scales = "free_y") +
  ggsidekick::theme_sleek() +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red")

png(here::here("figs", "ms_figs", "fes.png"), res = 250, units = "in", 
    height = 3.25, width = 6.5)
fe_plot
dev.off()



## predictions

# fl-stage 
fl_stage <- chin %>% 
  filter(stage == "mature") %>%
  # group_by(stage) %>% 
  summarize(
    fl = mean(fl, na.rm = T)
  ) 

# yday preds
new_yday <- expand.grid(
  year_day = seq(120, 260, by = 5),
  stage = unique(chin$stage),
  origin = unique(chin$origin),
  # since predicting only on fixed effects factor levels don't matter
  agg = "ECVI",
  year = "2022",
  fl = fl_stage$fl
) 

fl_preds <- predict(
  fit_fl, newdata = new_yday, se.fit = TRUE,
  # include only FE terms
  terms = c("(Intercept)", "stageimmature", "originwild", "originunknown",
            "s(year_day)", "s(year_day):stagemature", 
            "s(year_day):stageimmature")
)
lipid_preds_yday <- predict(
  fit_lipid, newdata = new_yday, se.fit = TRUE, type = "response",
  # include only FE terms
  terms = c("(Intercept)", "stageimmature", "originwild", "originunknown",
            "s(year_day)", "s(year_day):stagemature", 
            "s(year_day):stageimmature", "s(fl)")
)

new_dat <- new_yday %>% 
  mutate(
    pred_fl = as.numeric(fl_preds$fit),
    lo_fl = pred_fl + (stats::qnorm(0.025) * as.numeric(fl_preds$se.fit)), 
    up_fl = pred_fl + (stats::qnorm(0.975) * as.numeric(fl_preds$se.fit)),
    pred_lipid = as.numeric(lipid_preds_yday$fit),
    lo_lipid = pred_lipid + (stats::qnorm(0.025) * 
                               as.numeric(lipid_preds_yday$se.fit)), 
    up_lipid = pred_lipid + (stats::qnorm(0.975) * 
                               as.numeric(lipid_preds_yday$se.fit))
  )


# main
stage_pal <- c("#1f78b4", "#a6cee3")
names(stage_pal) <- unique(new_dat$stage)
  
yday_fl <- ggplot(new_dat %>% filter(origin == "hatchery")) +
  geom_line(aes(x = year_day, y = pred_fl, colour = stage)) +
  geom_ribbon(aes(x = year_day, ymin = lo_fl, ymax = up_fl, fill = stage), 
              alpha = 0.3) +
  scale_fill_manual(values = stage_pal, guide = "none") +
  scale_colour_manual(values = stage_pal, guide = "none") +
  labs(y = "Predicted Fork Length") +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(121, 182, 244),
                     labels = c("May 1", "July 1", "Sep 1")) +
  scale_y_continuous(expand = c(0, 0)) 

yday_lipid <- ggplot(new_dat %>% filter(origin == "hatchery")) +
  geom_line(aes(x = year_day, y = pred_lipid, colour = stage)) +
  geom_ribbon(aes(x = year_day, ymin = lo_lipid, ymax = up_lipid, fill = stage), 
              alpha = 0.3) +
  scale_fill_manual(values = stage_pal) +
  scale_colour_manual(values = stage_pal) +
  labs(y = "Predicted Lipid Content") +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(121, 182, 244),
                     labels = c("May 1", "July 1", "Sep 1")) +
  scale_y_continuous(expand = c(0, 0)) 



yday_legend <- cowplot::get_legend(yday_lipid)
yday_plots1 <- cowplot::plot_grid(
  yday_fl, 
  yday_lipid +
    theme(legend.position = "none"),
  ncol = 2
)
# x_grob_yday <- grid::textGrob("Date")
# yday_plots2 <- gridExtra::grid.arrange(
#   gridExtra::arrangeGrob(yday_plots1, bottom = x_grob_yday))

png(here::here("figs", "ms_figs", "yday_smooths.png"), res = 250, units = "in",
    height = 3.5, width = 6.5)
cowplot::plot_grid(
  yday_legend, 
  yday_plots1,
  ncol = 1,
  rel_heights = c(0.075, 1)
)
dev.off()


# fl preds
new_fl <- expand.grid(
  fl = seq(65, 95, by = 1),
  # one stage since maturity confounded with size
  stage = unique(chin$stage),
  origin = unique(chin$origin),#unique(chin$origin),
  # since predicting only on fixed effects factor levels don't matter
  agg = "ECVI",
  year = "2022",
  year_day = mean(chin$year_day)
) 

lipid_preds_fl <- predict(
  fit_lipid, newdata = new_fl, se.fit = TRUE, type = "response",
  # include only FE terms
  terms = c("(Intercept)", "stageimmature", "originwild", "originunknown",
            "s(year_day)", "s(year_day):stagemature", 
            "s(fl)")
)

new_dat2 <- new_fl %>% 
  mutate(
    pred_lipid = as.numeric(lipid_preds_fl$fit),
    lo_lipid = pred_lipid + (stats::qnorm(0.025) * 
                               as.numeric(lipid_preds_fl$se.fit)), 
    up_lipid = pred_lipid + (stats::qnorm(0.975) * 
                               as.numeric(lipid_preds_fl$se.fit))
  )

fl_lipid <- ggplot(new_dat2 %>% filter(stage == "mature", origin == "hatchery")) +
  geom_line(aes(x = fl, y = pred_lipid, colour = stage)) +
  geom_ribbon(aes(x = fl, ymin = lo_lipid, ymax = up_lipid, fill = stage), 
              alpha = 0.3) +
  scale_fill_manual(values = stage_pal, guide = "none") +
  scale_colour_manual(values = stage_pal, guide = "none") +
  ggsidekick::theme_sleek() +
  labs(x = "Fork Length", y = "Lipid Content") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) 


png(here::here("figs", "ms_figs", "fl_smooths.png"), res = 250, units = "in",
    height = 4.25, width = 4.25)
fl_lipid
dev.off()


# stock preds
mean_fl <- chin %>% 
  filter(stage == "mature") %>% 
  summarize(
    fl = mean(fl, na.rm = T)
  ) %>% 
  pull(fl)
new_stock <- expand.grid(
  year_day = mean(chin$year_day),
  stage = unique(chin$stage),
  origin = unique(chin$origin),
  # since predicting only on fixed effects factor levels don't matter
  agg = unique(chin$agg),
  year = "2022",
  fl = mean_fl
) %>% 
  filter(stage == "mature", origin == "hatchery", !is.na(agg)) 
new_stock$agg2 <- str_replace(new_stock$agg, " ", "\n")

fl_preds_agg <- predict(
  fit_fl, newdata = new_stock, se.fit = TRUE,
  exclude = "s(year)"
)
lipid_preds_agg <- predict(
  fit_lipid, newdata = new_stock, se.fit = TRUE, type = "response",
  exclude = "s(year)"
  )


new_dat3 <- new_stock %>% 
  mutate(
    pred_fl = as.numeric(fl_preds_agg$fit),
    lo_fl = pred_fl + (stats::qnorm(0.025) * 
                               as.numeric(fl_preds_agg$se.fit)), 
    up_fl = pred_fl + (stats::qnorm(0.975) * 
                               as.numeric(fl_preds_agg$se.fit)),
    pred_lipid = as.numeric(lipid_preds_agg$fit),
    lo_lipid = pred_lipid + (stats::qnorm(0.025) * 
                               as.numeric(lipid_preds_agg$se.fit)), 
    up_lipid = pred_lipid + (stats::qnorm(0.975) * 
                               as.numeric(lipid_preds_agg$se.fit))
  )

fl_agg <- ggplot(new_dat3) +
  geom_pointrange(
    aes(x = agg2, y = pred_fl, ymin = lo_fl, ymax = up_fl, fill = agg2),
    shape = 21
  ) +
  scale_fill_viridis_d(guide = "none") +
  labs(y = "Predicted Fork Length") +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank()
  )

lipid_agg <- ggplot(new_dat3) +
  geom_pointrange(
    aes(x = agg2, y = pred_lipid, ymin = lo_lipid, ymax = up_lipid, fill = agg2),
    shape = 21
  ) +
  scale_fill_viridis_d(guide = "none") +
  labs(y = "Predicted Lipid Content") +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank()
  )

png(here::here("figs", "ms_figs", "agg_cond_ests.png"), res = 250, units = "in",
    height = 4.25, width = 4.25)
cowplot::plot_grid(
  fl_agg, lipid_agg, ncol = 1
)
dev.off()
