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
           agg_name %in% c("Fraser", "WCVI", "ECVI") & genetic_source == "GSI" ~ 
             "wild",
           # cu %in% c("LFR-fall", "SWVI", "QP-fall", "CWCH-KOK", "STh-SHUR",
           #   "EVIGStr-sum") & genetic_source == "GSI" ~ "wild",
           agg_name %in% c("PugetSo", "Col") & clip == "N" ~ "wild",
           TRUE ~ "unknown"
         )
         ) %>% 
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
    agg_name = case_when(
      cu == "LCR" ~ "Low Col.",
      cu == "UWR" ~ "Low Col.",
      cu %in% c("LFR-fall", "STh-0.3", "STh-SHUR") ~ "Fraser Sub.",
      cu %in% c("MFR-summer", "NTh-sum", "STh-1.3", "LTh") ~ "Fraser Year.",
      agg == "Fraser" ~ "Fraser Sub.",
      agg == "Col" ~ "Up Col.",
      TRUE ~ agg
    ),
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
    agg_name = fct_relevel(
      agg_name, "WCVI", "ECVI", "Fraser Year.", "Fraser Sub.", "PugetSo",
      "Low Col.", "Up Col.", "WA_OR", "Cali"
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
    fl, size_bin, lipid, origin, genetic_source, stage, stock, cu, agg_name
  ) %>% 
  # remove 2023 data until GSI available
  filter(!year == "2023")


# summary of sample sizes
chin %>% 
  group_by(year) %>% 
  summarize(n_chin = length(unique(fish)),
            n_tag = length(unique(vemco_code)),
            start_date = min(deployment_time, na.rm = T),
            end_date = max(deployment_time, na.rm = T))

saveRDS(chin, here::here("data", "clean_catch.RDS"))


# STOCK COMPOSITION ------------------------------------------------------------

# table of stock breakdown
total_comp_table <- chin %>% 
  group_by(agg_name, cu, stock) %>% 
  tally() %>% 
  print(n = Inf)

# table of clipping rates
clip_table <- chin %>% 
  group_by(agg_name, cu, origin) %>% 
  tally() %>% 
  print(n = Inf)


# overall aggregate comp
full_comp <- chin %>% 
  group_by(size_bin) %>%
  mutate(size_n = n()) %>% 
  ungroup() %>% 
  group_by(agg_name, size_bin, size_n) %>%
  tally() %>%
  mutate(prop = n / size_n) 

full_comp_stacked <- ggplot(full_comp, 
                            aes(fill = agg_name, y = prop, x = size_bin)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d(name = "Stock Aggregate", na.value = "grey60") +
  labs(y = "Proportion Stock Composition", x = "Size Class") +
  ggsidekick::theme_sleek()


# stock composition by month
month_comp <- chin %>% 
  group_by(month) %>%
  mutate(month_n = n()) %>% 
  ungroup() %>% 
  group_by(agg_name, month, month_n) %>%
  tally() %>%
  mutate(prop = n / month_n) 

month_comp_stacked <- ggplot(month_comp, 
                            aes(fill = agg_name, y = prop, x = month)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d(name = "Stock Aggregate", na.value = "grey60" ) +
  labs(y = "Proportion Stock Composition", x = "Month") +
  ggsidekick::theme_sleek() +
  theme(axis.title.y = element_blank())


# composition by adipose clip
origin_comp_stacked <- ggplot(chin, 
                            aes(fill = agg_name, #y = prop, 
                                x = origin)) + 
  geom_bar(position="stack") +#, stat="identity") +
  scale_fill_viridis_d(guide = "none", na.value = "grey60") +
  labs(y = "Proportion Stock Composition", x = "Origin") +
  ggsidekick::theme_sleek() +
  theme(axis.title.y = element_blank())


# composition by maturity stage
stage_comp <- chin %>% 
  group_by(stage) %>%
  mutate(stage_n = n()) %>% 
  ungroup() %>% 
  group_by(agg_name, stage, stage_n) %>%
  tally() %>%
  mutate(prop = n / stage_n) 

stage_comp_stacked <- ggplot(stage_comp, 
                            aes(fill = agg_name, y = prop, 
                                x = stage)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d(guide = "none", na.value = "grey60") +
  labs(x = "Maturity Stage") +
  ggsidekick::theme_sleek() +
  theme(axis.title.y = element_blank())

p1 <- cowplot::plot_grid(
  origin_comp_stacked, stage_comp_stacked,
  ncol = 2
)
y_grob <- grid::textGrob("Proportion Stock Composition", 
                         rot = 90)
p2 <- cowplot::plot_grid(
  month_comp_stacked, p1, nrow = 2
)

png(here::here("figs", "ms_figs", "stock_comp.png"), res = 250, units = "in", 
    height = 5.5, width = 5.5)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(p2, left = y_grob))
dev.off()



# CONDITION --------------------------------------------------------------------

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
lipid_fl <- ggplot(chin %>% filter(!is.na(lipid), !is.na(agg_name))) +
  geom_point(aes(x = fl, y = lipid, fill = stage), shape = 21, alpha = 0.6) +
  scale_fill_viridis_d(option = "B", guide = "none") +
  ggsidekick::theme_sleek() +
  facet_wrap(~agg_name) +
  labs(x = "Fork Length", y = "Lipid Content") +
  theme(legend.position = "top", legend.title = element_blank())

ggplot(chin %>% filter(!is.na(lipid), !is.na(agg_name))) +
  geom_point(aes(x = year_day, y = fl, fill = stage), shape = 21, alpha = 0.6) +
  scale_fill_viridis_d(option = "B", guide = "none") +
  ggsidekick::theme_sleek() +
  facet_wrap(~agg_name) +
  theme(legend.position = "top", legend.title = element_blank())


## size and lipid boxplots 
fl_bp <- chin %>%
  filter(!is.na(agg_name)) %>%
  ggplot(., aes(x = agg_name, y = fl, fill = agg_name)) +
  geom_boxplot() +
  scale_fill_viridis_d(guide = "none") +
  labs(y = "Fork Length") +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

lipid_bp <- chin %>%
  filter(!is.na(agg_name),
         !is.na(lipid)) %>%
  ggplot(., aes(x = agg_name, y = lipid, fill = agg_name)) +
  geom_boxplot() +
  scale_fill_viridis_d(guide = "none") +
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
# lipid_fl
cowplot::plot_grid(
  lipid_fl, p3, ncol = 2
)
dev.off()


## FIT FL/CONDITION MODELs -----------------------------------------------------

library(mgcv)


## fit
chin_fl <- chin %>% 
  filter(!is.na(fl),
         !is.na(agg_name))

fit_fl <- gam(
  fl ~ s(year_day, bs = "tp", m = 2, k = 4) + 
    s(year_day, by = stage, bs = "tp", m = 1, k = 4) + stage + origin + 
    s(agg_name, bs = "re") + s(year, bs = "re"),
  data = chin_fl,
  method = "REML",
  link =
)

chin_lipid <- chin %>% 
  filter(!is.na(lipid),
         !is.na(agg_name))

fit_lipid <- gam(
  lipid ~ s(fl, bs = "tp", m = 2, k = 4) + 
    s(year_day, bs = "tp", m = 2, k = 4) + 
    s(year_day, by = stage, bs = "tp", m = 1, k = 4) + stage + origin + 
    s(agg_name, bs = "re") + s(year, bs = "re"),
  data = chin_fl,
  method = "REML"
  # family = Gamma(link = "log")
)
fit_lipid2 <- gam(
  lipid ~ s(fl, bs = "tp", m = 2, k = 4) + 
    s(year_day, bs = "tp", m = 2, k = 4) + 
    s(year_day, by = stage, bs = "tp", m = 1, k = 4) + stage + origin + 
    s(agg_name, bs = "re") + s(year, bs = "re"),
  data = chin_fl,
  method = "REML"
)

## checks
sim_lipid <- simulate(fit_lipid2, newdata = chin_lipid, nsim = 50) %>% 
  as.matrix()
fix_pred <- predict(fit_lipid2, newdata = chin_lipid, type = "response") %>% 
  as.numeric()
dharma_res <- DHARMa::createDHARMa(
      simulatedResponse = sim_lipid,
      observedResponse = chin_lipid$lipid,
      fittedPredictedResponse = fix_pred
    )
plot(dharma_res)


## predictions

# fl-stage 
fl_stage <- chin %>% 
  group_by(stage) %>% 
  summarize(
    fl = mean(fl, na.rm = T)
  )

# yday preds
new_yday <- expand.grid(
  year_day = seq(120, 260, by = 5),
  stage = unique(chin$stage)
) %>% 
  left_join(., fl_stage, by = "stage") 

fl_preds <- predict(fit_fl, newdata = fit_fl)
