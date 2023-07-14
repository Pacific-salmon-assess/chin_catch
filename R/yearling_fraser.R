## Early Run Fraser Timing
# July 14, 2023
# Use a) model preds from Freshwater et al. 2021 and b) available GSI data
# from tagging program to explore changes in stock composition.

library(tidyverse)


## HISTORICAL DATA -------------------------------------------------------------

# composition model predictions originally published in PeerJ paper
comp_pred_dat <- readRDS(here::here("data",
                                    "peerj_composition_model_predictions.RDS"))

preds <- comp_pred_dat$comp_pred_ci[[3]]

trim_preds <- preds %>% 
  filter(stock %in% c("Fraser_Summer_5.2", "Fraser_Spring_4.2", 
                      "Fraser_Spring_5.2"))

max_month <- trim_preds %>% 
  group_by(stock, region_c) %>% 
  mutate(
    max_month = ifelse(
      pred_prob_est == max(pred_prob_est), TRUE, FALSE
    )
  ) %>% 
  ungroup() %>% 
  filter(max_month == TRUE) %>% 
  select(max_month = month_n, region, stock)


low_prop <- trim_preds %>%
  left_join(., max_month, by = c("region", "stock")) %>% 
  filter(month_n > max_month & pred_prob_est < 0.01,
         !stock == "Fraser_Spring_4.2") %>% 
  group_by(region, stock) %>% 
  summarize(min_month = min(month_n)) 

trim_preds %>%
  filter(month_n == "7.5") %>% 
  select(region, stock, pred_prob_est)


line_comp <- left_join(trim_preds, low_prop, by = c("region", "stock")) %>% 
  filter(month_n > 4 & month_n < 9.3) %>% 
  ggplot(.,
       aes(x = month_n, y = pred_prob_est, colour = region_c)) +
  geom_line() +
  geom_vline(aes(xintercept = min_month, colour = region_c), lty = 2) +
  facet_wrap(~stock) +
  ggsidekick::theme_sleek() + 
  labs(y = "Predicted Proportion of Catch", x = "Month")

png(here::here("figs", "fraser_yearling", "line_preds.png"))
line_comp
dev.off()


## HISTORICAL DATA -------------------------------------------------------------

set_dat <- readRDS(here::here("data", "cleanSetData.RDS")) %>% 
  mutate(nearshore = ifelse(coast_dist_km < 9.26, TRUE, FALSE))

chin_catch <- readRDS(here::here("data", "clean_catch.RDS")) %>% 
  left_join(., set_dat %>% select(event, nearshore), by = "event") %>% 
  mutate(
    mu = case_when(
      agg_name %in% c("WCVI", "ECVI") ~ agg_name,
      cu == "LFR-fall" ~ "Fraser\nFall",
      cu == "LTh" ~ "Fraser\nSpr. 4.2",
      cu %in% c("MFR-spring", "UFR-spring") ~ "Fraser\nSpr. 5.2",
      cu %in% c("MFR-summer", "NTh-sum", "STh-1.3") ~ "Fraser\nSum. 5.2",
      cu %in% c("STh-0.3", "STh-SHUR") ~ "Fraser\nSum. 0.3",
      TRUE ~ cu
    ))

calc_ppn <- function(dat) {
 dat %>%
    group_by(year) %>%
    mutate(year_n = n()) %>% 
    ungroup() %>% 
    group_by(mu2, year, year_n) %>%
    tally() %>%
    mutate(prop = n / year_n) 
}


# total catch
total <- calc_ppn(chin_catch) %>% 
  mutate(data = "all data")
total_nearshore <- chin_catch %>% 
  filter(nearshore == FALSE) %>% 
  calc_ppn() %>% 
  mutate(data = "offshore only")
july <- chin_catch %>% 
  filter(month == "7") %>% 
  calc_ppn() %>% 
  mutate(data = "july only")
july_nearshore <- chin_catch %>% 
  filter(nearshore == FALSE & month == "7") %>% 
  calc_ppn() %>% 
  mutate(data = "july offshore only")
pred_all <- do.call(rbind, list(total, total_nearshore, july, july_nearshore)) %>% 
  ungroup() %>% 
  mutate(
    data = fct_relevel(data, "all data", "offshore only", "july only", 
                       "july offshore only"))

plot_dat <- expand.grid(
  mu2 = unique(chin_catch$mu2),
  year = unique(chin_catch$year),
  data = unique(pred_all$data)
) %>% 
  filter(!mu2 == "other") %>% 
  left_join(
    ., 
    pred_all %>% 
      select(year, year_n, data) %>% 
      distinct(),
    by = c("year", "data")
  ) %>% 
  left_join(., pred_all %>% select(-year_n), by = c("mu2", "year", "data")) %>% 
  mutate(
    n = replace_na(n, 0),
    prop = replace_na(prop, 0)
  ) 

box_ppn <- ggplot(plot_dat, aes(x = data, y = prop)) +
  geom_boxplot() +
  facet_wrap(~mu2, ncol = 1) +
  ggsidekick::theme_sleek() +
  labs(y = "Observed Stock Composition", x = "Dataset")

png(here::here("figs", "fraser_yearling", "obs_comp.png"))
box_ppn
dev.off()

