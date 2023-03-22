### Non-tag analyses
# Initial data exploration of non tag-related data (i.e. gear, injuries, stock
# composition)
# Nov. 1, 2021

library(tidyverse)

# Import all fish tagdata
chin <- readRDS(here::here("data", "tagging_data", "cleanTagData_GSI.RDS")) %>% 
  mutate(
    mature = ifelse(fl > 65, "mature", "immature"),
    year = as.factor(year),
    agg_name = case_when(
      cu == "LCR" ~ "Low Col.",
      agg_name == "Col" ~ "Up Col.",
      TRUE ~ agg_name
    ),
    mu = case_when(
      agg_name %in% c("WCVI", "ECVI") ~ agg_name,
      cu == "LFR-fall" ~ "Fraser\nFall",
      cu == "LTh" ~ "Fraser\nSpr. 1.2",
      cu %in% c("MFR-summer", "NTh-sum", "STh-1.3") ~ "Fraser\nSum. 1.3",
      cu %in% c("STh-0.3", "STh-SHUR") ~ "Fraser\nSum. 0.3",
      TRUE ~ cu
    )
  )

# summary of sample sizes
chin %>% 
  group_by(year) %>% 
  summarize(n_chin = length(unique(fish)),
            n_tag = length(unique(acoustic_year)),
            start_date = min(deployment_time, na.rm = T),
            end_date = max(deployment_time, na.rm = T))


# STOCK COMPOSITION ------------------------------------------------------------

# table of stock breakdown
total_comp_table <- chin %>% 
  group_by(agg_name, cu, stock) %>% 
  tally() %>% 
  print(n = Inf)

# table of tagged stock breakdown
tag_comp_table <- chin %>% 
  filter(!is.na(acoustic_year)) %>% 
  group_by(agg_name, cu, stock) %>% 
  tally() %>% 
  print(n = Inf)

# table of clipping rates
clip_table <- chin %>% 
  group_by(agg_name, cu, clip) %>% 
  tally() %>% 
  print(n = Inf)

# tagged fish aggregate comp
tag_comp <- chin %>% 
  filter(!is.na(acoustic),
         !is.na(agg_name),
         mature == "mature") %>%
  group_by(year) %>%
  mutate(year_n = n()) %>% 
  ungroup() %>% 
  group_by(agg_name, year, year_n) %>%
  tally() %>%
  ungroup() %>%
  mutate(prop = n / year_n,
         group = "tagged mature")

# overall aggregate comp
full_comp <- chin %>% 
  filter(!is.na(agg_name),
         mature == "mature") %>% 
  group_by(year) %>%
  mutate(year_n = n()) %>% 
  ungroup() %>% 
  group_by(agg_name, year, year_n) %>%
  tally() %>%
  mutate(prop = n / year_n,
         group = "overall mature") %>% 
  rbind(., tag_comp)

full_comp_stacked <- ggplot(full_comp, aes(fill=agg_name, y=prop, x=group)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d() +
  facet_wrap(~year) +
  ggsidekick::theme_sleek()

## composition in early part of season
early_comp_stacked <- chin %>% 
  filter(year_day < 167,
         !is.na(agg_name)) %>% 
  group_by(mature) %>%
  mutate(mature_n = n()) %>%
  ungroup() %>%
  group_by(agg_name, year, mature, mature_n) %>%
  tally() %>%
  mutate(prop = n / mature_n) %>% 
  ggplot(., aes(fill=agg_name, y=prop, x=mature)) + 
  geom_bar(position="stack", stat = "identity") +
  scale_fill_viridis_d() 


full_tally <- chin %>%
  filter(!is.na(agg_name)) %>% 
  mutate(
    country = case_when(
      agg_name %in% c("Fraser", "WCVI", "ECVI") ~ "Canada",
      TRUE ~ "USA"
    ),
    agg_name = recode(agg_name, PugetSo = "PS", Col = "Col.", 
                      Cali = "CA", WA_OR = "WA/OR", Fraser = "Fra.")
  ) %>% 
  group_by(country, agg_name, year, mature) %>% 
  tally()

comp_bar <- ggplot(full_tally) +
  geom_bar(aes(x = reorder(agg_name, -n), y = n, fill = country), 
           stat = "identity") +
  scale_fill_manual(values = c("#e41a1c", "#377eb8"),
                    name = "Country\nof Origin") +
  labs(x = "Stock Aggregate", y = "Count") +
  ggsidekick::theme_sleek() +
  facet_grid(mature~year) +
  theme(
    legend.position = "top"
  )

png(here::here("figs", "exp", "stock_comp_bar.png"), height = 5.5, 
    width = 11, units = "in", res = 200)
comp_bar
dev.off()


# BC composition
bc_tally <- chin %>%
  filter(!is.na(agg_name),
         agg_name %in% c("Fraser", "WCVI", "ECVI")) %>% 
  group_by(agg_name, mature, mu, year) %>% 
  tally()

bc_comp_bar <- ggplot(bc_tally) +
  geom_bar(aes(x = reorder(mu, -n), y = n, fill = agg_name), 
           stat = "identity") +
  scale_fill_manual(values = c("#de2d26", "#fc9272", "#fee0d2"), 
                    name = "Stock\nAggregate") +
  labs(x = "Management Unit", y = "Count") +
  ggsidekick::theme_sleek() +
  facet_grid(mature~year)

pdf(here::here("figs", "exp", "stock_comp.pdf"))
full_comp_stacked
bc_comp_bar
dev.off()


# maturity composition
png(here::here("figs", "exp", "comp_maturity_bar.png"), height = 5.5, 
    width = 10, units = "in", res = 300)
chin %>% 
  filter(!is.na(agg_name)) %>% 
  group_by(year, mature) %>%
  mutate(year_n = n()) %>% 
  ungroup() %>% 
  group_by(agg_name, mature, year, year_n) %>%
  tally() %>%
  mutate(prop = n / year_n) %>% 
  ggplot(., aes(fill=agg_name, y=prop, x=mature)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d() +
  facet_wrap(~year) +
  ggsidekick::theme_sleek()
dev.off()


# stock composition of immature
agg_list <- chin %>%
  filter(!is.na(cu),
         mature == "immature") %>% 
  group_by(cu, agg_name, clip) %>% 
  tally() %>% 
  group_by(agg_name) %>% 
  group_nest()

imm_list <- purrr::map2(
  agg_list$agg_name, agg_list$data,
  function (title, x) {
    ggplot(x %>% droplevels()) +
      geom_bar(aes(x = reorder(cu, -n), y = n, fill = cu), 
               stat = "identity") +
      scale_fill_brewer(palette = which(agg_list$agg_name == title)) +
      labs(x = "Stock", y = "Count", main = title) +
      ggsidekick::theme_sleek() +
      facet_wrap(~clip)
  }
)


png(here::here("figs", "exp", "stock_comp_immature_clip.png"), height = 7.5, 
    width = 13, units = "in", res = 300)
cowplot::plot_grid(imm_list[[1]], imm_list[[2]], imm_list[[3]], imm_list[[4]],
                   imm_list[[5]])
dev.off()
 

# CONDITION --------------------------------------------------------------------

## Body size by year
ggplot(chin %>% filter(!is.na(agg_name))) +
  geom_boxplot(aes(x = year, y = fl, fill = mature)) +
  ggsidekick::theme_sleek() +
  facet_grid(mature~agg_name)


## Clip and body size
ggplot(chin) +
  geom_boxplot(aes(x = clip, y = fl)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~agg_name)


## Energy density
# Fork length
ggplot(chin %>% filter(!is.na(mean_log_e))) +
  geom_point(aes(x = fl, y = mean_log_e, fill = year), shape = 21) +
  scale_fill_viridis_d() +
  ggsidekick::theme_sleek() +
  facet_wrap(~agg_name)

# Stocks
ggplot(chin %>% filter(!is.na(mean_log_e))) +
  geom_boxplot(aes(x = year, y = mean_log_e, fill = mature)) +
  scale_fill_viridis_d(option = "A") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none") +
  facet_grid(mature~agg_name)

chin %>% 
  filter(agg_name %in% c("ECVI", "WCVI", "Fraser"),
         !is.na(mean_log_e),
         mature == "mature") %>% 
  ggplot(.) +
  geom_boxplot(aes(x = mu, y = mean_log_e, fill = mu)) +
  scale_fill_viridis_d() +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none") 


# Capture date
ggplot(chin) +
  geom_boxplot(aes(x = agg_name, y = year_day, fill = agg_name)) +
  scale_fill_viridis_d(option = "A") +
  ggsidekick::theme_sleek() +
  facet_wrap(mature~year)


# Injuries 
inj <- chin %>% 
  filter(!is.na(injury)) %>% 
  select(fish, injury:fin_dam, year, gear_type) %>% 
  pivot_longer(., cols = injury:fin_dam, names_to = "condition", 
               values_to = "score") 

ggplot(inj) +
  geom_bar(aes(x = condition, y = as.factor(score), fill = as.factor(score)), 
           stat = "identity") +
  facet_wrap(~year)

inj %>% 
  group_by(year, gear_type, condition) %>% 
  summarize(
    n = length(fish),
    mean_score = mean(score),
    se_score = sd(score) / sqrt(n),
    up = mean_score + 1.96 * se_score,
    lo = mean_score - 1.96 * se_score
  ) %>% 
  ggplot(.) +
  geom_pointrange(aes(x = condition, y = mean_score, ymin = lo, ymax = up,
                      fill = gear_type), 
                  shape = 21,
                  position = position_dodge(width=0.75)) +
  facet_wrap(~year)


no_tag <- chin %>% 
  filter(
    year == "2022",
    # fl > 70,
    is.na(acoustic)#,
    # injury == "3"
  )

no_tag %>% 
  group_by(hook_loc) %>% 
  tally()


# SIZE BY STOCK ----------------------------------------------------------------

# colour by tagging year
ggplot(chin) +
  geom_histogram(aes(x = fl, fill = year)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~mu)

# colour by whether tagged or not
chin %>% 
  mutate(surv_tag = case_when(
    acoustic_type == "V13P" ~ "yes",
    TRUE ~ "no"
  )) %>% 
  ggplot() +
  geom_histogram(aes(x = fl, fill = surv_tag)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~agg_name)
