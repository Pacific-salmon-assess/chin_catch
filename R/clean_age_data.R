### Clean aging data
## Oct 26, 2023
## 2019-21 scales aged by David Huff contractor
## Clean and compare to stock ID data, then export for catch analysis

library(tidyverse)


chin <- readRDS(here::here("data", "cleanTagData_GSI.RDS"))

age_in <- read.csv(here::here("data", "aging_data", "age_estimates.csv"),
                   stringsAsFactors = FALSE, na.strings=c("", "NA")) %>% 
  janitor::clean_names() %>% 
  select(fish = fish_i_d, european_age, age_confidence = age_confidence_1_4) %>% 
  # remove blank rows, resorbed scale samples and low confidence ages
  filter(!is.na(fish),
         !european_age == "R",
         !age_confidence == "1") 

age_in %>% 
  group_by(european_age) %>% 
  tally()

chin_age <- chin %>% 
  select(fish, year, fl, stock, cu, agg_name, stock_id) %>% 
  left_join(., age_in, by = "fish") %>% 
  filter(!is.na(european_age)) %>% 
  mutate(
    fw_age = ifelse(grepl("0.", european_age), "0", "1"),
    ocean_age = gsub("^.*\\.", "", european_age)
  )

labs <- chin_age %>% 
  group_by(ocean_age, fw_age) %>% 
  tally()


png(here::here("figs", "ms_figs", "age_size.png"), res = 250, units = "in", 
    height = 5.5, width = 5.5)
ggplot() +
  geom_boxplot(data = chin_age, aes(x = ocean_age, y = fl, fill = fw_age)) +
  ggsidekick::theme_sleek() +
  geom_text(data = labs, 
            aes(x = ocean_age , y = 37, label = n, group = fw_age),
            position = position_dodge(width = .7)) +
  labs(x = "Ocean Age", y = "Fork Length") +
  scale_fill_discrete(name = "Freshwater Age")
dev.off()


#age comp by CU
chin_age %>% 
  group_by(cu, agg_name) %>%
  mutate(nn = n()) %>% 
  ungroup() %>% 
  group_by(cu, agg_name, european_age, nn) %>%
  tally() %>%
  mutate(prop = n / nn) %>% 
  ggplot(., aes(x = cu, y = prop, fill = european_age)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d(name = "Age Composition", na.value = "grey60") +
  labs(y = "Proportion Age", x = "CU") +
  ggsidekick::theme_sleek() +
  facet_wrap(~agg_name, scales = "free_x")
  