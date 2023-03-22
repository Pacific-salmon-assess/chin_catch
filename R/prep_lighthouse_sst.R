sst <- read.csv(
  here::here(
    "data", "spatial_data", 
    "Amphitrite_Point_-_Average_Monthly_Sea_Surface_Temperatures_1934-2022.csv"
  )
) %>% 
  janitor::clean_names() %>% 
  pivot_longer(., cols = c(jan:dec), names_to = "month", values_to = "sst") %>% 
  mutate(month =  match(month, tolower(month.abb))) %>% 
  filter(year >= 2019,
         month > 3 & month < 10) %>% 
  glimpse()

saveRDS(sst, 
        here::here("data", "spatial_data", "trim_amphitrite.rds"))

ggplot(sst) +
  geom_line(aes(x = month_n, y = sst, colour = as.factor(year)))
