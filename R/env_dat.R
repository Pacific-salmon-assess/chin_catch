## Sea surface temperature 
# 1) Amphitrite lighthouse
# 2) La Perouse weather buoy


library(tidyverse)


amph_dat <- read.csv(
  here::here(
    "data", "amphitrite_daily_sst.csv"
  ),
  skip = 1
) %>% 
  janitor::clean_names() %>% 
  mutate(
    date = as.POSIXct(date_yyyy_mm_dd),
    year = lubridate::year(date),
    month = lubridate::month(date)
  ) 

amph_dat %>%
  filter(year > 2018, 
         month > 4 & month < 10, 
         !temperature_c > 100) %>% 
  mutate(
    month = as.factor(month),
    year = as.factor(year)
  ) %>% 
  ggplot(.) +
  geom_boxplot(aes(x = month, y = temperature_c, fill = year))


buoy_dat <- read.csv(
  here::here(
    "data", "la_perouse_buoy.csv"
  )
) %>% 
  janitor::clean_names() %>% 
  mutate(
    date_time = as.POSIXct(date, format = "%m/%d/%Y %H:%M"),
    year = lubridate::year(date_time),
    month = lubridate::month(date_time)
  ) 

buoy_dat_trim <- buoy_dat %>% 
  filter(
    year > 2018, 
    month > 4 & month < 10
  ) %>% 
  mutate(
    month = as.factor(month),
    year = as.factor(year)
  ) 

ggplot(buoy_dat_trim) +
  geom_boxplot(aes(x = month, y = sstp, fill = year))


buoy_dat_trim %>% 
  group_by(year) %>% 
  summarize(sd(sstp))

buoy_dat_trim %>% 
  group_by(month) %>% 
  summarize(sd(sstp))
