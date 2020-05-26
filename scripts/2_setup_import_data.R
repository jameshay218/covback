library(tidyverse)
setwd("~/Documents/GitHub/covback/")
devtools::load_all()

travel_probs <- read_csv("data/raw/extracted_import_proportions.csv")
province_key <- read_csv("data/raw/extracted_data_key.csv")

## Change to 4000000 for lower travel scenario
total_travellers <- 5000000
wuhan_pop_ini <- 9785388
#wuhan_pop_ini <- 11080000

## Change to 2020-01-23 for lower travel scenario
index_date_end <- as.POSIXct("2020-01-25", 
           format = "%Y-%m-%d", tz = "UTC")

tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")

## Get names that match our data
colnames(province_key)[1] <- "province"
travel_probs <- right_join(travel_probs, province_key)
travel_probs$date <- as.Date(travel_probs$date, format="%m/%d/%y")
travel_probs <- travel_probs %>% arrange(order, date)
travel_probs <- travel_probs %>% select(date, province_use, percentage)
travel_probs <- travel_probs %>% filter(date <= "2020-01-25")

## What % of travellers leave Hubei? 
travel_probs_hubei <- travel_probs %>% filter(province_use == "Hubei")
travel_probs_hubei$percentage_true <- travel_probs_hubei$percentage
travel_probs_hubei$percentage <- -1
travel_probs_hubei$prop_not_hubei <- 1 - travel_probs_hubei$percentage

## Get date range
dates <- travel_probs %>% filter(province_use == "Hubei") %>% arrange(date) %>% pull(date)

## What % of travellers that leave Wuhan go to the provinces considered here?
travel_probs <- travel_probs %>% filter(province_use != "Hubei")
travel_probs$province_use <- factor(travel_probs$province_use, levels=province_key$province_use[2:nrow(province_key)])
## Percentage is "given that someone has left Wuhan, what is the probability that they went to province X?"

ggplot(travel_probs) + geom_line(aes(x=date,y=percentage,col=province_use))

import_probs <- travel_probs %>% bind_rows(travel_probs_hubei) %>% 
  select(date, province_use, percentage) %>% 
  pivot_wider(names_from=date, values_from=percentage)
import_probs$province_use <- factor(import_probs$province_use, levels=province_key$province_use)
import_probs <- import_probs %>% arrange(province_use)

## Average first 7 days 
times_subset_early <- seq(as.Date("2019-11-01",origin="2019-11-01"), as.Date("2019-12-31",origin="2019-11-01"),by="1 day")
times_subet_late <- seq(as.Date("2020-01-26",origin="2019-11-01"), as.Date("2020-03-03",origin="2019-11-01"),by="1 day")
import_probs_first_week <- import_probs[,1:8]
import_probs_start <- import_probs_first_week %>% pivot_longer(-province_use) %>%
  mutate(name = as.Date(name, origin="2019-11-01")) %>%
  group_by(province_use) %>%
  summarise(mean_prob = mean(value)) %>%
  right_join(expand_grid(province_use=unique(import_probs$province_use),dates=times_subset)) %>%
  pivot_wider(province_use, names_from=dates,values_from=mean_prob)

import_probs_end <- expand_grid(province_use=unique(import_probs$province_use),dates=times_subet_late, mean_prob=0) %>%
  pivot_wider(province_use, names_from=dates,values_from=mean_prob)
import_probs_final <- left_join(import_probs_start,import_probs) %>% left_join(import_probs_end)


## Probability of export is the number of travellers that move out of Wuhan that
## go to other Chinese provinces
export_prob_dat <- read_csv("data/raw/export_probs_raw.csv")
export_prob_dat$Date <- as.POSIXct(export_prob_dat$Date,format="%m/%d/%Y",tz="UTC")

export_probs <- create_export_prob_matrix(total_travellers, wuhan_pop_ini, export_prob_dat, tmin, tmax,
                                          index_date_end=index_date_end)

times <- seq(tmin, tmax, by="1 day")

## Export prob is "what is the probability of someone leaving Wuhan for a non-Hubei province?
export_probs_final <- tibble(date=times,export_prob=export_probs$probs)

write_csv(export_probs_final, "data/export_probs_matched.csv")
#write_csv(export_probs_final, "data/export_probs_lower.csv")
write_csv(import_probs_final, "data/import_probs_matched.csv")