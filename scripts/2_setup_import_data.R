library(tidyverse)
setwd("~/Documents/GitHub/covback/")
devtools::load_all()
travel_probs <- read_csv("data/raw/extracted_import_proportions.csv")
province_key <- read_csv("data/raw/extracted_data_key.csv")

total_travellers <- 5000000
wuhan_pop_ini <- 11080000

tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")

## Get names that match our data
colnames(province_key)[1] <- "province"
travel_probs <- right_join(travel_probs, province_key)
travel_probs$date <- as.Date(travel_probs$date, format="%m/%d/%y")
travel_probs <- travel_probs %>% arrange(order, date)
travel_probs <- travel_probs %>% select(date, province_use, percentage)
travel_probs <- travel_probs %>% filter(date <= "2020-01-25")

export_probs <- travel_probs %>% 
  filter(province_use == "Hubei") %>% 
  arrange(date) %>% mutate(x=1-percentage)

## First deal with imports
## Get date range
dates <- travel_probs %>% filter(province_use == "Hubei") %>% arrange(date) %>% pull(date)

## Get proportion of people leaving Wuhan that go to provinces other than Hubei
travel_probs_hubei <- travel_probs %>% filter(province_use == "Hubei")
travel_probs_hubei$percentage_scaled <- 1
travel_probs <- travel_probs %>% filter(province_use != "Hubei")
export_frac <- travel_probs %>% group_by(date) %>% summarise(export_prop=sum(percentage))
travel_probs <- travel_probs %>% right_join(export_frac)

## Convert to relative proportion of people going to other Chinese provinces that go to each province
travel_probs$percentage_scaled <- travel_probs$percentage/travel_probs$export_prop
travel_probs$province_use <- factor(travel_probs$province_use, levels=province_key$province_use[2:nrow(province_key)])

ggplot(travel_probs) + geom_line(aes(x=date,y=percentage_scaled,col=province_use))


import_probs <- travel_probs %>% bind_rows(travel_probs_hubei) %>% 
  select(date, province_use, percentage_scaled) %>% 
  pivot_wider(names_from=date, values_from=percentage_scaled)

## Probability of export is the number of travellers that move out of Wuhan that
## go to other Chinese provinces
export_prob_dat <- read_csv("data/raw/export_probs_raw.csv")
#export_prob_dat$frac_leave_hubei <- export_prob_dat$x
export_prob_dat$Date <- as.POSIXct(export_prob_dat$Date,format="%m/%d/%Y",tz="UTC")
export_probs <- create_export_prob_matrix(total_travellers, wuhan_pop_ini, export_prob_dat, tmin, tmax,
                                          index_date_end=as.POSIXct("2020-01-23", 
                                                                    format = "%Y-%m-%d", tz = "UTC"))

times <- seq(tmin, tmax, by="1 day")
export_probs_final <- tibble(date=times,export_prob=export_probs$probs)

write_csv(export_probs_final, "export_probs_lower.csv")
write_csv(import_probs, "import_probs.csv")