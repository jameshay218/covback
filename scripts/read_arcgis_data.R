
library(tidyverse)

setwd("~/Documents/COVID-19/csse_covid_19_data/csse_covid_19_time_series/")

confirmed <- read_csv("time_series_19-covid-Confirmed.csv")
colnames(confirmed)[1:4] <- c("province","country_region","lat","long")
confirmed <- confirmed %>% pivot_longer(cols=-c("province","country_region","lat","long"), names_to="date",values_to="n")
confirmed$date <- as.POSIXct(confirmed$date,format="%m/%d/%Y")

confirmed %>% filter(country_region == "Mainland China") %>%
  ggplot() + geom_line(aes(x=date,y=n),stat="identity") + facet_wrap(~province,scales="free_y")

#confirmed %>% filter(province=="Hubei") %>% ggplot() + geom_line(aes(x=date,y=n))

confirmed <- confirmed %>% 
  group_by(province) %>%
  mutate(prev=dplyr::lag(n,n=1)) %>%
  mutate(diff=n-prev)  %>% ungroup()

write_csv(confirmed,"~/Documents/covback/data/confirmed_data.csv")

## The 13th is when the mega new cases came out
confirmed %>% filter(country_region == "Mainland China") %>% filter(date <= "20-02-20") %>%
  ggplot() + geom_line(aes(x=date,y=diff),stat="identity") + facet_wrap(~province,scales="free_y")


first_date <- as.POSIXct("20-01-01",format="%Y-%m-%d")
last_date <- as.POSIXct("20-02-20",format="%Y-%m-%d")

confirmed1 <- confirmed %>% filter(country_region == "Mainland China") %>% filter(date <= "20-02-20") %>% select(province, date, n)

incidence <- est_daily_incidence(confirmed1, first_date, last_date)
ggplot(incidence) + geom_point(aes(x=date, y=incidence)) + facet_wrap(~province,scales="free_y")

incidence <- incidence %>% filter(Date >= "20-01-23")
