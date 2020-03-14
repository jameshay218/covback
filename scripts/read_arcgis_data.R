
library(tidyverse)

setwd("~/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series/")

confirmed <- read_csv("time_series_19-covid-Confirmed.csv")
colnames(confirmed)[1:4] <- c("province","country_region","lat","long")
confirmed <- confirmed %>% pivot_longer(cols=-c("province","country_region","lat","long"), names_to="date",values_to="n")
confirmed$date <- as.POSIXct(confirmed$date,format="%m/%d/%y", origin="01/01/2019",tz="UTC")

confirmed %>% filter(country_region == "Mainland China") %>%
  ggplot() + geom_line(aes(x=date,y=n),stat="identity") + facet_wrap(~province,scales="free_y")

#confirmed %>% filter(province=="Hubei") %>% ggplot() + geom_line(aes(x=date,y=n))

confirmed <- confirmed %>% 
  group_by(province) %>%
  mutate(prev=dplyr::lag(n,n=1)) %>%
  mutate(diff=n-prev)  %>% ungroup()

## Remove Shandong after 20th Feb as big spike
confirmed <- confirmed %>% mutate(diff = ifelse(province == "Shandong" & date > "2020-02-20",NA, diff))

## First day is just starting confirmed case number
#confirmed <- confirmed %>% mutate(diff=ifelse(is.na(diff), n, diff))
#confirmed <- confirmed %>% mutate(diff=ifelse(diff < 0, 0, diff))
write_csv(confirmed,"~/Documents/GitHub/covback/data/confirmed_data.csv")

## The 13th is when the mega new cases came out
confirmed %>% filter(country_region == "Mainland China") %>% #filter(date <= "20-02-20") %>%
  ggplot() + geom_line(aes(x=date,y=diff),stat="identity") + facet_wrap(~province,scales="free_y")


first_date <- as.POSIXct("2020-01-01",format="%Y-%m-%d")
last_date <- as.POSIXct("2020-02-20",format="%Y-%m-%d")

confirmed1 <- confirmed %>% filter(country_region == "Mainland China") %>% filter(date <= "20-02-20") %>% select(province, date, n)

incidence <- est_daily_incidence(confirmed1, first_date, last_date)
ggplot(incidence) + geom_point(aes(x=date, y=incidence)) + facet_wrap(~province,scales="free_y")

incidence <- incidence %>% filter(Date >= "20-01-23")
