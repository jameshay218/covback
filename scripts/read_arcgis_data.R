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

## First day is just starting confirmed case number
#confirmed <- confirmed %>% mutate(diff=ifelse(is.na(diff), n, diff))
#confirmed <- confirmed %>% mutate(diff=ifelse(diff < 0, 0, diff))

## Get factor order
confirmed <- confirmed %>% filter(country_region=="Mainland China")
factor_levels <- confirmed %>% group_by(province) %>% 
  filter(n == max(n,na.rm=TRUE)) %>% select(province, n) %>%
  unique() %>%
  arrange(-n) %>% pull(province)

## Rearrange province key
province_key <- read_csv("~/Documents/GitHub/covback/data/raw/extracted_data_key.csv")
province_key$order <- match(province_key$province_use, factor_levels)
province_key <- province_key %>% arrange(order)
write_csv(province_key, "~/Documents/GitHub/covback/data/raw/extracted_data_key.csv")

confirmed$province <- factor(confirmed$province, levels=province_key$province_use)
confirmed <- confirmed %>% arrange(province) %>% drop_na()

## Remove Shandong after 20th Feb as big spike
confirmed1 <- confirmed %>% mutate(diff = ifelse(province == "Shandong" & date > "2020-02-20",NA, diff))

write_csv(confirmed1,"~/Documents/GitHub/covback/data/real/confirmed_data_raw.csv")

confirmed_matched <- confirmed %>% select(province, date, diff) %>% arrange(province) %>% drop_na()
confirmed_matched <- confirmed_matched %>% mutate(diff = ifelse(province == "Shandong" & date > "2020-02-20",NA, diff))

tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")

confirmed_matched$date <- match(confirmed_matched$date, times) - 1
extra_dats <- expand.grid(province = unique(confirmed_matched$province),
                          date=0:(min(confirmed_matched$date)-1),
                          diff=NA)
confirmed_matched <- bind_rows(confirmed_matched, extra_dats)
confirmed_matched <- confirmed_matched %>% arrange(province, date)
colnames(confirmed_matched) <- c("province","date","n")
confirmed_matched$province_name <- confirmed_matched$province
confirmed_matched$province <- as.numeric(confirmed_matched$province)
## Remove Shandong after 20th Feb as big spike
write_csv(confirmed_matched, "~/Documents/GitHub/covback/data/real/confirmed_data_matched_final.csv")

## The 13th is when the mega new cases came out
confirmed %>% filter(country_region == "Mainland China") %>% #filter(date <= "20-02-20") %>%
  ggplot() + geom_line(aes(x=date,y=diff),stat="identity") + facet_wrap(~province,scales="free_y")

first_date <- as.POSIXct("2020-01-01",format="%Y-%m-%d")
last_date <- as.POSIXct("2020-02-20",format="%Y-%m-%d")

confirmed1 <- confirmed %>% filter(country_region == "Mainland China") %>% filter(date <= "20-02-20") %>% select(province, date, n)

incidence <- est_daily_incidence(confirmed1, first_date, last_date)
ggplot(incidence) + geom_point(aes(x=date, y=incidence)) + facet_wrap(~province,scales="free_y")

incidence <- incidence %>% filter(Date >= "20-01-23")
