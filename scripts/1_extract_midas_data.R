library(tidyverse)
setwd("~/Documents/GitHub/midas_covid/data/cases/china/daily_cases_chinacdc_EN/")

## List all files with data
files <- list.files()
files <- files[!(files %in% c("collection_metadata.csv","data_guide.csv"))]

all_dat <- NULL
all_dats <- NULL
## Only want to store province name and number of new diagnoses
subset_cols <- c("province","New diagnosis")

## REad in each file
for(file in files){
  tmp_dat <- read_csv(file)
  
  ## Case sensitive
  colnames(tmp_dat)[which(colnames(tmp_dat) == "Province")] <- "province"  
  tmp_dat <- tmp_dat[,subset_cols]
  
  ## Get date from filename
  tmp_dat$date <- str_split(file, pattern="_")[[1]][1]
  all_dats[[file]] <- tmp_dat
  all_dat <- bind_rows(all_dat, tmp_dat)
}
all_dat$date <- as.Date(all_dat$date)

## Need to match names given here to names that I'm using
name_keys <- read_csv("~/Documents/GitHub/covback/data/real/province_names.csv")

## Merge the data with name keys to get the name I need
colnames(all_dat) <- c("province_raw","n","date")
all_dat <- left_join(all_dat, name_keys)

## Check against arcgis data, see generally the same but some small numerical differences
confirmed_dat1 <- read_csv("~/Documents/GitHub/covback/data/real/confirmed_data_matched_final.csv")
confirmed_dat1$date <- as.Date(confirmed_dat1$date, origin="2019-11-01")
confirmed_dat1 <- confirmed_dat1 %>% select(province_name, date,n)
colnames(confirmed_dat1)[1] <- "province"

#all_dat %>% ggplot() + geom_line(aes(x=date, y=n)) + 
#  geom_point(data=confirmed_dat1,aes(x=date,y=n),col="red") + facet_wrap(~province,scales="free_y")

## Only want province I'm using
all_dat <- all_dat %>% filter(all_dat$province %in% unique(confirmed_dat1$province))

## Change factor levels and reorder by total number of cases
all_dat$province <- factor(all_dat$province, levels=unique(confirmed_dat1$province))
all_dat <- all_dat %>% select(date, province, n)
colnames(all_dat)[2] <- "province_raw"
## Remove Shandong after 19th Feb as big spike
all_dat <- all_dat %>% mutate(n = ifelse(province_raw == "Shandong" & date > "2020-02-19",NA, n))

## Cutoff of 3rd march, as want before second increase
all_dat <- all_dat %>% filter(date <= "2020-03-03")

## Get data from 1st November 2019 to now, filling with zeros
## Also the 15th is missing, so add zeros for that day. Pretty much fine as way before majority of cases
## For most provinces
times_fill <- c(seq(as.Date("2019-11-01"),as.Date("2020-01-09"),by="1 day"),as.Date("2020-01-15"))
all_dat_fill <- as_tibble(expand.grid(date=times_fill,province_raw=unique(all_dat$province_raw),n=0))
all_dat_fill$province_raw <- factor(all_dat_fill$province_raw, levels=unique(confirmed_dat1$province))


## Merge with dummy data, arrange and visualise
all_dat_final <- bind_rows(all_dat, all_dat_fill)
all_dat_final$province <- as.numeric(all_dat_final$province_raw)
all_dat_final <- all_dat_final  %>% arrange(province, date)
all_dat_final %>% ggplot() + geom_line(aes(x=date, y=n)) + facet_wrap(~province_raw, scales="free_y")

times <- seq(as.Date("2019-11-01"),as.Date("2020-03-03"), by="1 day")

all_dat_final$date <- match(all_dat_final$date, times) - 1

## Save
write_csv(all_dat_final, "~/Documents/GitHub/covback/data/real/midas_data_final.csv")



