setwd("~/Documents/GitHub/covback/")

## Fig 1:  incidence in China
## Case confirmation data
tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")

confirmed_data <- read_csv("data/real/confirmed_data_matched_final.csv")
confirmed_data$time <- times[confirmed_data$date + 1]
confirmed_data <- confirmed_data %>% select(province, date, n)
colnames(confirmed_data)[3] <- "n"
all_reports <- expand.grid(province=unique(confirmed_data$province),
                           date=match(times, times))
confirmed_data <- confirmed_data %>% right_join(all_reports)
confirmed_data$date <- confirmed_data$date - 1

confirmed_data %>% 
  ggplot() + geom_line(aes(x=date,y=n)) + facet_wrap(~province,scales="free_y")

## Real import probs
import_probs <- read.csv("data/import_probs_matched.csv",header = TRUE)
import_probs <- t(import_probs[,2:ncol(import_probs)])
import_probs <- rbind(rep(1, ncol(import_probs)),import_probs)
row.names(import_probs)[1] <- "Hubei"

use_provinces <- row.names(import_probs)
use_provinces <- intersect(use_provinces, unique(confirmed_data$province))

#confirmed_data <- confirmed_data %>% filter(province %in% use_provinces)
#import_probs <- import_probs[match(use_provinces, row.names(import_probs)),]
#confirmed_data$province <- match(confirmed_data$province, use_provinces)
#confirmed_data <- confirmed_data %>% arrange(province, date)


factor_levels <- confirmed_data %>% filter(province != "Hubei") %>% 
  group_by(province) %>%
  summarise(x=sum(n,na.rm=TRUE)) %>%
  arrange(-x) %>% pull(province)

quants <- quants[!(quants$province == "Hubei" & quants$var == "infections"),]
quants1 <- quants[!(quants$province == "Hubei" & quants$date > "2020-01-23"),]
quants1 <- quants1 %>% filter(province != "Hubei")
confirmed_data2 <- confirmed_data %>% filter(province != "Hubei")

confirmed_data2$province <- factor(confirmed_data2$province, levels=factor_levels)
quants1$province <- factor(quants1$province, levels=factor_levels)

confirmed_data2 <- confirmed_data1
confirmed_data2$date <- as.Date(confirmed_data2$date, origin="2019-11-01")
confirmed_data2$province <- as.numeric(confirmed_data2$province)
confirmed_data2$province <- provinces[confirmed_data2$province]


quants1 <- quants_summary

confirmed_data2$province <- factor(confirmed_data2$province,levels=factor_levels_names)
confirmed_data2 <- confirmed_data2 %>% drop_na()
quants1$province <- factor(quants1$province,levels=factor_levels_names)
quants1 <- quants1 %>% drop_na()

p <- ggplot(quants1[quants1$var %in% c("infections","observations","onsets") & quants1$date >= "2020-01-01",]) +
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
  geom_line(aes(x=date,y=median,col=var)) +
  geom_point(data=confirmed_data2[confirmed_data2$date >= "2020-01-01",],aes(x=date,y=n),size=0.5) +
  scale_x_date(breaks="7 days") +
  geom_vline(xintercept=as.Date("2020-01-23",origin="2019-11-01"),linetype="dashed") +
  theme_pubr() + 
  theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~province,ncol=5,scales="free_y")
pdf("incidence_all.pdf",height=7,width=8)
p
dev.off()
