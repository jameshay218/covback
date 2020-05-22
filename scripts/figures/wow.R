library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggthemes)

export_theme <- theme_tufte() + 
  theme(
    ## Axis text and titles
    axis.text.x = element_text(size=8,angle = 45, hjust = 1,family="sans"),
    axis.text.y=element_text(size=8,family="sans"),
    axis.title.x=element_text(size=8,family="sans",vjust=-1),
    axis.title.y=element_text(size=8,family="sans"),
    ## Axis lines
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),
    ## Legends
    legend.title=element_text(size=8,family="sans",face="italic"),
    legend.text=element_text(size=8,family="sans"),
    legend.key.size= unit(0.2, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    ## Strips for facet_wrap
    strip.text=element_text(size=8,family="sans",face="bold"),
    strip.background=element_blank(),
    ## Overall plot
    panel.grid.major=element_line(size=0.25,color="#f0f0f0"),
    plot.tag=element_text(size=12,family="sans",face="bold")
    
  )
cols <- ggsci::pal_npg()(10)

setwd("~/Google Drive/nCoV/backcalculation_paper/figures_check/")
dirs <- list.files()
topwd <- getwd()
all_prevs <- list()
all_prevs_hubei <- list()
all_incidences <- list()
for(tmp in dirs){
  print(tmp)
  setwd(topwd)
  setwd(tmp)
  
  prev_file <- list.files(pattern="hubei_prevalence.csv")
  prev <- read_csv(prev_file)
  all_prevs_hubei[[tmp]] <- prev
  
}

#dirs <- c("diffuse_r_local","main","t_switch_6")
for(tmp in dirs){
  print(tmp)
  setwd(topwd)
  setwd(tmp)
  
  prev_file <- list.files(pattern="all_prevalence.csv")
  prev <- read_csv(prev_file)
  prev$runname <- tmp
  all_prevs[[tmp]] <- prev
  
  inc_file <- list.files(pattern="all_incidence.csv")
  inc <- read_csv(inc_file)
  inc$version <- tmp
  all_incidences[[tmp]] <- inc
}

## Get range of peak times by province
all_incidences_comb <- do.call("rbind", all_incidences)
range_province_peaks <- all_incidences_comb %>% 
  filter(province %in% c("Beijing","Hunan","Shanghai","Wuhan") & var == "Infection incidence") %>% 
  group_by(version, province) %>%
  filter(median == max(median)) %>% group_by(province) %>% 
  summarise(early_date=min(date),late_date=max(date))
range_ylim <- tibble(province=c("Hunan","Beijing","Shanghai"),ymin=-1,ymax=c(90,38,30))
range_province_peaks <- range_province_peaks %>% left_join(range_ylim)

dirs <- c("diffuse_r_local","main","t_switch_6")
all_incidences <- all_incidences[dirs]

all_prevalences <- do.call("rbind", all_prevs)
all_prevalences %>% 
  filter(var == "Prevalence of infections eligible\n for international travel") %>% 
  group_by(province, runname) %>% 
  filter(median==max(median)) %>% #View()
  ungroup() %>% group_by(runname) %>% summarise(max_date=max(date),min_date=min(date))
  #select(date)

n_wuhan <- 11080000
all_prevs_comb <- do.call("rbind",all_prevs_hubei)
all_prevs_comb %>% filter(var %in% c("Prevalence of infections eligible\n for international travel")) %>% ggplot() + 
  geom_ribbon(aes(x=date,ymin=lower/n_wuhan,ymax=upper/n_wuhan,fill=runname),alpha=0.5) +
  geom_line(aes(x=date,y=median/n_wuhan,col=runname)) +
  scale_y_continuous(limits=c(0,0.03)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_pubr() + facet_wrap(~var) + 
  geom_vline(xintercept=as.Date("2020-01-23"),linetype="dashed") +
  geom_vline(xintercept=as.Date("2020-01-31"),linetype="dashed")


pop_wuhan <- 11080000

tmp_hubei <- all_prevs_comb %>% filter(date >= "2019-12-08" & 
                                         var %in% c("Prevalence of infections eligible\n for international travel") & 
                                         runname %in% c("diffuse_r_local","main","t_switch_6"))
runname_key <- c("diffuse_r_local"="Scenario 2", "t_switch_6"="Scenario 3", "main"="Scenarios 1, 4 and 5")
tmp_hubei$runname <- runname_key[tmp_hubei$runname]
tmp_hubei$runname <- factor(tmp_hubei$runname, levels=c("Scenarios 1, 4 and 5", "Scenario 2", "Scenario 3"))
colnames(tmp_hubei)[which(colnames(tmp_hubei) == "runname")] <- "Scenario"



tmp_hubei %>% group_by(Scenario) %>% filter(date == "2020-01-23") %>%
  mutate(median = median/pop_wuhan,
         lower=lower/pop_wuhan,
         upper=upper/pop_wuhan) %>% select(Scenario, median, lower, upper)
tmp_hubei$province <- "Wuhan"


hubei_plot <- tmp_hubei %>%
  ggplot() + 
  geom_rect(xmin=as.Date("2020-01-23"),xmax=as.Date("2021-01-23"),ymin=0,ymax=1,fill="grey90",alpha=0.25)+
  geom_vline(xintercept=(as.Date("2020-01-23", format="%Y-%m-%d", tz="LMT", origin="2019-11-01")),
             linetype="dashed") +
  geom_line(aes(x=date,y=median/pop_wuhan,col=Scenario)) + 
  geom_ribbon(aes(x=date,ymin=lower/pop_wuhan,ymax=upper/pop_wuhan,fill=Scenario),alpha=0.25) +
  scale_fill_manual(values=c("black","#4DBBD5FF","#7E6148FF"))+
  scale_color_manual(values=c("black","#4DBBD5FF","#7E6148FF"))+
  scale_x_date(breaks=seq(as.Date("2020-01-01"),as.Date("2020-01-31"),by="5 days"),limits=c(as.Date("2020-01-01"),as.Date("2020-01-31")),labels=seq(as.Date("2020-01-01"),as.Date("2020-01-31"),by="5 days")) + 
  scale_y_continuous(expand=c(0,0),limits=c(0,0.02),breaks=seq(0,0.02,by=0.005)) +
  theme_bw() +
  ylab("Daily prevalence (per capita)") +
  xlab("Date") +
  export_theme + 
  theme(legend.position=c(0.2,0.8),
        plot.margin = margin(0,0,0,0, "cm"))+
  facet_wrap(~province) +
  labs(tag="B")
hubei_plot


all_incidences_comb_overall <- do.call("rbind", all_incidences)
all_incidences_comb <- all_incidences[["main"]]
all_incidences_comb <- all_incidences_comb %>% filter(province %in% c("Hunan","Beijing","Shanghai")) 


import_probs <- read.csv("~/Documents/GitHub/covback/data/import_probs_matched.csv")
provinces <- as.character(import_probs[,1])
provinces[which(provinces == "Inner Mongolia")] <- "Inner_Mongolia"


confirmed_data1 <- as.data.frame(read_csv("~/Documents/GitHub/covback/data/real/midas_data_final.csv"))
confirmed_data1 <- confirmed_data1 %>% mutate(n=ifelse(province==1, NA, n))
confirmed_data1 <- confirmed_data1 %>% select(-province_raw)


factor_levels <- confirmed_data1  %>% filter(province != 1) %>% 
  group_by(province) %>%
  summarise(x=sum(n,na.rm=TRUE)) %>%
  arrange(-x) %>% pull(province)
factor_levels_names <- provinces[factor_levels]

confirmed_data2 <- confirmed_data1
confirmed_data2$date <- as.Date(confirmed_data2$date, origin="2019-11-01")
confirmed_data2$province <- as.numeric(confirmed_data2$province)
confirmed_data2$province <- provinces[confirmed_data2$province]
confirmed_data2$province <- factor(confirmed_data2$province,levels=factor_levels_names)
confirmed_data2 <- confirmed_data2 %>% drop_na()
confirmed_data2 <- confirmed_data2 %>% filter(province %in% c("Hunan","Beijing","Shanghai"))

all_incidences_comb$var  <- factor(all_incidences_comb$var, levels=c("Infection incidence", "Symptom onset incidence","Simulated incidence of case confirmations"))

p1 <- ggplot(all_incidences_comb[all_incidences_comb$date >= as.Date("2020-01-01"),]) +
  geom_rect(data=range_province_peaks,aes(xmin=early_date,xmax=late_date,ymin=ymin,ymax=ymax),fill="orange",alpha=0.25) +
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
  geom_line(aes(x=date,y=median,col=var)) +
  geom_point(data=confirmed_data2[confirmed_data2$date >= "2020-01-01",],aes(x=date,y=n),size=0.5) +
  scale_x_date(breaks="7 days") +
  scale_y_continuous(breaks=seq(0,100,by=10),expand=c(0,0)) +
  scale_fill_manual(values=cols[c(1,3,4)])+
  scale_color_manual(values=cols[c(1,3,4)])+
  geom_vline(xintercept=as.Date("2020-01-23",origin="2019-11-01"),linetype="dashed") +
  theme_bw() +
  xlab("") +
  ylab("Daily incidence (absolute numbers)") +
  facet_wrap(~province,scales="free_y") +
  export_theme +
  theme(
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position="top",
        legend.margin=margin(0,0,0,0),
        plot.margin = margin(0,0,-0.5,0, "cm"),
        legend.box.margin=margin(-10,-10,-5,-10)
  )+ 
  labs(tag="A")
pdf("~/Documents/Fig1.pdf",width=7,height = 6)
p1 /hubei_plot + plot_layout(heights=c(1,1.5))
dev.off()
png("~/Documents/Fig1.png",width=7,height = 6,res=300,units="in")
p1 /hubei_plot + plot_layout(heights=c(1,1.5))
dev.off()
