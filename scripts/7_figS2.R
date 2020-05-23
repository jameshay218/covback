library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggthemes)

export_theme <- theme_tufte() + 
  theme(
    ## Axis text and titles
    axis.text.x = element_text(size=8,angle = 45, hjust = 1,family="sans"),
    axis.text.y=element_text(size=8,family="sans"),
    axis.title.x=element_text(size=10,family="sans",vjust=-1),
    axis.title.y=element_text(size=10,family="sans"),
    ## Axis lines
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),
    ## Legends
    legend.title=element_text(size=10,family="sans",face="italic"),
    legend.text=element_text(size=10,family="sans"),
    legend.key.size= unit(0.2, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    ## Strips for facet_wrap
    strip.text=element_text(size=8,family="sans"),
    strip.background=element_rect(fill="#f0f0f0"),
    ## Overall plot
    panel.grid.major=element_line(size=0.25,color="#f0f0f0"),
    plot.tag=element_text(size=12,family="sans",face="bold")
    
  )
cols <- as.vector(pals::polychrome(11))

setwd("~/Google Drive/nCoV/backcalculation_paper/figures_final_raw//")
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
  all_incidences[[tmp]] <- inc
}

all_prevalences <- do.call("rbind", all_prevs)
all_prevalences <- all_prevalences %>% mutate(var=ifelse(var == "Prevalence of infections eligible\n for international travel", "Prevalence of infections eligible for international travel", var))
all_prevalences %>% 
  filter(var == "Prevalence of infections eligible for international travel") %>% 
  group_by(province, runname) %>% 
  filter(median==max(median)) %>% #View()
  ungroup() %>% group_by(runname) %>% summarise(max_date=max(date),min_date=min(date))
#select(date)


n_wuhan <- 9785388

runname_key <- c(
  "main"="Model variant 1 (baseline, intermediate Wuhan prevalence)",
  "diffuse_r_local"="Model variant 2 (increased local transmission, lower Wuhan prevalence)", 
  "t_switch_6"="Model variant 3 (symptom peak 29th Jan, higher Wuhan prevalence)",
                 "earliest_seed"="Model variant 4 (17th Nov seed)",
                 "early_seed"="Model variant 5 (1st Dec seed)",
                 "fewer_travellers"="Model variant 6 (4 million travellers)",
  "t_switch_1"="Model variant 7 (symptom peak 23rd Jan)",
  "t_switch_2"="Model variant 8 (symptom peak 24th Jan)",
  "t_switch_3"="Model variant 9 (symptom peak 26th Jan)",
  "t_switch_4"="Model variant 10 (symptom peak 27th Jan)",
  "t_switch_5"="Model variant 11 (symptom peak 28th Jan)",
  "fixed_serial_int"="Model variant 12 (fixed serial interval distribution)"
  
  )

all_prevs_comb <- do.call("rbind",all_prevs_hubei)
all_prevs_comb <- all_prevs_comb %>% mutate(var=ifelse(var == "Prevalence of infections eligible\n for international travel", "Prevalence of infections eligible for international travel in Wuhan", var))
all_prevs_comb <- all_prevs_comb %>% filter(runname != "r_local_unconstrained")

all_prevs_comb$runname1 <- runname_key[all_prevs_comb$runname]

runname_levels <- all_prevs_comb %>% 
  filter(var == "Prevalence of infections eligible for international travel in Wuhan") %>% 
  group_by(runname1) %>% 
  filter(median == max(median)) %>% 
  arrange(-median) %>% pull(runname1)


all_prevs_comb$runname1 <- factor(all_prevs_comb$runname1, levels=runname_levels)
colnames(all_prevs_comb)[13] <- "Sensitivity"
all_prevs_plot <- all_prevs_comb %>% filter(var %in% c("Prevalence of infections eligible for international travel in Wuhan")) %>% ggplot() + 
  geom_ribbon(aes(x=date,ymin=lower/n_wuhan,ymax=upper/n_wuhan,fill=Sensitivity),alpha=0.5) +
  geom_line(aes(x=date,y=median/n_wuhan,col=Sensitivity)) +
  scale_y_continuous(limits=c(0,0.02),breaks=seq(0,0.02,by=0.0025)) +
  scale_x_date(breaks="7 days",limits=c(as.Date("2019-12-08",origin="2019-11-01"),as.Date("2020-01-30",origin="2019-11-01"))) +
  scale_color_manual(values=c(cols,"black"))+
  scale_fill_manual(values=c(cols,"black"))+
  export_theme + facet_wrap(~var) + 
  theme(legend.position=c(0.4,0.7),strip.text=element_text(size=10))+
  geom_vline(xintercept=as.Date("2020-01-23"),linetype="dashed") +
ylab("Daily prevalence (per capita)") +
  xlab("Date")


pdf("~/Documents/all_prev.pdf",width=7,height = 5)
all_prevs_plot
dev.off()
png("~/Documents/all_prev.png",width=7,height = 5,res=300,units="in")
all_prevs_plot
dev.off()

