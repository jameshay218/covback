library(extrafont)
library(ggthemes)
library(ggsci)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(pals)
#install.packages("ggthemes")

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
    panel.grid.major=element_line(size=0.25,color="#f0f0f0")
  )

## Real export probs
export_probs <- read.csv("data/export_probs_matched.csv")
export_probs <- export_probs$export_prob

## Real import probs
import_probs <- read.csv("data/import_probs_matched.csv",header = TRUE)
#import_probs <- t(import_probs[,2:ncol(import_probs)])
#import_probs <- rbind(rep(1, ncol(import_probs)),import_probs)
#row.names(import_probs)[1] <- "Hubei"

tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")
export_probs_plot <- tibble(y=export_probs,x=times)
export_probs_plot <- export_probs_plot %>% filter(times >= "2019-12-08")

pA <- ggplot(data=export_probs_plot) +
  geom_line(aes(x=x,y=y))+
  geom_vline(xintercept=(as.POSIXct("2020-01-23", format="%Y-%m-%d")),
             linetype="dashed") +
  scale_x_datetime(breaks="5 days")+
  theme_bw() +
  scale_y_continuous(limits=c(0,0.02), breaks=seq(0,0.02,by=0.002)) +
  ylab("Daily probability of leaving Wuhan\n to another Chinese province") +
  xlab("Date") +
  export_theme

import_probs_melted <- reshape2::melt(import_probs,id="province_use")
import_probs_melted$variable <- as.numeric(as.factor(import_probs_melted$variable))
import_probs_melted$variable <- times[import_probs_melted$variable]
import_probs_melted <- import_probs_melted %>% 
  filter(variable >= "2019-12-08" & province_use != "Hubei")
colnames(import_probs_melted)[1] <- "Province"

import_order <- import_probs_melted %>% group_by(Province) %>% summarise(total_prob=sum(value)) %>%
  arrange(-total_prob) %>% pull(Province)
import_probs_melted$Province <- factor(import_probs_melted$Province, levels=import_order)


pB <- ggplot(data=import_probs_melted) +
  geom_line(aes(x=variable,y=value, col=Province))+
  geom_vline(xintercept=(as.POSIXct("2020-01-23", format="%Y-%m-%d")),
             linetype="dashed") +
  #scale_colour_manual(values=as.vector(polychrome(29)))+
  scale_x_datetime(breaks="5 days")+
  scale_y_continuous(limits=c(0,0.25),breaks=seq(0,0.25,by=0.05)) +
  guides(col=guide_legend(ncol=3)) +
  theme_bw() +
  ylab("Proportion of Wuhan travelers received") +
  xlab("Date") +
  export_theme + 
  theme(legend.position=c(0.8,0.5))
pB

figS3 <- pA / pB


pdf("figS5.pdf",height=8,width=8)
figS3
dev.off()
png("figS5.png",height=8,width=8,units="in",res=300)
figS3
dev.off()

