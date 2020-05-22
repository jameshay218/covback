library(grid)
library(tidyverse)
library(ggpubr)
library(patchwork)

setwd("~/Documents/GitHub/covback/")
devtools::load_all()

tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")

export_probs <- read_csv("data/export_probs_matched.csv")$export_prob

## Real import probs
import_probs <- read_csv("data/import_probs_matched.csv")
import_probs <- as.matrix(import_probs[,2:ncol(import_probs)])
colnames(import_probs) <- NULL


## Simulate data
N <- 1000000
parTab <- read.csv("pars/partab_single.csv",stringsAsFactors = FALSE)

parTab[parTab$names %in% c("confirm_delay_shape","confirm_delay_scale"),"fixed"] <- 0
parTab[parTab$names == "growth_rate","values"] <- 0.15
parTab[parTab$names == "K","upper_bound"] <- N
parTab[parTab$names == "K","values"] <- N*0.6
parTab[parTab$names == "K","fixed"] <- 1
parTab[parTab$names == "reporting_propn","fixed"] <- 1
parTab[parTab$names == "reporting_propn","values"] <- 1
pars <- parTab$values
names(pars) <- parTab$names
t_cutoff <- 365
xmin <- 110

dat <- simulate_observed_data_single(pars, 200, t_cutoff, TRUE, "poisson",N,N, model_ver = "seir")
confirmed_data <- dat$aggregated %>% filter(var == "date_confirm_observable") %>% select(date, n) %>% mutate(province="1")
confirmed_data <- confirmed_data %>% mutate(n = ifelse(date <= t_cutoff, n, NA))

startTab <- generate_start_tab(parTab)
startTab[startTab$names %in% c("confirm_delay_shape","confirm_delay_scale"),"values"] <- c(6,2)
## Make strong prior on alpha and sigma

#prior_func <- NULL
## MCMC

f <- create_model_func(parTab, confirmed_data, ver="model",model_ver="seir",incubation_ver = "lnorm")
dat <- f(parTab$values)

dat <- dat %>% filter(var %in% c("infections","onsets","confirmations"))
obs <- dat %>% filter(var=="confirmations") %>% mutate(observations=rpois(n(), n)) %>% filter(date <= xmin)


p_top_a <- ggplot(obs) + 
  geom_point(aes(x=date,y=n),size=1,col="grey10")+
  ylab("New case confirmations per day") +
  xlab("Day") +
  theme_pubr() +
  scale_x_continuous(limits=c(0,110),breaks=seq(0,110,by=10)) +
  theme(axis.title=element_text(size=10),
        plot.tag=element_text(size=12),
        axis.text=element_text(size=8))+
  labs(tag="A")
x <- seq(0,30,by=0.1)

meanLog <- as.numeric(parTab[parTab$names == "lnorm_incu_par1","values"])
sdLog <- as.numeric(parTab[parTab$names == "lnorm_incu_par2","values"])
prob_sympt <- dlnorm(x, meanLog, sdLog)
incu_period <- data.frame(x=x,y=prob_sympt)
p_top_b <- ggplot(incu_period) + 
  geom_ribbon(aes(x=x,ymax=y,min=0), fill="#CC79A7") +
  ylab("Probability of\n symptom onset") +
  xlab("Days since since infection") +
  theme_pubr() +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.title=element_text(size=10),
        plot.tag=element_text(size=12),
        axis.text=element_text(size=8)) +
  labs(tag="B")

confirm_delay_shape <- 2.72 #as.numeric(parTab[parTab$names == "confirm_delay_shape","values"])
confirm_delay_scale <- 2.70 #as.numeric(parTab[parTab$names == "confirm_delay_scale","values"])
prob_confirm <- dgamma(x, shape=confirm_delay_shape, scale=confirm_delay_scale)
confirm_period <- data.frame(x=x,y=prob_confirm)

p_top_c <- ggplot(confirm_period) + 
  geom_ribbon(aes(x=x,ymax=y,min=0), fill="#56B4E9")+
  ylab("Probability of\n test confirmation") +
  xlab("Days since since symptom onset") +
  theme_pubr() +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.title=element_text(size=10),
        plot.tag=element_text(size=12),
        axis.text=element_text(size=8))+
  labs(tag="C")

colnames(dat)[3] <- "Variable"

var_key <- c("infections"="Infections",
             "onsets"="Symptom onsets",
             "confirmations"="Case confirmations")

dat$Variable <- var_key[dat$Variable]
dat$Variable <- factor(dat$Variable,levels=var_key)

p_bottom <- ggplot(dat) + 
  geom_rect(xmin=xmin,xmax=max(dat$date),ymin=0,ymax=max(dat$n),
            fill="grey70",col="grey70") +  
  geom_line(data=dat[dat$date <= xmin,], aes(x=date,y=n, col=Variable),size=1) +
  geom_line(data=dat[dat$date >= xmin,], aes(x=date,y=n, col=Variable),size=1,linetype="dashed") +
  ylab("Daily incidence") +
  xlab("Day") +
  geom_point(data=obs,aes(x=date,y=observations),size=1,col="grey10") +
  #scale_x_continuous(breaks=seq(0,80,by=10)) +
  scale_color_manual(values=c("#009E73","#0072B2", "#E69F00")) +
  theme_pubr() + 
  theme(axis.title=element_text(size=10),
        axis.text=element_text(size=8),
        plot.tag=element_text(size=12),
        legend.position=c(0.2,0.8),
        plot.margin = unit(c(0,0,0,0),"cm"))+
  labs(tag="D")
top_p <- p_top_a + (p_top_b / p_top_c) + plot_layout(widths=c(1.4,1))

top_p / p_bottom + plot_layout(heights=c(1,1.5))
pdf("figures/fig1.pdf",height=6,width=7)
top_p / p_bottom + plot_layout(heights=c(1,1.5))
dev.off()
png("figures/fig1.png",height=6,width=7,units="in",res=300)
top_p / p_bottom + plot_layout(heights=c(1,1.5))
dev.off()

