library(grid)
library(lazymcmc)
library(tidyverse)
library(ggpubr)
library(patchwork)

setwd("~/Documents/GitHub/covback")
Rcpp::compileAttributes()
devtools::document()
devtools::load_all()
#set.seed(2)

tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")

confirmed_data <- read_csv("data/confirmed_data.csv")
confirmed_data <- confirmed_data %>% filter(country_region == "Mainland China")
confirmed_data$date <- match(confirmed_data$date, times)
confirmed_data <- confirmed_data %>% select(province, date, diff)
colnames(confirmed_data)[3] <- "n"
all_reports <- expand.grid(province=unique(confirmed_data$province),
                           date=match(times, times))
confirmed_data <- confirmed_data %>% right_join(all_reports)
confirmed_data$date <- confirmed_data$date - 1
#confirmed_data <- confirmed_data %>% filter(date <= 80)

confirmed_data %>% #filter((province != "Hubei") | (province == "Hubei" & date <= 86)) %>%
ggplot() + geom_line(aes(x=date,y=n)) + facet_wrap(~province,scales="free_y")

## Incubation period draws and parameter values
inc_period_draws <- read.csv("data/backer_draws.csv",stringsAsFactors=FALSE)
inc_period_draws1 <- read.csv("data/backer_weibull_draws.csv")
parTab <- read.csv("pars/partab_provinces.csv",stringsAsFactors=FALSE)

## Serial interval draws
#serial_interval_draws <- read.csv("~/Documents/Github/covback/data/lognormal-truncated.csv")
#wow <- fit_lnorm_normal_prior(serial_interval_draws)
## Use estimated serial interval from http://rs.yiigle.com/yufabiao/1183269.htm
## gamma with pars 5.23 and 0.87
serial_interval_par1 <- 5.23
serial_interval_par2 <- 0.87
serial_interval <- dgamma(0:20,serial_interval_par1,serial_interval_par2)

plot(calculate_serial_interval_probs(40, serial_interval_par1, serial_interval_par2),type='l')


## Real export probs
export_probs <- read.csv("data/export_probs.csv")
#export_probs <- export_probs[1:112,]
#tmax <- nrow(export_probs)-1
export_probs <- export_probs$prob_leaving

## Real import probs
import_probs <- read.csv("data/import_probs.csv",header = TRUE)
#import_probs <- import_probs[1:112,]
import_probs <- t(import_probs[,2:ncol(import_probs)])
import_probs <- rbind(rep(1, ncol(import_probs)),import_probs)
row.names(import_probs)[1] <- "Hubei"

use_provinces <- row.names(import_probs)
use_provinces <- intersect(use_provinces, unique(confirmed_data$province))

confirmed_data <- confirmed_data %>% filter(province %in% use_provinces)
import_probs <- import_probs[match(use_provinces, row.names(import_probs)),]
confirmed_data$province <- match(confirmed_data$province, use_provinces)
confirmed_data <- confirmed_data %>% arrange(province, date)

#times <- 0:tmax
parTab <- create_many_province_partab(parTab, length(use_provinces), tmax)
provinces <- unique(parTab$province)
provinces <- provinces[provinces != "all"]
n_provinces <- length(provinces)

## Make strong prior on alpha and sigma
#prior_func <- create_incubation_prior(inc_period_draws)
prior_func <- create_prior_startdate(parTab, inc_period_draws, 37, 5)
## Generate some fake data

confirm_delay_pars <- read.csv("data/fitted_confirm_delays.csv")
confirm_delay_pars <- confirm_delay_pars %>%  select(shape, scale)
confirm_delay_pars$date_onset <- 0:(nrow(confirm_delay_pars)-1)
confirm_delay_pars <- confirm_delay_pars %>% filter(date_onset <= tmax)
plot_reporting_landscape(confirm_delay_pars$shape, confirm_delay_pars$scale)

parTab[parTab$names == "t0" & parTab$province == "1","fixed"] <- 1
parTab[parTab$names == "t0" & parTab$province == "1","values"] <- 37
parTab[parTab$names == "growth_rate" & parTab$province == "1","values"] <- 0.25
parTab[parTab$names == "growth_rate", "upper_bound"] <- 0.5
parTab[parTab$names == "growth_rate" & parTab$province != "1","values"] <- 0
parTab[parTab$names == "growth_rate" & parTab$province != "1","fixed"] <- 1
parTab[parTab$names == "growth_rate" & parTab$province != "1","lower_bound"] <- 0
parTab[parTab$names == "t0" & parTab$province != "1","values"] <- 0.00001
parTab[parTab$names == "t0" & parTab$province != "1","fixed"] <- 1
parTab[parTab$names == "t0" & parTab$province == "1",c("lower_bound","upper_bound")] <- c(0, 60)
parTab[parTab$names == "t0" & parTab$province == "1",c("lower_start","upper_start")] <- c(31,40)
parTab[parTab$names == "i0" & parTab$province != "1","values"] <- 0

parTab[parTab$names %in% c("shape","scale"),"fixed"] <- 0
parTab[parTab$names %in% c("shape","scale"),"values"] <- c(3,3)
parTab[parTab$names %in% c("shape","scale"),"lower_start"] <- c(2.5,2.5)
parTab[parTab$names %in% c("shape","scale"),"upper_start"] <- c(3.5,3.5)
parTab[parTab$names %in% c("lnorm_mean","lnorm_sd"),"fixed"] <- 1

parTab[parTab$names =="export_prob","fixed"] <- 1
parTab[parTab$names =="export_prob","values"] <- 1

confirm_delay_pars <- NULL

## Check that posterior works
confirmed_data1 <- confirmed_data %>% mutate(n=ifelse(province=="1", NA, n))
f <- create_model_func_provinces(parTab,confirmed_data1, confirm_delay_pars = confirm_delay_pars, 
                                 daily_import_probs = import_probs, daily_export_probs = export_probs,
                                 PRIOR_FUNC=prior_func)
f(parTab$values)

## Check that model works
f <- create_model_func_provinces(parTab,confirmed_data1, confirm_delay_pars = confirm_delay_pars, 
                                 daily_import_probs = import_probs, daily_export_probs = export_probs,
                                 PRIOR_FUNC=prior_func, ver="model",model_ver=2)
pars <- parTab$values
pars[which(parTab$names == "K")] <- 100000
dat <- f(pars)

ggplot(dat) + geom_line(aes(x=date,y=n,col=var)) + facet_wrap(~province,scales="free_y")

parTab[parTab$names == "t0" & parTab$province == "1","values"] <- 37
parTab[parTab$names == "t0" & parTab$province == "1","upper_bound"] <- 62
parTab[parTab$names == "growth_rate" & parTab$province == "1","fixed"] <- 1
startTab <- generate_start_tab(parTab)
startTab[startTab$names %in% c("weibull_alpha","weibull_sigma"),"values"] <- c(2.5,6)
#startTab[startTab$names %in% c("weibull_alpha","weibull_sigma"),"fixed"] <- 1

## MCMC
## Run first chain
mcmcPars <- c("iterations"=20000,"popt"=0.44,"opt_freq"=2000,
              "thin"=10,"adaptive_period"=10000,"save_block"=1000)
#startTab[startTab$province=="all","fixed"] <- 1
startTab[startTab$names == "weibull_alpha","lower_bound"] <- 1.5

output <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars, filename="logistic2",
                   CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=NULL,
                   PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                   confirm_delay_pars=confirm_delay_pars,
                   daily_import_probs = import_probs, daily_export_probs = export_probs,
                   ver="posterior",model_ver=2)


## Use this as input to multivariate chain
chain <- read.csv(output$file)
best_pars <- get_best_pars(chain)
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,0.5,w=0.8)

## Start from best location of previous chain
startTab$values <- best_pars

## Run second chain
mcmcPars <- c("iterations"=200000,"popt"=0.234,"opt_freq"=1000,
              "thin"=10,"adaptive_period"=100000,"save_block"=1000)
output <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars, filename="chains/serial_interval",
                   CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=mvrPars,
                   PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                   confirm_delay_pars=confirm_delay_pars,
                   daily_import_probs = import_probs, daily_export_probs = export_probs,
                   ver="posterior")


## Check convergence
chain <- read.csv(output$file)
pdf("tmp.pdf")
plot(coda::as.mcmc(chain[,c("shape","scale","weibull_alpha","weibull_sigma","t0","growth_rate","K","lnlike")]))
dev.off()

chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]
quants <- generate_prediction_intervals(chain, parTab, confirmed_data, confirm_delay_pars=confirm_delay_pars,
                                        daily_import_probs = import_probs, daily_export_probs = export_probs,
                                        nsamp=100,return_draws = FALSE,model_ver=2)
#samps <- quants$samp_ids
#samps <- sort(samps)
chain_save <- chain[chain$sampno %in% samps,]
chain_save$sampno <- match(chain_save$sampno, samps)

#quants <- quants$draws

quants$province <- row.names(import_probs)[quants$province]
confirmed_data2 <- confirmed_data

confirmed_data2$province <- row.names(import_probs)[confirmed_data2$province]
confirmed_data2$date <- as.Date(confirmed_data2$date, origin="2019-11-01")
quants$date <- as.Date(quants$date, origin="2019-11-01")

#quants %>% filter(var %in% c("total_prevalence","infection_prev","onset_prev") & sampno == 1) %>% ggplot() +
#  geom_line(aes(x=date,y=n,col=var)) + facet_wrap(~province,scales="free_y")
#write_csv(quants, "prevalence_estimates_summary_logistic_growth.csv")
#write_csv(quants, "prevalence_estimates_logistic_growth.csv")
#write_csv(chain_save, "mcmc_chain_thinned_logistic_growth.csv")

hubei_plot <- ggplot(quants[quants$province == "Hubei" & #quants$date <= "2020-01-25" & 
                              quants$date >= "2019-12-01" &
                              quants$var %in% c("confirmations","onsets","total_prevalence"),]) + 
  geom_vline(xintercept=(as.Date("2020-01-23", format="%Y-%m-%d", tz="LMT", origin="2019-11-01")),linetype="dashed") +
  geom_line(aes(x=date,y=median,col=var)) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
  scale_x_date(breaks="2 days") +
  theme_bw() +
  ylab("Daily prevalence") +
  xlab("Date") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position=c(0.2,0.2),
        panel.grid.minor=element_blank())

hubei_plot_inset <- ggplot(quants[quants$province == "Hubei" & quants$date <= "2020-01-01" & 
                                    quants$date >= "2019-12-01" &
                                    quants$var %in% c("confirmations","onsets"),]) + 
  geom_line(aes(x=date,y=median,col=var)) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
  scale_x_date(breaks="2 days") +
  scale_y_continuous(breaks=seq(0,160,by=20)) +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="none", 
        panel.grid.minor=element_blank())
vp <- viewport(width=0.5,height=0.5, x=0.4,y=0.7)
png("hubei_incidence.png", height=5,width=8,res=300,units="in")
print(hubei_plot)
print(hubei_plot_inset, vp=vp)
dev.off()

factor_levels <- confirmed_data2 %>% filter(province != "Hubei") %>% 
  group_by(province) %>%
  summarise(x=sum(n,na.rm=TRUE)) %>%
  arrange(-x) %>% pull(province)

quants <- quants[!(quants$province == "Hubei" & quants$var == "infections"),]
quants1 <- quants[!(quants$province == "Hubei" & quants$date > "2020-01-23"),]
quants1 <- quants1 %>% filter(province != "Hubei")
confirmed_data2 <- confirmed_data2 %>% filter(province != "Hubei")

confirmed_data2$province <- factor(confirmed_data2$province, levels=factor_levels)
quants1$province <- factor(quants1$province, levels=factor_levels)

p <- ggplot(quants1[quants1$var %in% c("infections","observations","onsets") & quants1$date >= "2020-01-01",]) +
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
  geom_line(aes(x=date,y=median,col=var)) +
  geom_point(data=confirmed_data2[confirmed_data2$date >= "2020-01-01",],aes(x=date,y=n),size=0.5) +
  scale_x_date(breaks="7 days") +
  geom_vline(xintercept=as.Date("2020-01-23",origin="2019-11-01"),linetype="dashed") +
  theme_pubr() + 
  theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~province,ncol=5,scales="free_y")
p


png("all_fits.png",height=10,width=12,units="in",res=300)
p + theme(legend.position="top")
dev.off()

library(data.table)
tmp <- chain[,colnames(chain) %like% "local_r"]
tmp <- tmp[,2:ncol(tmp)]
colnames(tmp) <- unique(confirmed_data2$province)
tmp1 <- reshape2::melt(tmp)
import_probs_melt <- reshape2::melt(import_probs)
import_probs_melt <- import_probs_melt %>% group_by(Var1) %>% summarise(val=mean(value)) %>% filter(Var1 != "Hubei") %>% ungroup()
colnames(import_probs_melt) <- c("variable", "Mean import probability")
tmp1 <- tmp1 %>% left_join(import_probs_melt)
tmp_order <- tmp1 %>% group_by(variable) %>% summarise(x=mean(value)) %>% arrange(-x) %>% pull(variable)
tmp1$variable <- factor(tmp1$variable, levels=factor_levels)
local_r_plot <- ggplot(tmp1) + 
  geom_violin(aes(x=variable, y=value, fill=`Mean import probability`),
              draw_quantiles=c(0.025,0.5,0.975),scale="width") + 
  scale_fill_gradient(low="blue",high="red") +
  geom_hline(yintercept=1,linetype="dashed") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="bottom") + 
  ylab("Average number of local cases per import") + xlab("Province") + 
  scale_y_continuous(expand=c(0,0),breaks=seq(0,10,by=1),limits=c(0,10))

png("local_r_plot.png",height=5,width=7,units="in",res=300)
local_r_plot
dev.off()

r0_ests <- read_csv("data/r0_ests.csv")

r0_ests <- r0_ests %>% select(Location,Value)
colnames(r0_ests) <- c("province","r0")

tmp2 <- plyr::ddply(tmp1, ~variable, function(x) mean(x$value))
colnames(tmp2) <- c("province","local_r")
wow <- merge(r0_ests, tmp2)

inc_p <- plot_incubation_period(chain,nsamp=200,xmax=40, prior_pars=inc_period_draws) + ylab("Probability density") + xlab("Delay from infection to symptom onset")
conf_p <- plot_confirm_delay(chain,nsamp=200,xmax=40) + ylab("Probability density") + xlab("Delay from symptom onset to confirmation")
delay_plots <- inc_p / conf_p
png("delay_plots.png",height=5,width=7,res=300,units="in")
delay_plots
dev.off()
