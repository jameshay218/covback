setwd("~/Documents/covback")
Rcpp::compileAttributes()
devtools::document()
devtools::load_all()

library(lazymcmc)
library(tidyverse)
tmin <- as.POSIXct("19-11-01",format="%Y-%m-%d")
tmax <- as.POSIXct("20-02-20",format="%Y-%m-%d")
times <- seq(tmin, tmax, by="1 day")
confirmed_data <- read_csv("data/confirmed_data.csv")
confirmed_data <- confirmed_data %>% filter(country_region == "Mainland China")
confirmed_data$date <-as.POSIXct(as.character(confirmed_data$date))

confirmed_data$date <- match(confirmed_data$date, times)
confirmed_data <- confirmed_data %>% select(province, date, diff)
colnames(confirmed_data)[3] <- "n"

all_reports <- expand.grid(province=unique(confirmed_data$province),
                           date=match(times, times))
confirmed_data <- confirmed_data %>% right_join(all_reports)
confirmed_data$date <- confirmed_data$date - 1
#confirmed_data <- confirmed_data %>% filter(date <= 80)

ggplot(confirmed_data) + geom_line(aes(x=date,y=n)) + facet_wrap(~province,scales="free_y")

## Incubation period draws and parameter values
inc_period_draws <- read.csv("~/Documents/case_to_infection/data/backer_weibull_draws.csv",stringsAsFactors=FALSE)
parTab <- read.csv("pars/partab_provinces.csv",stringsAsFactors=FALSE)

## Real export probs
export_probs <- read.csv("data/export_probs.csv")
export_probs <- export_probs[1:112,]
tmax <- nrow(export_probs)-1
export_probs <- export_probs$prob_leaving

## Real import probs
import_probs <- read.csv("data/import_probs.csv",header = TRUE)
import_probs <- import_probs[1:112,]
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
prior_func <- create_incubation_prior(inc_period_draws)
## Generate some fake data

confirm_delay_pars <- read.csv("data/fitted_confirm_delays.csv")
confirm_delay_pars <- confirm_delay_pars %>%  select(shape, scale)
confirm_delay_pars$date_onset <- 0:(nrow(confirm_delay_pars)-1)
confirm_delay_pars <- confirm_delay_pars %>% filter(date_onset <= tmax)
plot_reporting_landscape(confirm_delay_pars$shape, confirm_delay_pars$scale)

parTab[parTab$names == "t0" & parTab$province == "1","fixed"] <- 0
parTab[parTab$names == "growth_rate", "upper_bound"] <- 0.3
parTab[parTab$names == "growth_rate" & parTab$province != "1","values"] <- 0.00001
parTab[parTab$names == "growth_rate" & parTab$province != "1","fixed"] <- 1
parTab[parTab$names == "t0" & parTab$province != "1","values"] <- 0.00001
parTab[parTab$names == "t0" & parTab$province != "1","fixed"] <- 1
parTab[parTab$names == "t0" & parTab$province == "1",c("lower_start","upper_start")] <- c(30,40)

parTab[parTab$names %in% c("shape","scale"),"fixed"] <- 0
confirm_delay_pars <- NULL

## Check that posterior works
f <- create_model_func_provinces(parTab,confirmed_data, confirm_delay_pars = confirm_delay_pars, 
                                 daily_import_probs = import_probs, daily_export_probs = export_probs,
                                 PRIOR_FUNC=prior_func)
f(parTab$values)

## Check that model works
f <- create_model_func_provinces(parTab,confirmed_data, confirm_delay_pars = confirm_delay_pars, 
                                 daily_import_probs = import_probs, daily_export_probs = export_probs,
                                 PRIOR_FUNC=prior_func, ver="model")
dat <- f(parTab$values)

ggplot(dat) + geom_line(aes(x=date,y=n,col=var)) + facet_wrap(~province,scales="free_y")

parTab[parTab$names == "t0" & parTab$province == "1","values"] <- 0
parTab[parTab$names == "t0" & parTab$province == "1","upper_bound"] <- 62
startTab <- generate_start_tab(parTab)
startTab[startTab$names %in% c("weibull_alpha","weibull_sigma"),"values"] <- c(2.5,6)
#startTab[startTab$names %in% c("weibull_alpha","weibull_sigma"),"fixed"] <- 1

## MCMC
## Run first chain
mcmcPars <- c("iterations"=30000,"popt"=0.44,"opt_freq"=1000,
              "thin"=10,"adaptive_period"=10000,"save_block"=100)
#startTab[startTab$province=="all","fixed"] <- 1

confirmed_data1 <- confirmed_data %>% mutate(n=ifelse(province=="1", NA, n))

output <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars, filename="test",
                   CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=NULL,
                   PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                   confirm_delay_pars=confirm_delay_pars,
                   daily_import_probs = import_probs, daily_export_probs = export_probs,
                   ver="posterior")


## Use this as input to multivariate chain
chain <- read.csv(output$file)
best_pars <- get_best_pars(chain)
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,0.01,w=0.8)

## Start from best location of previous chain
startTab$values <- best_pars

## Run second chain
mcmcPars <- c("iterations"=50000,"popt"=0.234,"opt_freq"=5000,
              "thin"=10,"adaptive_period"=20000,"save_block"=100)
output <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars, filename="test",
                   CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=mvrPars,
                   PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                   confirm_delay_pars=confirm_delay_pars,
                   daily_import_probs = import_probs, daily_export_probs = export_probs,
                   ver="posterior")


## Check convergence
chain <- read.csv(output$file)
pdf("tmp.pdf")
plot(coda::as.mcmc(chain[,c("weibull_alpha","weibull_sigma","t0","growth_rate","lnlike")]))
dev.off()

chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]
quants <- generate_prediction_intervals(chain, parTab, confirmed_data, confirm_delay_pars=confirm_delay_pars,
                                        daily_import_probs = import_probs, daily_export_probs = export_probs,
                                        nsamp=100)
quants <- quants[!(quants$province == "1" & quants$var == "infections"),]
quants$province <- row.names(import_probs)[quants$province]
confirmed_data2 <- confirmed_data

confirmed_data2$province <- row.names(import_probs)[confirmed_data2$province]
confirmed_data2$date <- as.Date(confirmed_data2$date, origin="2019-11-01")
quants$date <- as.Date(quants$date, origin="2019-11-01")
quants1 <- quants[!(quants$province == "Hubei" & quants$date > "2020-01-23"),]
p <- ggplot(quants1[quants1$var %in% c("infections","observations","onsets"),]) +
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
  geom_line(aes(x=date,y=median,col=var)) +
  geom_line(data=confirmed_data2,aes(x=date,y=n),size=0.5) +
  scale_x_date(breaks="7 days") +
  geom_vline(xintercept=as.Date("2020-01-23",origin="2019-11-01"),linetype="dashed") +
  theme_bw() + 
  theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~province,ncol=5,scales="free_y")
p

p <- plot_model_fit(chain, parTab, confirmed_data,confirm_delay_pars=confirm_delay_pars,
                    daily_import_probs = import_probs, daily_export_probs = export_probs,
                    nsamp=100)
pdf("tmp.pdf",height=12,width=10)
p + theme(legend.position=c(0.8,0.1))
dev.off()
