library(grid)
library(lazymcmc)
library(tidyverse)
library(ggpubr)
library(patchwork)

# setwd("~/Documents/GitHub/covback")
setwd("~/git_repos/covback")
Rcpp::compileAttributes()
devtools::document()
devtools::load_all()
#set.seed(2)

# Ada specify day 0
tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
# Ada specify day at which inference ends
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
# Ada specify times for inference
times <- seq(tmin, tmax, by="1 day")

# Ada read in incidence data from each province
# to do: actually generate this data from a forward simulation
confirmed_data <- read_csv("data/sim_data.csv")
# Ada plot data for posterity -- I think it needs more noise
ggplot(confirmed_data) + geom_line(aes(x=date,y=n)) + facet_wrap(~province,scales="free_y")

## Incubation period draws and parameter values from Backer et al.  
# Ada we actually estimate the incunation period distribution, the Backer et al.
# draws are only used to set a trong prior on said distribution 
inc_period_draws <- read.csv("data/backer_draws.csv",stringsAsFactors=FALSE)

# Ada start setting up parTab
parTab <- read.csv("pars/startTab_back_calc.csv",stringsAsFactors=FALSE)
startTab <- generate_start_tab(parTab)

## Make strong prior on alpha and sigma
prior_func <- create_prior_startdate(parTab, inc_period_draws, 37, 5)

## MCMC
# there is only one parameter to be fitted so it's a short chain
mcmcPars <- c("iterations"=2000,"popt"=0.234,"opt_freq"=10,
              "thin"=10,"adaptive_period"=1000,"save_block"=10)

# create_model_func_provinces is a closure that creates a function to evaluate the likelihood
# model_ver = 2 is logistic growth
set.seed(1)
stop()
output <- run_MCMC(parTab=startTab, data=confirmed_data, mcmcPars=mcmcPars, filename="temp2",
                   CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=NULL,
                   PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                   ver="posterior",model_ver=2)