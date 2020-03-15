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

## Incubation period draws and parameter values from Backer et al.  
# Ada we actually estimate the incunation period distribution, the Backer et al.
# draws are only used to set a trong prior on said distribution 
inc_period_draws <- read.csv("data/backer_draws.csv",stringsAsFactors=FALSE)

# Ada start setting up parTab
startTab <- parTab <- read.csv("pars/startTab.csv",stringsAsFactors=FALSE)

## Real export probs
export_probs <- read.csv("data/export_probs_matched.csv") %>% unlist %>% unname
import_probs <- read.csv("data/import_probs_matched.csv",row.names = 1) %>% as.matrix

# Ada: the below is unifying the incidence, export prob and import prob data
# to make sure they cover the same provinces
confirmed_data <- read_csv("data/confirmed_data_matched.csv", col_types = "idd")

## Make strong prior on alpha and sigma
prior_func <- create_prior_startdate(parTab, inc_period_draws, 37, 5)

# Ada remove Hubei data (we are not using Hubei data for inference)
confirmed_data1 <- confirmed_data %>% mutate(n=ifelse(province=="1", NA, n))
stop()
## MCMC
## Run first chain
# Ada model fitting settings
mcmcPars <- c("iterations"=20000,"popt"=0.44,"opt_freq"=2000,
              "thin"=10,"adaptive_period"=10000,"save_block"=1000)


# Ada run MCMC
# initial short run to find a good place in parameter space, updating one parameter at a time
# create_model_func_provinces is a closure that creates a function to evaluate the likelihood
# model_ver = 2 is logistic growth
set.seed(1)
output <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars, filename="logistic2",
                   CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=NULL,
                   PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                   daily_import_probs = import_probs, daily_export_probs = export_probs,
                   ver="posterior",model_ver=2)

## Use maximum likellihood parss from short chain as starting point for MCMC
# where all paramteres are updated simultaneously
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

output <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars, filename="28_provinces",
                   CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=mvrPars,
                   PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                   daily_import_probs = import_probs, daily_export_probs = export_probs,
                   ver="posterior", model_ver=2)
