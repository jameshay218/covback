library(grid)
library(lazymcmc)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(doParallel)

chain_save_wd <- "~/Documents/GitHub/covback_chains_final/final_20200522/"

## Either load package locally or install from github
setwd("~/Documents/GitHub/covback/")
#install.packages("~/Documents/GitHub/covback/",repos=NULL,type="source")
library(covback)

#mcmcPars1 <- c("iterations"=100000,"popt"=0.44,"opt_freq"=1000,
#               "thin"=100,"adaptive_period"=50000,"save_block"=100)
#mcmcPars2 <- c("iterations"=350000,"popt"=0.234,"opt_freq"=1000,
#               "thin"=100,"adaptive_period"=150000,"save_block"=100)
mcmcPars1 <- c("iterations"=80000,"popt"=0.44,"opt_freq"=1000,
               "thin"=100,"adaptive_period"=50000,"save_block"=1000)
mcmcPars2 <- c("iterations"=100000,"popt"=0.234,"opt_freq"=1000,
               "thin"=100,"adaptive_period"=100000,"save_block"=1000)

## Table giving all scenarios with parameter settings, enumerated out for each chain number
scenario_key <- read_csv("~/Documents/GitHub/covback/scripts/scenario_key.csv")
#scenario_key <- scenario_key %>% filter(chain_no == 1)

## Set up parallelisation
n_clusters <- 16
cl <- makeCluster(n_clusters)
registerDoParallel(cl)

## Timeframe of model fitting
## Day 0 is 2019-11-01
tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")

## Real export probs (within China)
export_probs <- read_csv("data/export_probs_matched.csv")$export_prob
export_probs_lower <- read_csv("data/export_probs_lower.csv")$export_prob

## Real import probs (within China)
import_probs <- read_csv("data/import_probs_matched.csv")
import_probs <- as.matrix(import_probs[,2:ncol(import_probs)])
colnames(import_probs) <- NULL

## Parameter table to control pars during MCMC
parTab <- read_csv("pars/partab_logistic_growth.csv")
parTab[parTab$names == "K","values"] <- log(parTab[parTab$names == "K","values"])
parTab[parTab$names == "K","upper_bound"] <- log(parTab[parTab$names == "K","upper_bound"])
parTab[parTab$names == "K","lower_bound"] <- log(parTab[parTab$names == "K","lower_bound"])
parTab[parTab$names == "K","lower_start"] <- log(parTab[parTab$names == "K","lower_start"])
parTab[parTab$names == "K","upper_start"] <- log(parTab[parTab$names == "K","upper_start"])

## Table giving confirmation delay distribution as it changes over time, as per Zhang et al. 2020
time_varying_report_pars <- data.frame(date=as.Date(times,origin="2019-11-01"),shape=3.18,scale=1/0.59)
time_varying_report_pars[time_varying_report_pars$date <= as.Date("2020-01-27",origin="2019-11-01"),"shape"] <- 3.72
time_varying_report_pars[time_varying_report_pars$date <= as.Date("2020-01-27",origin="2019-11-01"),"scale"] <- 1/0.42

## Confirmed case data
confirmed_data1 <- as.data.frame(read_csv("data/real/midas_data_final.csv"))
confirmed_data1 <- confirmed_data1 %>% select(-province_raw)
confirmed_data1 <- confirmed_data1 %>% mutate(n=ifelse(province==1, NA, n))

## Set up priors
subset_parTab <- parTab[parTab$names == "local_r",]
r_index <- which(subset_parTab$province != "1")
par_names <- parTab$names

####################
## RUN MCMC
####################
res <- foreach(i=1:nrow(scenario_key),.packages=c("covback","lazymcmc","tidyverse")) %dopar% {
  setwd(chain_save_wd)
  dir.create(scenario_key$runname[i])
  setwd(scenario_key$runname[i])
  filename_tmp <- paste0(scenario_key$runname[i], "_",scenario_key$chain_no[i])
  
  ## Update model parameters based on scenario
  parTab[parTab$names == "t0","values"] <- scenario_key$t0_val[i]
  parTab[parTab$names == "t0","fixed"] <- scenario_key$t0_fixed[i]
  parTab[parTab$names == "local_r_sd","values"] <- scenario_key$r_local_sd[i]
  parTab[parTab$names == "t_switch","fixed"] <- scenario_key$t_switch_fixed[i]
  parTab[parTab$names == "t_switch","values"] <- scenario_key$t_switch_val[i]
  parTab[parTab$names %in% c("serial_interval_mean","serial_interval_var"),"fixed"] <- scenario_key$serial_interval_fixed[i]
  
  ## Hyperprior on R_local parameters
  prior_func_rlocal <- function(pars){
    r_local_sd1 <- pars["local_r_sd"]
    r_local_mean1 <- pars["local_r_mean"]
    gamma_pars <- gamma_pars_from_mean_sd(r_local_mean1, r_local_sd1^2)
    r_locals <- pars[which(names(pars) == "local_r")]
    lik <- dgamma(r_locals[r_index],gamma_pars[[1]],scale=gamma_pars[[2]],log=TRUE)
    return(sum(lik))
  }
  
  ## Use gamma hyperprior for R_local if specified
  if(scenario_key$r_local_prior[i] == 1){
    prior_func <- prior_func_rlocal
  } else {
    prior_func <- NULL
  }
  
  ## 5mil or 4mil travellers scenario
  if(scenario_key$travellers[i] == 1){
    export_probs_use <- export_probs
  } else {
    export_probs_use <- export_probs_lower
  }
  
  ## Random starting locations
  startTab <- generate_start_tab(as.data.frame(parTab))
  f <- create_model_func_provinces_fixed(parTab, data=confirmed_data1, PRIOR_FUNC = prior_func, daily_import_probs = import_probs, 
                                       daily_export_probs = export_probs_use,
                                         time_varying_confirm_delay_pars=time_varying_report_pars,
                                         incubation_ver="lnorm",ver="model",
                                         noise_ver="poisson",model_ver="logistic",calculate_prevalence = TRUE,scale_reporting = FALSE)
  dat <- f(parTab$values)

  ## MCMC
  ## Run first chain
  output <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars1, filename=filename_tmp,
                     CREATE_POSTERIOR_FUNC=create_model_func_provinces_fixed, mvrPars=NULL,
                     PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                     daily_import_probs = import_probs, daily_export_probs = export_probs_use,
                     time_varying_confirm_delay_pars=time_varying_report_pars,
                     incubation_ver="lnorm",
                     noise_ver="poisson",model_ver="logistic",solve_prior=FALSE)
  
  
  ## Use this as input to multivariate chain
  chain <- read.csv(paste0(filename_tmp,"_univariate_chain.csv"))
  
  ## Check convergence
  pdf(paste0(filename_tmp,"_chain.pdf"))
  plot(coda::as.mcmc(chain[,c("serial_interval_mean","serial_interval_var","t0","K","local_r","lnlike")]))
  dev.off()
  
  best_pars <- get_best_pars(chain)
  chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
  
  covMat <- cov(chain)
  mvrPars <- list(covMat,0.5,w=0.8, 0.8)
  
  ## Start from best location of previous chain
  startTab$values <- best_pars
  
  ## Run second chain
  output <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars2, filename=filename_tmp,
                     CREATE_POSTERIOR_FUNC=create_model_func_provinces_fixed, mvrPars=mvrPars,
                     PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                     daily_import_probs = import_probs, daily_export_probs = export_probs_use,
                     time_varying_confirm_delay_pars=time_varying_report_pars,
                     incubation_ver="lnorm",
                     noise_ver="poisson",model_ver="logistic",solve_prior=FALSE)
    
  ## Check convergence
  chain <- read.csv(output$file)
  pdf(paste0(filename_tmp,"_chain.pdf"))
  plot(coda::as.mcmc(chain[,c("serial_interval_mean","serial_interval_var","t0","K","local_r","lnlike")]))
  dev.off()
}
