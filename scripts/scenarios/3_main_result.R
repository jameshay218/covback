library(grid)
library(lazymcmc)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(doParallel)

## Either load package locally or install from github
setwd("~/Documents/GitHub/covback/")
devtools::load_all()
#install.packages("~/Documents/GitHub/covback/",repos=NULL,type="source")
#library(covback)

mcmcPars1 <- c("iterations"=100000,"popt"=0.44,"opt_freq"=1000,
               "thin"=100,"adaptive_period"=200000,"save_block"=100)
mcmcPars2 <- c("iterations"=200000,"popt"=0.234,"opt_freq"=1000,
               "thin"=100,"adaptive_period"=400000,"save_block"=100)

## Table giving all scenarios with parameter settings, enumerated out for each chain number
scenario_key <- read_csv("~/Documents/GitHub/covback/scripts/scenarios/scenario_key.csv")
scenario_key <- scenario_key[scenario_key$runname %in% c("diffuse_r_local","fewer_travellers","early_seed") & scenario_key$chain_no %in% c(3),]

## Set up parallelisation
n_clusters <- 4
cl <- makeCluster(n_clusters)
registerDoParallel(cl)

## Timeframe of model fitting
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

## Table giving 
time_varying_report_pars <- data.frame(date=as.Date(times,origin="2019-11-01"),shape=3.18,scale=1/0.59)
time_varying_report_pars[time_varying_report_pars$date <= as.Date("2020-01-27",origin="2019-11-01"),"shape"] <- 3.72
time_varying_report_pars[time_varying_report_pars$date <= as.Date("2020-01-27",origin="2019-11-01"),"scale"] <- 1/0.42

confirmed_data1 <- as.data.frame(read_csv("data/real/midas_data_final.csv"))
confirmed_data1 <- confirmed_data1 %>% select(-province_raw)
confirmed_data1 <- confirmed_data1 %>% mutate(n=ifelse(province==1, NA, n))

## Set up priors
subset_parTab <- parTab[parTab$names == "local_r",]
r_index <- which(subset_parTab$province != "1")
par_names <- parTab$names

parTab[parTab$names == "t_switch","fixed"] <- 0
parTab[parTab$names == "growth_rate" & parTab$province == 1,"fixed"] <- 0

####################
## NOT USED
####################
## Create priors on delay distribution parameters
# prior_table <- read.csv("pars/prior_quantiles.csv",stringsAsFactors = FALSE)
# prior_sds <- find_all_prior_sds(prior_table)
# 
# prior_func_delays <- function(pars){
#   serial_pars <- gamma_pars_from_mean_sd(pars["serial_interval_mean"],pars["serial_interval_sd"]^2)
#   serial_interval_1 <- dnorm(serial_pars[[1]], prior_table$mean[1], prior_sds[1], 1)
#   serial_interval_2 <- dnorm(1/serial_pars[[2]], prior_table$mean[2], prior_sds[2], 1)
#   
#   incubation_period_1 <- dnorm(pars["lnorm_incu_par1"], prior_table$mean[3], prior_sds[3], 1)
#   incubation_period_2 <- dnorm(pars["lnorm_incu_par2"], prior_table$mean[4], prior_sds[4], 1)
#   
#   return(serial_interval_1 + serial_interval_2 + incubation_period_1 + incubation_period_2)
# }
# pars <- parTab$values
# names(pars) <- parTab$names
# prior_func_delays(pars)

res <- foreach(i=1:nrow(scenario_key),.packages=c("covback","lazymcmc","tidyverse")) %dopar% {
  setwd("../covback_chains_final/main_results_final_maybe")
  filename_tmp <- paste0(scenario_key$runname[i], "_",scenario_key$chain_no[i])
  
  ## Update model parameters based on scenario
  parTab[parTab$names == "t0","values"] <- scenario_key$t0_val[i]
  parTab[parTab$names == "t0","fixed"] <- scenario_key$t0_fixed[i]
  parTab[parTab$names == "local_r_sd","values"] <- scenario_key$r_local_sd[i]
  parTab[parTab$names == "t_switch","fixed"] <- scenario_key$t_switch_fixed[i]
  parTab[parTab$names == "t_switch","values"] <- scenario_key$t_switch_val[i]
  
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
  plot(coda::as.mcmc(chain[,c("serial_interval_mean","serial_interval_sd","t0","K","local_r","lnlike")]))
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
  plot(coda::as.mcmc(chain[,c("serial_interval_mean","serial_interval_sd","t0","K","local_r","lnlike")]))
  dev.off()
}