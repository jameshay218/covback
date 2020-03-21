library(grid)
library(lazymcmc)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(doParallel)
#library(covback)
setwd("~/Documents/GitHub/covback/")
install.packages("~/Documents/GitHub/covback/",repos=NULL,type="source")
library(covback)
#devtools::document()
devtools::load_all()

mcmcPars1 <- c("iterations"=100000,"popt"=0.44,"opt_freq"=1000,
              "thin"=10,"adaptive_period"=50000,"save_block"=1000)
mcmcPars2 <- c("iterations"=200000,"popt"=0.234,"opt_freq"=1000,
              "thin"=10,"adaptive_period"=100000,"save_block"=1000)


nchains <- 4
filenames <- c("seed_times_1","seed_times_2",
               "travellers","serial_interval",
               "time_to_peak",
               "rlocal_prior_diffuse")
chain_nos <- rep(1:nchains, length(filenames))
filenames <- rep(filenames, each=nchains)
## Set up parallelisation
n_clusters <- 8
#registerDoParallel(cores=n_clusters)
cl <- makeCluster(n_clusters)
registerDoParallel(cl)

tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")

## Incubation period draws and parameter values
inc_period_draws <- NULL

## Real export probs
export_probs <- read_csv("data/export_probs_matched.csv")$export_prob
export_probs_lower <- read_csv("data/export_probs_lower.csv")$export_prob
## Real import probs
import_probs <- read_csv("data/import_probs_matched.csv")
import_probs <- as.matrix(import_probs[,2:ncol(import_probs)])

#parTab <- create_many_province_partab(parTab, length(use_provinces), tmax)
parTab <- read_csv("pars/partab_logistic_growth.csv")

provinces <- unique(parTab$province)
provinces <- provinces[provinces != "all"]
n_provinces <- length(provinces)

## Make strong prior on alpha and sigma
prior_func <- create_random_effects_rlocal(parTab,province_prior="2",r_local_mean=0.4,
                                           r_local_sd=0.05,func_ver="gamma")

## Check that posterior works
confirmed_data1 <- as.data.frame(read_csv("data/real/confirmed_data_matched_final.csv"))
confirmed_data1 <- confirmed_data1 %>% select(-province_name)
confirmed_data1 <- confirmed_data1 %>% mutate(n=ifelse(province==1, NA, n))
parTab[parTab$names == "t_switch","fixed"] <- 1
parTab[parTab$names == "t_switch","lower_bound"] <- 46
parTab[parTab$names == "t_switch","values"] <- 50
parTab[parTab$names == "t_switch","upper_bound"] <- 53
## Where to save all chains
setwd("../covback_chains_final")

res <- foreach(i=seq_along(filenames),.packages=c("covback","lazymcmc","tidyverse")) %dopar% {
  #for(i in seq_along(filenames)){
  setwd("~/Documents/GitHub/covback_chains_final/sensitivity")
  filename <- filenames[i]
  chain_no <- chain_nos[i]
  filename_tmp <- paste0(filename, "_",chain_no)
  ## Set up properties for this scenario
  parTab1 <- parTab
  import_probs1 <- import_probs
  export_probs1 <- export_probs
  prior_func1 <- prior_func
  
  if(filename == "seed_times_1"){
    print("seed_time_1")
    parTab1[parTab1$names == "t0" & parTab1$province == "1","values"] <- 16
    parTab1[parTab1$names == "t_switch","values"] <- 69
    parTab1[parTab1$names == "t_switch","lower_bound"] <- 60
    parTab1[parTab1$names == "t_switch","upper_bound"] <- 80
  } else if(filename == "seed_times_2") {
    print("seed_time_2")
    parTab1[parTab1$names == "t0" & parTab1$province == "1","values"] <- 30
    parTab1[parTab1$names == "t_switch","values"] <- 55
    parTab1[parTab1$names == "t_switch","lower_bound"] <- 50
    parTab1[parTab1$names == "t_switch","upper_bound"] <- 80
  } else if(filename == "travellers") {
    print("travellers")
    export_probs1 <- export_probs_lower
  } else if(filename=="serial_interval"){
    print("serial_interval")
    parTab1[parTab1$names %in% c("confirm_delay_shape","confirm_delay_scale"), "fixed"] <- 1
    parTab1[parTab1$names %in% c("serial_interval_gamma_alpha","serial_interval_gamma_scale"), "fixed"] <- 0
  } else if(filename=="time_to_peak"){
    print("time_to_peak")
    parTab1[parTab1$names == "t_switch","fixed"] <- 0
  } else if(filename=="rlocal_prior"){
    print("rlocal_prior")
    prior_func1 <- create_random_effects_rlocal(parTab,province_prior="2",r_local_mean=0.4,
                                               r_local_sd=0.25,func_ver="gamma")
  } else {
    parTab1 <- parTab
  }
  #prior_func <- NULL
  startTab <- generate_start_tab(as.data.frame(parTab1))
  
  ## MCMC
  ## Run first chain
  output <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars1, filename=filename_tmp,
                     CREATE_POSTERIOR_FUNC=create_model_func_provinces_fixed, mvrPars=NULL,
                     PRIOR_FUNC = prior_func1, OPT_TUNING=0.2,
                     daily_import_probs = import_probs1, daily_export_probs = export_probs1,
                     incubation_ver="lnorm",
                     noise_ver="poisson",model_ver=2)
  
  
  ## Use this as input to multivariate chain
  chain <- read.csv(paste0(filename_tmp,"_univariate_chain.csv"))
  pdf(paste0(filename_tmp,"_chain.pdf"))
  plot(coda::as.mcmc(chain[,c("confirm_delay_shape","confirm_delay_scale","local_r.1","K")]))
  dev.off()
  #chain <- read.csv(output$file)
  best_pars <- get_best_pars(chain)
  chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
  covMat <- cov(chain)
  mvrPars <- list(covMat,0.5,w=0.8)
  
  ## Start from best location of previous chain
  startTab$values <- best_pars
  
  ## Run second chain
  output <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars2, filename=filename_tmp,
                     CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=mvrPars,
                     PRIOR_FUNC = prior_func1, OPT_TUNING=0.2,
                     daily_import_probs = import_probs1, daily_export_probs = export_probs1,
                     incubation_ver="lnorm",
                     noise_ver="poisson",model_ver=2)
  
  
  ## Check convergence
  chain <- read.csv(output$file)
  pdf(paste0(filename_tmp,"_chain.pdf"))
  plot(coda::as.mcmc(chain[,c("confirm_delay_shape","confirm_delay_scale","local_r.1","K")]))
  dev.off()
}

