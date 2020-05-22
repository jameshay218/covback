library(grid)
library(lazymcmc)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(doParallel)
#library(covback)
setwd("~/Documents/GitHub/covback/")
#install.packages("~/Documents/GitHub/covback/",repos=NULL,type="source")
#library(covback)
#devtools::document()
devtools::load_all()

mcmcPars1 <- c("iterations"=100000,"popt"=0.44,"opt_freq"=1000,
               "thin"=1,"adaptive_period"=50000,"save_block"=100)
mcmcPars2 <- c("iterations"=350000,"popt"=0.234,"opt_freq"=1000,
               "thin"=1,"adaptive_period"=150000,"save_block"=100)

nchains <- 2
filenames <- c("unconstrained","tight_sd","unknown_sd",
               "relaxed_sd","unconstrained_fixed","tight_sd_fixed",
               "unknown_sd_fixed",
               "relaxed_sd_fixed")
chain_nos <- rep(1:nchains, length(filenames))
filenames <- rep(filenames, each=nchains)

## Set up parallelisation
n_clusters <- 16

cl <- makeCluster(n_clusters)
registerDoParallel(cl)

tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")

## Incubation period draws and parameter values
inc_period_draws <- NULL

## Real export probs
export_probs <- read_csv("data/export_probs_matched.csv")$export_prob

## Real import probs
import_probs <- read_csv("data/import_probs_matched.csv")
import_probs <- as.matrix(import_probs[,2:ncol(import_probs)])
colnames(import_probs) <- NULL

parTab <- read_csv("pars/partab_logistic_growth.csv")
#parTab[parTab$names == "growth_rate","fixed"] <- 1
#parTab[parTab$names == "t_switch","fixed"] <- 0

confirmed_data1 <- as.data.frame(read_csv("data/real/midas_data_final.csv"))
confirmed_data1 <- confirmed_data1 %>% select(-province_raw)
confirmed_data1 <- confirmed_data1 %>% mutate(n=ifelse(province==1, NA, n))


confirmed_data1 <- confirmed_data1[confirmed_data1$province %in% 1:6,]
parTab <- parTab[parTab$province %in% c("all", 1:6),]
import_probs <- import_probs[1:6,]

## Set up priors
subset_parTab <- parTab[parTab$names == "local_r",]
r_index <- which(subset_parTab$province != "1")
par_names <- parTab$names

## Base prior is to have Gaussian with mean 25th Jan and sd = 2 on
## time of peak onset incidence. The +5 days is the mode of the 
## fixed incubation period
prior_func_base <- function(pars){
  names(pars) <- parTab$names
  t_switch <- log(pars["K"]-1)/pars["growth_rate"]  + pars[which(names(pars)=="t0")[1]] + 5
  prior_prob <- dnorm(t_switch, 86, 2, 1)
}

prior_func_rlocal <- function(pars){
  r_local_sd1 <- pars["local_r_sd"]
  r_local_mean1 <- pars["local_r_mean"]
  gamma_pars <- gamma_pars_from_mean_sd(r_local_mean1, r_local_sd1^2)
  r_locals <- pars[which(names(pars) == "local_r")]
  lik <- dgamma(r_locals[r_index],gamma_pars[[1]],scale=gamma_pars[[2]],log=TRUE)
  lik2 <- prior_func_base(pars)
  return(sum(lik) + lik2)
}


f <- create_model_func_provinces_fixed(parTab, data=confirmed_data1, daily_import_probs = import_probs, daily_export_probs = export_probs,
                                       incubation_ver="lnorm",PRIOR_FUNC=NULL,
                                       noise_ver="poisson",model_ver=2)
f(parTab$values)


parTab1 <- parTab
res <- foreach(i=seq_along(filenames),.packages=c("covback","lazymcmc","tidyverse")) %dopar% {
  setwd("../covback_chains_final/testing_still_works")
  filename_tmp <- paste0(filenames[i], "_",chain_nos[i])
  parTab <- parTab1
  parTab[parTab$names == "size","fixed"] <- 1
  parTab[parTab$names %in% c("confirm_delay_shape","confirm_delay_scale"), "fixed"] <- 0

  if(filenames[i] %in% c("relaxed_sd","relaxed_sd_fixed")){
    print("relaxed_sd")
    parTab[parTab$names == "local_r_sd","values"] <- 5
    prior_func <- prior_func_rlocal
    if(filenames[i] == "relaxed_sd_fixed"){
      print("relaxed_sd_fixed")
      parTab[parTab$names == "t0","fixed"] <- 1
    }
  } else if (filenames[i] %in% c("tight_sd","tight_sd_fixed")) {
    print("tight_sd")
    parTab[parTab$names == "local_r_sd","values"] <- 0.5
    prior_func <- prior_func_rlocal
    if(filenames[i] == "tight_sd_fixed"){
      print("tight_sd_fixed")
      parTab[parTab$names == "t0","fixed"] <- 1
    }
  } else if (filenames[i] %in% c("unknown_sd","unknown_sd_fixed")){
    print("unknown_sd")
    parTab[parTab$names == "local_r_sd","fixed"] <- 0
    prior_func <- prior_func_rlocal
    if(filenames[i] == "unknown_sd_fixed"){
      print("unknown_sd_fixed")
      parTab[parTab$names == "t0","fixed"] <- 1
    }
  } else {
    print("unconstrained")
    prior_func <- prior_func_base
    if(filenames[i] == "unconstrained_fixed"){
      print("unconstrained_fixed")
      parTab[parTab$names == "t0","fixed"] <- 1
    }
  }
  startTab <- generate_start_tab(as.data.frame(parTab))
  
  ## MCMC
  ## Run first chain
  output <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars1, filename=filename_tmp,
                     CREATE_POSTERIOR_FUNC=create_model_func_provinces_fixed, mvrPars=NULL,
                     PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                     daily_import_probs = import_probs, daily_export_probs = export_probs,
                     incubation_ver="lnorm",
                     noise_ver="poisson",model_ver=2)
  
  
  ## Use this as input to multivariate chain
  chain <- read.csv(paste0(filename_tmp,"_univariate_chain.csv"))
  
  pdf(paste0(filename_tmp,"_chain.pdf"))
  plot(coda::as.mcmc(chain[,c("confirm_delay_shape","confirm_delay_scale","weibull_alpha","weibull_sigma","t0","growth_rate","K","lnlike")]))
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
                     CREATE_POSTERIOR_FUNC=create_model_func_provinces_fixed, mvrPars=mvrPars,
                     PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                     daily_import_probs = import_probs, daily_export_probs = export_probs,
                     incubation_ver="lnorm",
                     noise_ver="poisson",model_ver=2)
    
  ## Check convergence
  chain <- read.csv(output$file)
  pdf(paste0(filename_tmp,"_chain.pdf"))
  plot(coda::as.mcmc(chain[,c("confirm_delay_shape","confirm_delay_scale","weibull_alpha","weibull_sigma","t0","growth_rate","K","lnlike")]))
  dev.off()
}

parTab <- parTab1
parTab[parTab$names == "growth_rate","fixed"] <- 1
parTab[parTab$names == "t_switch","fixed"] <- 1
parTab[parTab$names == "t0","fixed"] <- 1
parTab1 <- parTab

## Base prior is to have Gaussian with mean 25th Jan and sd = 2 on
## time of peak onset incidence. The +5 days is the mode of the 
## fixed incubation period
prior_func_base <- function(pars){
  names(pars) <- parTab$names
  #t_switch <- log(pars["K"]-1)/pars["growth_rate"]  + pars[which(names(pars)=="t0")[1]] + 5
  #prior_prob <- dnorm(t_switch, 86, 2, 1)
  prior_prob <- 0
  prior_prob
}

prior_func_rlocal2 <- function(pars){
  r_local_sd1 <- pars["local_r_sd"]
  r_local_mean1 <- pars["local_r_mean"]
  gamma_pars <- gamma_pars_from_mean_sd(r_local_mean1, r_local_sd1^2)
  r_locals <- pars[which(names(pars) == "local_r")]
  lik <- dgamma(r_locals[r_index],gamma_pars[[1]],scale=gamma_pars[[2]],log=TRUE)
  return(sum(lik) )
}



setwd("~/Documents/GitHub/covback/")
res <- foreach(i=seq_along(filenames),.packages=c("covback","lazymcmc","tidyverse")) %dopar% {
#  i <- 2
  setwd("../covback_chains_final/testing_2")
  filename_tmp <- paste0(filenames[i], "_",chain_nos[i])
  parTab <- parTab1
  parTab[parTab$names == "size","fixed"] <- 1
  parTab[parTab$names %in% c("confirm_delay_shape","confirm_delay_scale"), "fixed"] <- 0
  parTab[parTab$names == "local_r" & parTab$province == "3","fixed"] <- 1
  parTab[parTab$names == "local_r" & parTab$province == "3","values"] <- 0.001
  
  if(filenames[i] %in% c("relaxed_sd","relaxed_sd_fixed")){
    print("relaxed_sd")
    parTab[parTab$names == "local_r_sd","values"] <- 5
    prior_func <- prior_func_rlocal
    if(filenames[i] == "relaxed_sd_fixed"){
      print("relaxed_sd_fixed")
      parTab[parTab$names == "t0","fixed"] <- 1
    }
  } else if (filenames[i] %in% c("tight_sd","tight_sd_fixed")) {
    print("tight_sd")
    parTab[parTab$names == "local_r_sd","values"] <- 0.5
    prior_func <- prior_func_rlocal
    if(filenames[i] == "tight_sd_fixed"){
      print("tight_sd_fixed")
      parTab[parTab$names == "t0","fixed"] <- 1
    }
  } else if (filenames[i] %in% c("unknown_sd","unknown_sd_fixed")){
    print("unknown_sd")
    parTab[parTab$names == "local_r_sd","fixed"] <- 0
    prior_func <- prior_func_rlocal
    if(filenames[i] == "unknown_sd_fixed"){
      print("unknown_sd_fixed")
      parTab[parTab$names == "t0","fixed"] <- 1
    }
  } else {
    print("unconstrained")
    prior_func <- prior_func_base
    if(filenames[i] == "unconstrained_fixed"){
      print("unconstrained_fixed")
      parTab[parTab$names == "t0","fixed"] <- 1
    }
  }
  startTab <- generate_start_tab(as.data.frame(parTab))
  
  ## MCMC
  ## Run first chain
  output <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars1, filename=filename_tmp,
                     CREATE_POSTERIOR_FUNC=create_model_func_provinces_fixed, mvrPars=NULL,
                     PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                     daily_import_probs = import_probs, daily_export_probs = export_probs,
                     incubation_ver="lnorm",
                     noise_ver="poisson",model_ver=2)
  
  
  ## Use this as input to multivariate chain
  chain <- read.csv(paste0(filename_tmp,"_univariate_chain.csv"))
  
  pdf(paste0(filename_tmp,"_chain.pdf"))
  plot(coda::as.mcmc(chain[,c("confirm_delay_shape","confirm_delay_scale","weibull_alpha","weibull_sigma","t0","growth_rate","K","lnlike")]))
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
                     CREATE_POSTERIOR_FUNC=create_model_func_provinces_fixed, mvrPars=mvrPars,
                     PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                     daily_import_probs = import_probs, daily_export_probs = export_probs,
                     incubation_ver="lnorm",
                     noise_ver="poisson",model_ver=2)
  
  ## Check convergence
  chain <- read.csv(output$file)
  pdf(paste0(filename_tmp,"_chain.pdf"))
  plot(coda::as.mcmc(chain[,c("confirm_delay_shape","confirm_delay_scale","weibull_alpha","weibull_sigma","t0","growth_rate","K","lnlike")]))
  dev.off()
}

