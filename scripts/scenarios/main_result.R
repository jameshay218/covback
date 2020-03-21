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
devtools::document()
devtools::load_all()

mcmcPars1 <- c("iterations"=100000,"popt"=0.44,"opt_freq"=1000,
              "thin"=10,"adaptive_period"=50000,"save_block"=1000)
mcmcPars2 <- c("iterations"=350000,"popt"=0.234,"opt_freq"=1000,
              "thin"=10,"adaptive_period"=150000,"save_block"=1000)

nchains <- 4
filename1 <- "main_free_switch"
## Set up parallelisation
n_clusters <- 4
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

## Real import probs
import_probs <- read_csv("data/import_probs_matched.csv")
import_probs <- as.matrix(import_probs[,2:ncol(import_probs)])
colnames(import_probs) <- NULL

#parTab <- create_many_province_partab(parTab, length(use_provinces), tmax)
parTab <- read_csv("pars/partab_logistic_growth.csv")

provinces <- unique(parTab$province)
provinces <- provinces[provinces != "all"]
n_provinces <- length(provinces)

## Make strong prior on alpha and sigma
prior_func <- create_random_effects_rlocal(parTab,province_prior="2",r_local_mean=0.4,
                                           r_local_sd=0.05,func_ver="gamma")
#prior_func <- create_random_effects_rlocal2(parTab,r_local_mean=1,r_local_sd=0.5, prior_ver="gamma")
## Generate some fake data

## Check that posterior works
confirmed_data1 <- as.data.frame(read_csv("data/real/confirmed_data_matched_final.csv"))
confirmed_data1 <- confirmed_data1 %>% select(-province_name)
confirmed_data1 <- confirmed_data1 %>% mutate(n=ifelse(province==1, NA, n))
parTab[parTab$names == "t_switch","fixed"] <- 1
parTab[parTab$names == "t_switch","lower_bound"] <- 46
parTab[parTab$names == "t_switch","values"] <- 50
parTab[parTab$names == "t_switch","upper_bound"] <- 53

res <- foreach(i=1:nchains,.packages=c("covback","lazymcmc","tidyverse")) %dopar% {
  setwd("../covback_chains_final/tswitch")
  filename_tmp <- paste0(filename1, "_",i)
  parTab[parTab$names == "size","fixed"] <- 1
  parTab[parTab$names %in% c("confirm_delay_shape","confirm_delay_scale"), "fixed"] <- 0
  
  #prior_func <- NULL
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
                     CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=mvrPars,
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
