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

mcmcPars1 <- c("iterations"=80000,"popt"=0.44,"opt_freq"=1000,
              "thin"=10,"adaptive_period"=20000,"save_block"=1000)
mcmcPars2 <- c("iterations"=200000,"popt"=0.234,"opt_freq"=1000,
              "thin"=10,"adaptive_period"=50000,"save_block"=1000)
nchains <- 4
filename <- "nbinom_rerun"
## Set up parallelisation
n_clusters <- 1
#registerDoParallel(cores=n_clusters)
cl <- makeCluster(n_clusters)
registerDoParallel(cl)

tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")

## Incubation period draws and parameter values
inc_period_draws <- read.csv("data/backer_draws.csv",stringsAsFactors=FALSE)

## Real export probs
export_probs <- read_csv("data/export_probs_final.csv")$export_prob

## Real import probs
import_probs <- read_csv("data/import_probs_matched.csv")
import_probs <- as.matrix(import_probs[,2:ncol(import_probs)])

#parTab <- create_many_province_partab(parTab, length(use_provinces), tmax)
parTab <- read_csv("pars/partab_logistic_growth.csv")

parTab[parTab$names %in% c("confirm_delay_shape", "confirm_delay_scale"),"fixed"] <- 0
parTab[parTab$names %in% c("serial_interval_gamma_alpha", "serial_interval_gamma_scale"),"fixed"] <- 1
parTab[parTab$names == "local_r","upper_bound"] <- 10

provinces <- unique(parTab$province)
provinces <- provinces[provinces != "all"]
n_provinces <- length(provinces)

## Make strong prior on alpha and sigma
prior_func <- NULL
prior_func <- create_random_effects_rlocal()
## Generate some fake data

## Check that posterior works
confirmed_data1 <- as.data.frame(read_csv("data/confirmed_data_matched.csv"))
confirmed_data1 <- confirmed_data1 %>% mutate(n=ifelse(province==1, NA, n))
#confirmed_data1 <- confirmed_data1 %>% mutate(n=ifelse(province==25, NA, n))
f <- create_model_func_provinces_fixed(parTab,confirmed_data1, daily_import_probs = import_probs, 
                                 daily_export_probs = export_probs,
                                 PRIOR_FUNC=prior_func,ver="model",
                                 noise_ver="poisson",incubation_ver="lnorm",
                                 model_ver=2)
f(parTab$values)
setwd("../covback_chains_final")
res <- foreach(i=1:nchains,.packages=c("covback","lazymcmc","tidyverse")) %dopar% {
  filename_tmp <- paste0(filename, "_",i)
  parTab[parTab$names == "size","fixed"] <- 0
  parTab[parTab$names == "local_r_mean","values"] <- 1
  parTab[parTab$names == "local_r_mean","fixed"] <- 1
  parTab[parTab$names == "local_r_sd","fixed"] <- 1
  parTab[parTab$names %in% c("confirm_delay_shape","confirm_delay_scale"), "fixed"] <- 0
  parTab[parTab$names == "local_r","upper_bound"] <- 25
  parTab[parTab$names == "local_r","upper_start"] <- 10
  parTab[parTab$province == "3" & parTab$names == "local_r","fixed"] <- 0
  parTab[parTab$province == "3" & parTab$names == "local_r","values"] <- 0.4
  
  #prior_func <- NULL
  startTab <- generate_start_tab(as.data.frame(parTab))
  
  ## MCMC
  ## Run first chain
  output <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars1, filename=filename_tmp,
                     CREATE_POSTERIOR_FUNC=create_model_func_provinces_fixed, mvrPars=NULL,
                     PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                     daily_import_probs = import_probs, daily_export_probs = export_probs,
                     incubation_ver="lnorm",
                     noise_ver="nbinom",model_ver=2)
  
  
  ## Use this as input to multivariate chain
  chain <- read.csv(paste0(filename_tmp,"_univariate_chain.csv"))
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
                     noise_ver="nbinom",model_ver=2)
  
  
  ## Check convergence
  chain <- read.csv(output$file)
  pdf(paste0(filename_tmp,"_chain.pdf"))
  plot(coda::as.mcmc(chain[,c("confirm_delay_shape","confirm_delay_scale","weibull_alpha","weibull_sigma","t0","growth_rate","K","lnlike")]))
  dev.off()
}

