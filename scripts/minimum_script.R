library(grid)
library(lazymcmc)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(doParallel)

## Set up parallelisation
n_clusters <- 4
#registerDoParallel(cores=n_clusters)
cl <- makeCluster(n_clusters)
registerDoParallel(cl)

#library(covback)
setwd("~/Documents/GitHub/covback/")
install.packages("~/Documents/GitHub/covback/",repos=NULL,type="source")
library(covback)

test_pars <- read.csv("data/test_pars.csv",stringsAsFactors=FALSE)

#################################
## SETUP PARAMETERS FOR EVERYONE
#################################
mcmcPars1 <- c("iterations"=2000,"popt"=0.44,"opt_freq"=2000,
              "thin"=1,"adaptive_period"=1000,"save_block"=1000)
mcmcPars2 <- c("iterations"=2000,"popt"=0.234,"opt_freq"=1000,
              "thin"=1,"adaptive_period"=1000,"save_block"=1000)

## Get times to solve over
tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")

## Read in data to fit
confirmed_data <- read_csv("data/confirmed_data.csv")
confirmed_data <- confirmed_data %>% filter(country_region == "Mainland China")
confirmed_data$date <- match(confirmed_data$date, times)
confirmed_data <- confirmed_data %>% select(province, date, diff)
colnames(confirmed_data)[3] <- "n"
all_reports <- expand.grid(province=unique(confirmed_data$province),
                           date=match(times, times))
confirmed_data <- confirmed_data %>% right_join(all_reports)
confirmed_data$date <- confirmed_data$date - 1

## Incubation period draws and parameter values
inc_period_draws <- read.csv("data/backer_draws.csv",stringsAsFactors=FALSE)
parTab <- read.csv("data/partab_to_run.csv",stringsAsFactors=FALSE)

## Serial interval draws
## Use estimated serial interval from http://rs.yiigle.com/yufabiao/1183269.htm
## gamma with pars 5.23 and 0.87
serial_interval_par1 <- 5.23
serial_interval_par2 <- 0.87

## Real export probs
export_dat <- read_csv("data/export_probs_raw.csv")
export_dat$Date <- as.POSIXct(export_dat$Date,format="%m/%d/%Y", tz="UTC")


## Real import probs
import_probs <- read.csv("data/import_probs.csv",header = TRUE)
import_probs <- t(import_probs[,2:ncol(import_probs)])
import_probs <- rbind(rep(1, ncol(import_probs)),import_probs)
row.names(import_probs)[1] <- "Hubei"

use_provinces <- row.names(import_probs)
use_provinces <- intersect(use_provinces, unique(confirmed_data$province))

confirmed_data <- confirmed_data %>% filter(province %in% use_provinces)
import_probs <- import_probs[match(use_provinces, row.names(import_probs)),]
confirmed_data$province <- match(confirmed_data$province, use_provinces)
confirmed_data <- confirmed_data %>% arrange(province, date)

setwd("../chains")

#####################################
## WHERE THE LOOP WILL START
#####################################
res <- foreach(i=1:nrow(test_pars),.packages=c("covback","lazymcmc","tidyverse")) %dopar% {
#for(i in 1:nrow(test_pars)) {  
  t0 <- test_pars$t0[i]
  model_ver <- test_pars$model_ver[i]
  chain_no <- test_pars$chain_no[i]
  local_growth <- as.logical(test_pars$local_growth[i])
  wuhan_pop_ini <- 13972843
  total_travellers <- test_pars$total_travellers[i]
  
  export_probs <- create_export_prob_matrix(total_travellers, wuhan_pop_ini, export_dat, tmin, tmax)
  filename1 <- paste0("chain_start",t0,"_ver",model_ver,"_local",local_growth,"_travellers",total_travellers,"_chain",chain_no)
  parTab[parTab$names == "t0" & parTab$province == 1,"values"] <- t0
  
  ## Make strong prior on alpha and sigma
  prior_func <- create_prior_startdate(parTab, inc_period_draws, t0, 5)
  startTab <- generate_start_tab(parTab)
  startTab[startTab$names %in% c("weibull_alpha","weibull_sigma"),"values"] <- c(2.5,6)
  
  if(!local_growth) {
    startTab[startTab$names == "local_r","values"] <- 0
    startTab[startTab$names == "local_r","fixed"] <- 1
  }
  
  ## MCMC
  ## Run first chain
  
  confirmed_data1 <- confirmed_data %>% mutate(n=ifelse(province=="1", NA, n))
  output_univariate <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars1, 
                                filename=filename1,
                     CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=NULL,
                     PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                     confirm_delay_pars=NULL,
                     daily_import_probs = import_probs, daily_export_probs = export_probs,
                     ver="posterior",model_ver=model_ver)
  
  ## Use this as input to multivariate chain
  chain <- read.csv(output_univariate$file)
  best_pars <- get_best_pars(chain)
  chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
  covMat <- cov(chain)
  mvrPars <- list(covMat,0.5,w=0.8)
  
  ## Start from best location of previous chain
  startTab$values <- best_pars
  
  ## Run second chain
  output_multivariate <- run_MCMC(parTab=startTab, data=confirmed_data1, mcmcPars=mcmcPars2, 
                                  filename=filename1,
                     CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=mvrPars,
                     PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                     confirm_delay_pars=NULL,
                     daily_import_probs = import_probs, daily_export_probs = export_probs,
                     ver="posterior")
  chain <- read.csv(output_multivariate$file)
  chain <- chain[chain$sampno >= mcmcPars2["adaptive_period"],]
  pdf(paste0(filename1, ".pdf"))
  plot(coda::as.mcmc(chain))
  dev.off()
  
  #################
  ## DONE
  #################
}
