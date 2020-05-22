library(grid)
library(lazymcmc)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(doParallel)
setwd("~/Documents/GitHub/covback")
#setwd("~/git_repos/covback")
Rcpp::compileAttributes()
#devtools::document()
devtools::load_all()
#library(covback)
#set.seed(2)

tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")

## Set up parallelisation
#n_clusters <- 3
#cl <- makeCluster(n_clusters)
#registerDoParallel(cl)

## Simulate data
N <- 1000000
parTab <- read.csv("pars/partab_single.csv",stringsAsFactors = FALSE)

parTab[parTab$names %in% c("confirm_delay_shape","confirm_delay_scale"),"fixed"] <- 0
parTab[parTab$names == "growth_rate","values"] <- 0.15
parTab[parTab$names == "K","upper_bound"] <- N
parTab[parTab$names == "K","values"] <- N*0.6
parTab[parTab$names == "K","fixed"] <- 1
parTab[parTab$names == "reporting_propn","fixed"] <- 1
parTab[parTab$names == "reporting_propn","values"] <- 1
pars <- parTab$values
names(pars) <- parTab$names
t_cutoff <- 365

dat <- simulate_observed_data_single(pars, 365, t_cutoff, TRUE, "poisson",N,N, model_ver = "seir")
confirmed_data <- dat$aggregated %>% filter(var == "date_confirm_observable") %>% select(date, n) %>% mutate(province="1")
confirmed_data <- confirmed_data %>% mutate(n = ifelse(date <= t_cutoff, n, NA))
p1 <- ggplot() + 
  geom_line(data=confirmed_data[confirmed_data$date <= 30,], aes(x=date,y=n),col="red") +
  scale_x_continuous(limits=c(0,200)) +
  scale_y_continuous(limits=c(0,2500))
p2 <- ggplot() +
  geom_line(data=confirmed_data[confirmed_data$date <= 60,], aes(x=date,y=n),col="blue")+
  scale_x_continuous(limits=c(0,200)) +
  scale_y_continuous(limits=c(0,2500))
p3 <- ggplot() +
  geom_line(data=confirmed_data[confirmed_data$date <= 90,], aes(x=date,y=n),col="green")+
  scale_x_continuous(limits=c(0,200)) +
  scale_y_continuous(limits=c(0,2500))
p1 / p2 / p3


startTab <- generate_start_tab(parTab)
startTab[startTab$names %in% c("confirm_delay_shape","confirm_delay_scale"),"values"] <- c(6,2)
## Make strong prior on alpha and sigma

#prior_func <- NULL
## MCMC

f <- create_model_func(parTab, confirmed_data, ver="model",model_ver="seir",incubation_ver = "lnorm")
f(parTab$values) %>% ggplot() + geom_line(aes(x=date,y=n,col=var))

# there is only one parameter to be fitted so it's a short chain
mcmcPars <- c("iterations"=50000,"popt"=0.234,"opt_freq"=1000,
              "thin"=10,"adaptive_period"=50000,"save_block"=1000)

# create_model_func_provinces is a closure that creates a function to evaluate the likelihood
# model_ver = 2 is logistic growth
#set.seed(1)
f <- create_model_func(parTab, confirmed_data, NULL, ver="posterior",incubation_ver="lnorm",model_ver=2)
f(parTab$values)
covMat <- diag(x=1:nrow(parTab)*parTab$values)
mvrPars <- list(covMat,0.01,w=0.9,0.8)

beta_mean <- 0.2
beta_var <- 0.005

beta_alpha <- (((1-beta_mean)/beta_var) - 1/beta_mean)*beta_mean^2
beta_beta <- beta_alpha*(1/beta_mean -1)

hist(rbeta(10000, beta_alpha, beta_beta))

prior_func <- function(pars){
  names(pars) <- parTab$names
  #prior_prob <- dbeta(pars["reporting_propn"],beta_alpha, beta_beta, 1)
  #prior_prob <- 0
  lik1 <- dnorm(pars["confirm_delay_shape"],6, 0.5,1)
  lik2 <- dnorm(pars["confirm_delay_scale"],2,0.5,1)
  
  lik3 <- dnorm(pars["gamma"],5, 2, 1)
  lik4 <- dnorm(pars["sigma"],5,2,1)
  lik5 <- dnorm(pars["R0"],2.2,2,1)
  
  return(lik1+ lik2 + lik3 + lik4 + lik5)
}
setwd("~/Documents/GitHub/covback_chains_final/sim_recovery_single/")

confirmed_data1 <- confirmed_data %>% mutate(n=ifelse(date >= 30, NA, n))
confirmed_data2 <- confirmed_data %>% mutate(n=ifelse(date >= 60, NA, n))
confirmed_data3 <- confirmed_data %>% mutate(n=ifelse(date >= 90, NA, n))

f <- create_model_func(startTab, confirmed_data1, prior_func, ver="posterior",incubation_ver="lnorm",model_ver="seir")
f(startTab$values)


#res <- foreach(i=1:3,.packages=c("covback","lazymcmc","tidyverse")) %dopar% {
for(i in 3:3){
  setwd("~/Documents/GitHub/covback_chains_final/sim_recovery_single/")
  if(i == 1){
    dat <- confirmed_data1
  } else if(i == 2){
    dat <- confirmed_data2
  } else {
    dat <- confirmed_data3
  }
  for(j in 1:1){
    output <- run_MCMC(parTab=startTab, data=dat, mcmcPars=mcmcPars, filename=paste0("data_",i,"_chain_",j),
                     CREATE_POSTERIOR_FUNC=create_model_func, mvrPars=mvrPars,
                     PRIOR_FUNC = prior_func, OPT_TUNING=0.2,incubation_ver="lnorm",
                     ver="posterior",model_ver="seir")
  }
}
chain <- read.csv(output$file)
chain <- read.csv("data_1_chain_1_multivariate_chain.csv")
chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]
#plot(coda::as.mcmc(chain))
chain$R0_less <- 1
chain$switch_time <- 400
quants <- generate_prediction_intervals(chain, parTab, confirmed_data,daily_import_probs = NULL, daily_export_probs=NULL, nsamp=1000, 
                                        incubation_ver="lnorm",model_ver="seir",single_province=TRUE)

chain$switch_time <- 100
quants1 <- generate_prediction_intervals(chain, parTab, confirmed_data,daily_import_probs = NULL, daily_export_probs=NULL, nsamp=1000, 
                                         incubation_ver="lnorm",model_ver="seir",single_province=TRUE)


quants %>% filter(var %in% c("observations","infections")) %>% 
  ggplot() + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
  geom_ribbon(data=quants1[quants1$var %in% c("infections","observations"),],aes(x=date,ymin=lower,ymax=upper,fill=var), alpha=0.2) +
  geom_line(data=quants1[quants1$var %in% c("infections","observations"),],aes(x=date,y=median,col=var)) +
  geom_vline(xintercept=100,linetype="dashed") +
  geom_vline(xintercept=100+18,linetype="dashed")+
  #coord_cartesian(ylim=c(0,100)) +
  geom_line(aes(x=date,y=median,col=var)) +
  geom_point(data=confirmed_data3,aes(x=date,y=n))
