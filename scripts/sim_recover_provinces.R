setwd("~/GitHub/covback")
Rcpp::compileAttributes()
devtools::document()
devtools::load_all()

library(lazymcmc)
library(tidyverse)

inc_period_draws <- read.csv("~/Github/case_to_infection/data/backer_weibull_draws.csv",stringsAsFactors=FALSE)
parTab <- read.csv("pars/partab_provinces.csv",stringsAsFactors=FALSE)

## Make strong prior on alpha and sigma
prior_func <- create_incubation_prior(inc_period_draws)
#prior_func <- NULL
## Generate some fake data
tmax <- 75
times <- 0:tmax
gamma_changes <- generate_confirmation_delays(5, 5, 15, 15, 0.1,0.1,250)
confirm_delay_pars <- gamma_changes$pars
sim_dat <- simulate_observed_data_provinces(parTab,tmax, tmax, gamma_changes$par,TRUE,"poisson")
ggplot(sim_dat$aggregated) + 
  geom_line(aes(x=date,y=n,col=var)) +
  facet_wrap(~province,scales="free_y")
dat1 <- sim_dat$aggregated %>% filter(var=="date_report_observable") %>% select(-var)

## Check that posterior works
f <- create_model_func_provinces(parTab,dat1, confirm_delay_pars = confirm_delay_pars, PRIOR_FUNC=prior_func)
f(parTab$values)


## Check that posterior works
f <- create_model_func_provinces(parTab,data=NULL,tmax=tmax, confirm_delay_pars = confirm_delay_pars, PRIOR_FUNC=prior_func,ver="model")
f(parTab$values)


startTab <- generate_start_tab(parTab)

## MCMC
## Run first chain
mcmcPars <- c("iterations"=50000,"popt"=0.44,"opt_freq"=1000,
              "thin"=1,"adaptive_period"=20000,"save_block"=1000)
output <- run_MCMC(parTab=startTab, data=dat1, mcmcPars=mcmcPars, filename="test",
                   CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=NULL,
                   PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                   confirm_delay_pars=gamma_changes$pars,ver="posterior")


## Use this as input to multivariate chain
chain <- read.csv(output$file)
best_pars <- get_best_pars(chain)
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,0.01,w=0.8)

## Start from best location of previous chain
startTab$values <- best_pars

## Run second chain
mcmcPars <- c("iterations"=30000,"popt"=0.234,"opt_freq"=250,
              "thin"=1,"adaptive_period"=10000,"save_block"=1000)

output <- run_MCMC(parTab=startTab, data=dat1, mcmcPars=mcmcPars, filename="test",
                   CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=mvrPars,
                   PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                   confirm_delay_pars=gamma_changes$pars,ver="posterior")


## Check convergence
chain <- read.csv(output$file)
pdf("tmp.pdf")
plot(coda::as.mcmc(chain[chain$sampno > 10000,]))
dev.off()

quants <- generate_prediction_intervals(chain, parTab, dat1, confirm_delay_pars,nsamp=100)

plot_model_fit(chain, parTab, sim_dat$aggregated,confirm_delay_pars,nsamp=100)
