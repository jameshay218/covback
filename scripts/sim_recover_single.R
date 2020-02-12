setwd("~/Documents/covback")
#Rcpp::compileAttributes()
#devtools::document()
devtools::load_all()

library(lazymcmc)
library(tidyverse)

inc_period_draws <- read.csv("~/Documents/case_to_infection/data/backer_weibull_draws.csv",stringsAsFactors=FALSE)
parTab <- read.csv("pars/partab_single.csv",stringsAsFactors=FALSE)

## Make strong prior on alpha and sigma
prior_func <- create_incubation_prior(inc_period_draws)

## Generate some fake data
tmax <- 75
times <- 0:tmax
gamma_changes <- generate_confirmation_delays(5, 5, 15, 15, 0.1,0.1,250)
sim_dat <- simulate_observed_data_single(parTab, 100, treport=100, confirm_delay_pars=gamma_changes$pars,noise_ver = "poisson")
plot_simulations(sim_dat$aggregated)
dat1 <- sim_dat$aggregated %>% filter(var == "date_report_observable") %>% select(-var)

## Check that posterior solving works
f <- create_model_func(parTab,data=dat1,confirm_delay_pars=gamma_changes$pars,
                       PRIOR_FUNC=prior_func)
f(parTab$values)

## MCMC
## Run first chain
mcmcPars <- c("iterations"=10000,"popt"=0.44,"opt_freq"=250,
              "thin"=1,"adaptive_period"=5000,"save_block"=1000)
output <- run_MCMC(parTab=parTab, data=dat1, mcmcPars=mcmcPars, filename="test",
                   CREATE_POSTERIOR_FUNC=create_model_func, mvrPars=NULL,
                   PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                   confirm_delay_pars=gamma_changes$pars,ver="posterior")

## Use this as input to multivariate chain
chain <- read.csv(output$file)
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,0.1,w=0.8)

## Run second chain
mcmcPars <- c("iterations"=30000,"popt"=0.234,"opt_freq"=250,
              "thin"=1,"adaptive_period"=10000,"save_block"=1000)
output <- run_MCMC(parTab=parTab, data=dat1, mcmcPars=mcmcPars, filename="test",
                   CREATE_POSTERIOR_FUNC=create_model_func, mvrPars=mvrPars,
                   PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                   confirm_delay_pars=gamma_changes$pars,ver="posterior")


## Check convergence
chain <- read.csv(output$file)
pdf("tmp.pdf")
plot(coda::as.mcmc(chain[chain$sampno > 10000,]))
dev.off()

p_fit <- plot_model_fit(chain, parTab, sim_dat$aggregated, confirm_delay_pars=gamma_changes$pars)
png("fit_plot.png",height=5,width = 8,units="in",res=90)
p_fit
dev.off()

p_estimates <- plot_posteriors(chain[chain$sampno > 10000,],parTab)
png("sim_recover.png",height=5,width = 8,units="in",res=90)
p_estimates
dev.off()
