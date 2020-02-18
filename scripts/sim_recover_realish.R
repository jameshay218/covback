setwd("~/Documents/covback")
Rcpp::compileAttributes()
devtools::document()
devtools::load_all()

library(lazymcmc)
library(tidyverse)

tmax <- 75

## Incubation period draws and parameter values
inc_period_draws <- read.csv("~/Documents/case_to_infection/data/backer_weibull_draws.csv",stringsAsFactors=FALSE)
parTab <- read.csv("pars/partab_provinces.csv",stringsAsFactors=FALSE)

## Generate some fake travel probabilities
## Daily probability of leaving the seed province
## Assume exportations stop halfway through

# export_probs <- c(runif(54,0.01,0.05),rep(0,75-54+1))
# 
# ## Let's say 5 provinces - exportations get spread across these
# import_probs <- matrix(0, nrow=n_provinces-1, ncol=tmax+1)
# import_prob_averages <- runif(n_provinces-1, 0.1,0.8)
# for(i in 1:nrow(import_probs)){
#   import_probs[i,] <- runif(tmax+1, import_prob_averages[i]-0.05, import_prob_averages[i]+0.05)
# }
# ## Normalize
# import_probs_sums <- colSums(import_probs)
# import_probs <- t(apply(import_probs, 1, function(x) x/import_probs_sums))
# import_probs <- rbind(rep(1, tmax+1), import_probs)
## Real export probs
export_probs <- read.csv("data/export_probs.csv")
export_probs <- export_probs[1:60,]
tmax <- nrow(export_probs)-1
export_probs <- export_probs$prob_leaving

## Real import probs
import_probs <- read.csv("data/import_probs.csv",header = TRUE)
import_probs <- import_probs[1:60,]
import_probs <- t(import_probs[,2:ncol(import_probs)])
import_probs <- rbind(rep(1, ncol(import_probs)),import_probs)

times <- 0:tmax
parTab <- create_many_province_partab(parTab, 27, tmax)
provinces <- unique(parTab$province)
provinces <- provinces[provinces != "all"]
n_provinces <- length(provinces)
#parTab[parTab$province == "1" & parTab$names == "t0","values"] <- 0
#parTab[parTab$province == "1" & parTab$names == "growth_rate","values"] <- 0.2
#parTab[parTab$province != "1" & parTab$names == "growth_rate","values"] <- 0.001
## Make strong prior on alpha and sigma
prior_func <- create_incubation_prior(inc_period_draws)
#prior_func <- NULL
## Generate some fake data

#gamma_changes <- generate_confirmation_delays(3, 1, 5, 10, 0.05,0.05,250)
#confirm_delay_pars <- gamma_changes$pars
confirm_delay_pars <- read.csv("data/fitted_confirm_delays.csv")
confirm_delay_pars <- confirm_delay_pars %>%  select(shape, scale)
confirm_delay_pars$date_onset <- 0:(nrow(confirm_delay_pars)-1)
confirm_delay_pars <- confirm_delay_pars %>% filter(date_onset <= tmax)
plot_reporting_landscape(confirm_delay_pars$shape, confirm_delay_pars$scale)

sim_dat <- simulate_observed_data_provinces(parTab,tmax, tmax, confirm_delay_pars,
                                            import_probs, export_probs, TRUE,"poisson")

dat2 <- sim_dat$aggregated
dat2$date <- as.Date(dat2$date, origin="12/01/2019", format="%m/%d/%Y")
p_dat <- ggplot(dat2) + 
  geom_line(aes(x=date,y=n,col=var)) +
  #scale_y_log10() +
  scale_x_date(breaks="7 days") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  facet_wrap(~province,scales="free_y", ncol=5)# +
#  scale_x_date()

pdf("sim_dat.pdf",height=12,width=15)
p_dat
dev.off()

dat1 <- sim_dat$aggregated %>% filter(var=="date_report_observable") %>% select(-var)

## Check that posterior works
f <- create_model_func_provinces(parTab,dat1, confirm_delay_pars = confirm_delay_pars, 
                                 daily_import_probs = import_probs, daily_export_probs = export_probs,
                                 PRIOR_FUNC=prior_func)
Rprof(tmp<-tempfile())
f(parTab$values)
Rprof(NULL)
summaryRprof(tmp)

## Check that model works
f <- create_model_func_provinces(parTab,dat1, confirm_delay_pars = confirm_delay_pars, 
                                 daily_import_probs = import_probs, daily_export_probs = export_probs,
                                 PRIOR_FUNC=prior_func, ver-"model")
dat <- f(parTab$values)

startTab <- generate_start_tab(parTab)
startTab[startTab$names %in% c("weibull_alpha","weibull_sigma"),"values"] <- c(2.5,6)

## MCMC
## Run first chain
mcmcPars <- c("iterations"=50000,"popt"=0.44,"opt_freq"=5000,
              "thin"=10,"adaptive_period"=20000,"save_block"=100)
startTab[startTab$province=="all","fixed"] <- 1
output <- run_MCMC(parTab=startTab, data=dat1, mcmcPars=mcmcPars, filename="test",
                   CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=NULL,
                   PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                   confirm_delay_pars=confirm_delay_pars,
                   daily_import_probs = import_probs, daily_export_probs = export_probs,
                   ver="posterior")


## Use this as input to multivariate chain
chain <- read.csv(output$file)
best_pars <- get_best_pars(chain)
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,0.01,w=0.8)

## Start from best location of previous chain
startTab$values <- best_pars

## Run second chain
mcmcPars <- c("iterations"=50000,"popt"=0.234,"opt_freq"=5000,
              "thin"=10,"adaptive_period"=20000,"save_block"=100)
output <- run_MCMC(parTab=startTab, data=dat1, mcmcPars=mcmcPars, filename="test",
                   CREATE_POSTERIOR_FUNC=create_model_func_provinces, mvrPars=mvrPars,
                   PRIOR_FUNC = prior_func, OPT_TUNING=0.2,
                   confirm_delay_pars=gamma_changes$pars,
                   daily_import_probs = import_probs, daily_export_probs = export_probs,
                   ver="posterior")


## Check convergence
chain <- read.csv(output$file)
pdf("tmp.pdf")
plot(coda::as.mcmc(chain[chain$sampno > 10000,]))
dev.off()

chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]
#quants <- generate_prediction_intervals(chain, parTab, dat1, confirm_delay_pars=confirm_delay_pars,
#                                        daily_import_probs = import_probs, daily_export_probs = export_probs,
#                                        nsamp=100)

p <- plot_model_fit(chain, parTab, sim_dat$aggregated,confirm_delay_pars=confirm_delay_pars,
                    daily_import_probs = import_probs, daily_export_probs = export_probs,
                    nsamp=100)
pdf("tmp.pdf",height=12,width=10)
p + theme(legend.position=c(0.8,0.1))
dev.off()
