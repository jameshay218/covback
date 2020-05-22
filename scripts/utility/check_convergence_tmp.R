setwd("~/Documents/GitHub/covback_chains_final/main_results_rerun/")
chain <- read.csv("free_tswitch_1_multivariate_chain.csv")
chain$growth_rate <- log(exp(chain$K)-1)/(chain$t_switch - chain$t0 - chain$incu_mode)
chain$K <- exp(chain$K)
plot(coda::as.mcmc(chain[chain$sampno > 2000,]))

library(data.table)

setwd("~/Documents/GitHub/covback_chains_final/testing_still_works//")
chain <- read.csv("relaxed_sd_1_multivariate_chain.csv")
prior_func_base <- function(pars){
  names(pars) <- parTab$names
  t_switch <- log(pars["K"]-1)/pars["growth_rate"]  + pars[which(names(pars)=="t0")[1]] + 5
  prior_prob <- dnorm(t_switch, 86, 2, 1)
}

chain$t_switch <- log(chain$K-1)/chain$growth_rate  + chain$t0 + 5
chain$prior_prob <- dnorm(chain$t_switch, 86, 2, 1)
which_use <- colnames(chain)[colnames(chain) %like% "local_r."]

chain$prior_rlocal <- apply(chain[,which_use[4:length(which_use)]], 1, function(x) sum(dnorm(x, 1, 0.5, 1)))
chain$likelihood <- chain$lnlike - chain$prior_prob - chain$prior_rlocal
plot(coda::as.mcmc(chain[chain$sampno > 200000,]))

library(covback)

gamma_pars_from_mean_sd(9.5,0.131^2)
hist(rgamma(1000, shape=5259, scale=0.00181))


chain <- read.csv("~/Documents/GitHub/covback_chains_final/sim_recovery_single/data_2_chain_1_multivariate_chain.csv")
plot(coda::as.mcmc(chain[chain$sampno > 3000,]))


library(odin)
library(reshape2)

seir_ode <- odin::odin({
  ## Derivatives
  deriv(S) <- -b*S*I
  deriv(E) <- b*S*I - s*E
  deriv(I) <- s*E-g*I
  deriv(R) <- g*I
  deriv(inc) <- b*S*I
  
  ## Initial conditions
  initial(S) <- 0.99
  initial(E) <- 0
  initial(I) <- 0.01
  initial(R) <- 0.00
  initial(inc) <- 0
  
  ## parameters
  b <- 0.1
  g <- 0.05
  s <- 0.2
  config(base) <- "seir"
})

seir_mod <- seir_ode()
times <- seq(0,365,length.out=366)
seir_out <- seir_mod$run(times)

microbenchmark::microbenchmark(seir_mod$run(times))


seir_out_long <- melt(as.data.frame(seir_out),"t")
ggplot(seir_out_long,aes(x=t,y=value,colour=variable,group=variable))+
  # Add line
  geom_line(lwd=2)+
  #Add labels
  xlab("Time")+ylab("Number")
