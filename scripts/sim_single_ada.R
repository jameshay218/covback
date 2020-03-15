chain <- read.csv("mcmc_chain_thinned_logistic_growth.csv", stringsAsFactors = FALSE)
best_pars <- get_best_pars(chain)
par_names <- c("shape", "scale", "size", "weibull_alpha", "weibull_sigma", "K")
startTab_back_calc <- read.csv("data/startTab_back_calc.csv", stringsAsFactors = FALSE)
startTab_back_calc[match(par_names, startTab_back_calc$names), "values"] <- best_pars[par_names]
f <- create_model_func_provinces(startTab_back_calc, PRIOR_FUNC=prior_func, ver="model",model_ver=2, tmax = 123)
pars <- startTab_back_calc$values
sim_no_noise <- f(pars)
ggplot(sim_no_noise) + geom_line(aes(x=date,y=n,col=var))
set.seed(1)
confirmations_noise <- sim_no_noise %>% 
  filter(var == "confirmations") %>% 
  mutate(x = vapply(n, function(x) rpois(1, x), numeric(1))) %>%
  select(-var, -n) %>%
  rename(n = x)
write.csv(confirmations_noise, file = "data/sim_data.csv", row.names = FALSE)
