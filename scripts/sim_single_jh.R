setwd("/Users/james/Documents/GitHub/covback")
parTab <- read.csv("pars/partab_logistic_growth.csv",stringsAsFactors = FALSE)
parTab <- parTab[parTab$province %in% c("all","1"),]



f <- create_model_func_provinces(parTab, PRIOR_FUNC=NULL, ver="model",model_ver=2, tmax = 123)
pars <- parTab$values
names(pars) <- parTab$names

pars["t_switch"] <- 100
pars["t0"] <- 0
pars["K"] <- 200000
sim_no_noise <- f(pars)

set.seed(1)
confirmations_noise <- sim_no_noise %>% 
  filter(var == "confirmations") %>% 
  mutate(n = vapply(n, function(x) rpois(1, x), numeric(1)))  %>%
  mutate(var="observations")

all_dat <- bind_rows(confirmations_noise, sim_no_noise)



ggplot(all_dat) + geom_line(aes(x=date,y=n,col=var)) +
  scale_y_continuous(limits=c(0,20000))

