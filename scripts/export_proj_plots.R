library(grid)
library(lazymcmc)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(data.table)
#library(covback)
setwd("~/Documents/GitHub/covback/")
devtools::document()
devtools::load_all()

library(coda)

nchains <- 3
filename <- "nbinom_rerun"

mcmcPars2 <- c("iterations"=200000,"popt"=0.234,"opt_freq"=1000,
               "thin"=10,"adaptive_period"=100000,"save_block"=1000)
tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")

## Incubation period draws and parameter values
inc_period_draws <- NULL

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

parTab[parTab$names %in% c("confirm_delay_shape", "confirm_delay_scale"),"fixed"] <- 0
parTab[parTab$names %in% c("serial_interval_gamma_alpha", "serial_interval_gamma_scale"),"fixed"] <- 1
parTab[parTab$names == "local_r","upper_bound"] <- 10
parTab[parTab$names == "size","fixed"] <- 0
parTab[parTab$names == "local_r_mean","values"] <- 1
parTab[parTab$names == "local_r_mean","fixed"] <- 1
parTab[parTab$names == "local_r_sd","fixed"] <- 1
parTab[parTab$names %in% c("confirm_delay_shape","confirm_delay_scale"), "fixed"] <- 0
parTab[parTab$names == "local_r","upper_bound"] <- 25
parTab[parTab$names == "local_r","upper_start"] <- 10
parTab[parTab$province == "3" & parTab$names == "local_r","fixed"] <- 0
parTab[parTab$province == "3" & parTab$names == "local_r","values"] <- 0.4



#chain1 <- read.csv("unknown_confirm_delay_lowerR_nohei_1_multivariate_chain.csv")
#chain2 <- read.csv("unknown_confirm_delay_lowerR_nohei_2_multivariate_chain.csv")
#chain3 <- read.csv("unknown_confirm_delay_lowerR_nohei_3_multivariate_chain.csv")
chain1 <- read.csv("nbinom_rerun_1_multivariate_chain.csv")
chain2 <- read.csv("nbinom_rerun_2_multivariate_chain.csv")
chain3 <- read.csv("nbinom_rerun_3_multivariate_chain.csv")

setwd("chains")
chains <- list(chain1, chain2, chain3)

free_indices <- which(parTab$fixed == 0) + 1

wow <- lapply(chains, function(x){
  x <-x[x$sampno > mcmcPars2["adaptive_period"],c(1,free_indices)]
  x <- x[seq(1,nrow(x),by=10),]
  #x <- as.mcmc(x)
})
for(i in 1:length(wow)){
  wow[[i]]$chain <- i
}
wow <- do.call("rbind",wow)
wow$sampno_all <- 1:nrow(wow)
melt_wow <- reshape2::melt(wow,id.vars=c("sampno","chain"))
import_probs <- read.csv("../data/import_probs_matched.csv")
provinces <- as.character(import_probs[,1])


labels <- paste0("R[",provinces[2:length(provinces)],"]^{local}")
names(labels) <- paste0("local_r.",as.numeric(parTab[parTab$names == "local_r" &
                                                       parTab$province != 1,]$province)-1)

labels <- c(labels, c("K"="K","size"="phi", "confirm_delay_shape"="k^s", 
                      "confirm_delay_scale"="theta^s"))

labels

melt_wow$variable <- as.character(melt_wow$variable)
melt_wow <- melt_wow %>% mutate(group=ifelse(variable %like% "local_r","R","other"))
#melt_wow <- melt_wow[melt_wow$variable != "size",]
melt_wow$variable <- labels[melt_wow$variable]
melt_wow$variable <- factor(melt_wow$variable, levels=labels)
melt_wow <- melt_wow %>% drop_na()
min_samp <- min(melt_wow$sampno)
max_samp <- max(melt_wow$sampno)
p_r <- melt_wow %>% filter(group == "R") %>% 
  ggplot() + geom_line(aes(x=sampno, y=value,col=as.factor(chain))) + 
  facet_wrap(~variable,labeller=label_parsed, ncol=5) +
  xlab("Iteration") +
  ylab("")+
  scale_x_continuous(limits=c(min_samp, max_samp),breaks=seq(min_samp, max_samp, by=500000)) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6)) + labs(tag = "B")
p_other <- melt_wow %>% filter(group == "other") %>% 
  ggplot() + geom_line(aes(x=sampno, y=value,col=as.factor(chain))) + 
  facet_wrap(~variable,labeller=label_parsed, nrow=1,scales="free_y") +
  xlab("") +
  ylab("") +
  scale_x_continuous(limits=c(min_samp, max_samp),
                     breaks=seq(min_samp, max_samp, by=50000)) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6)) + labs(tag = "A")

p_all <- (p_other + p_r) + plot_layout(ncol=1,heights=c(0.2,1))
p_all
pdf("test_plot.pdf",height=8,width=8)
print(p_all)
dev.off()


p_r <- melt_wow %>% filter(group == "R") %>% 
  ggplot() + geom_density(aes(x=value,col=as.factor(chain))) + 
  facet_wrap(~variable,labeller=label_parsed, ncol=5,scales="free") +
  xlab("Iteration") +
  ylab("")+
  #scale_x_continuous(limits=c(min_samp, max_samp),breaks=seq(min_samp, max_samp, by=500000)) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6)) 
p_r
## Real import probs
import_probs <- read_csv("../data/import_probs_matched.csv")
import_probs <- as.matrix(import_probs[,2:ncol(import_probs)])
quants_summary <- generate_prediction_intervals(chain1, parTab, confirmed_data1, 
                                                daily_import_probs = import_probs, daily_export_probs = export_probs,
                                                nsamp=100,return_draws = FALSE,model_ver=2,noise_ver="poisson",
                                                incubation_ver="lnorm")
quants_summary$chain <- i
quants_summary$date <- as.Date(quants_summary$date, origin="2019-11-01")
quants_summary$province <- provinces[quants_summary$province]


factor_levels <- confirmed_data1 %>% filter(province != 1) %>% 
  group_by(province) %>%
  summarise(x=sum(n,na.rm=TRUE)) %>%
  arrange(-x) %>% pull(province)
factor_levels_names <- provinces[factor_levels]

quants_summary_hubei <- quants_summary %>% filter(province == "Hubei")
quants_summary_other <- quants_summary %>% filter(province != "Hubei")

quants_summary_other$province <- factor(quants_summary_other$province, levels=factor_levels_names)


var_key <- c("cantravel_prevalence"="Prevalence of pre-confirmation infections",
             "presymptomatic_prevalence"="Prevalence of pre-symptomatic infections")
tmp <- quants_summary_other[quants_summary_other$var %in% names(var_key),]
tmp$var <- var_key[tmp$var]
prev_all_p <- ggplot(tmp[tmp$date >= "2020-01-01",]) +
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
  geom_line(aes(x=date,y=median,col=var)) +
  #geom_point(data=confirmed_data2[confirmed_data1$date >= "2020-01-01",],aes(x=date,y=n),size=0.5) +
  scale_x_date(breaks="7 days") +
  scale_fill_manual(values=c("#E69F00","#0072B2"))+
  scale_color_manual(values=c("#E69F00","#0072B2"))+
  geom_vline(xintercept=as.Date("2020-01-23",origin="2019-11-01"),linetype="dashed") +
  theme_pubr() + 
  ylab("Daily prevalence (absolute numbers)") +
  xlab("Date") +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
        axis.text.y=element_text(size=7),
        axis.title = element_text(size=8),
        legend.text = element_text(size=7),
        legend.title = element_blank(),
        legend.position=c(0.7,0.025),
        panel.grid.minor=element_blank()) +
  facet_wrap(~province,ncol=5,scales="free_y")
pdf("prev_all_p.pdf",height=7,width=8)
plot(prev_all_p)
dev.off()


tmp_hubei <- quants_summary_hubei %>% filter(date >= "2019-12-08" & var %in% names(var_key))
tmp_hubei$var <- var_key[tmp_hubei$var]
hubei_plot <- tmp_hubei %>%
  ggplot() + 
  geom_vline(xintercept=(as.Date("2020-01-23", format="%Y-%m-%d", tz="LMT", origin="2019-11-01")),
             linetype="dashed") +
  geom_line(aes(x=date,y=median,col=var)) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
  scale_fill_manual(values=c("#E69F00","#0072B2"))+
  scale_color_manual(values=c("#E69F00","#0072B2"))+
  scale_x_date(breaks="7 days") + 
  scale_y_continuous(limits=c(0,200000),breaks=seq(0,200000,by=25000)) +
  theme_bw() +
  ylab("Daily prevalence (absolute numbers)") +
  xlab("Date") +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
        axis.text.y=element_text(size=7),
        axis.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_blank(),
        legend.position="bottom",
        panel.grid.minor=element_blank())
hubei_plot

tmp_hubei2 <- quants_summary_hubei %>% filter(date >= "2019-12-08" & date <= "2020-01-01" & var %in% names(var_key))
tmp_hubei2$var <- var_key[tmp_hubei2$var]
hubei_plot_inset <- tmp_hubei2 %>%
  ggplot() + 
  geom_vline(xintercept=(as.Date("2020-01-23", format="%Y-%m-%d", tz="LMT", origin="2019-11-01")),linetype="dashed") +
  geom_line(aes(x=date,y=median,col=var)) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
  scale_fill_manual(values=c("#E69F00","#0072B2"))+
  scale_color_manual(values=c("#E69F00","#0072B2"))+
  scale_x_date(breaks="2 days") + 
  scale_y_continuous(limits=c(0,1400),breaks=seq(0,1400,by=200)) +
  theme_bw() +
  ylab("") +
  xlab("") +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=6),
        panel.background = element_blank(),
        axis.text.y=element_text(size=6),
        legend.position="none",
        panel.grid.minor=element_blank())
vp <- viewport(width=0.4,height=0.5, x=0.28,y=0.73)

pdf("hubei_prevalence.pdf",height=6,width=8)
print(hubei_plot)
print(hubei_plot_inset, vp=vp)
dev.off()

















quants <- generate_prediction_intervals(chain1, parTab, confirmed_data1, 
                                        daily_import_probs = import_probs, daily_export_probs = export_probs,
                                        nsamp=1000,return_draws = TRUE,model_ver=2,noise_ver="poisson",
                                        incubation_ver="lnorm")

samps <- quants$samp_ids
samps <- sort(samps)
chain_save <- chain1[chain1$sampno %in% samps,]
chain_save$sampno <- match(chain_save$sampno, samps)
chain_save$chain <- i

quants <- quants$draws
quants$chain <- i
quants$province <- provinces[quants$province]
#quants$province <- row.names(import_probs)[quants$province]
quants$date <- as.Date(quants$date, origin="2019-11-01")

write_csv(quants_summary, "prevalence_estimates_summary.csv")
write_csv(quants, "prevalence_estimates_draws.csv")
write_csv(chain_save, "mcmc_chain_thinned_logistic_growth.csv")


## 444 cases confirmed by 23rd Jan as per JHU data
number_conf <- quants %>% filter(province == "Hubei" & date <= "2020-01-23" & var == "observations") %>%
  group_by(sampno, chain) %>% 
  summarise(x=sum(n)) %>% pull(x)
mean(444/number_conf)
quantile(444/number_conf, c(0.025,0.5,0.975))

## Or 30th, 4903
number_conf2 <- quants %>% filter(province == "Hubei" & date <= "2020-01-30" & var == "observations") %>%
  group_by(sampno, chain) %>% 
  summarise(x=sum(n)) %>% pull(x)
mean(4903/number_conf2)
quantile(4903/number_conf2, c(0.025,0.5,0.975))

inc_p <- plot_incubation_period(chain1,nsamp=200,xmax=40, prior_pars=inc_period_draws) + ylab("Probability density") + xlab("Delay from infection to symptom onset")
conf_p <- plot_confirm_delay(chain1,nsamp=100,xmax=40) + ylab("Probability density") + xlab("Delay from symptom onset to confirmation")
delay_plots <- inc_p / conf_p
#png("delay_plots.png",height=5,width=7,res=300,units="in")
#delay_plots
#dev.off()
