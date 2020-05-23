library(grid)
library(lazymcmc)
library(tidyverse)
library(ggpubr)
library(patchwork)
#library(covback)
library(data.table)
library(coda)
library(doParallel)
library(data.table)
library(ggsci)
library(ggthemes)
#n_clusters <- 12
#library(covback)
#cl <- makeCluster(n_clusters)
#registerDoParallel(cl)
export_theme <- theme_tufte() + 
  theme(
    ## Axis text and titles
    axis.text.x = element_text(size=8,angle = 45, hjust = 1,family="sans"),
    axis.text.y=element_text(size=8,family="sans"),
    axis.title.x=element_text(size=10,family="sans",vjust=-1),
    axis.title.y=element_text(size=10,family="sans"),
    ## Axis lines
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),
    ## Legends
    legend.title=element_text(size=10,family="sans",face="italic"),
    legend.text=element_text(size=10,family="sans"),
    legend.key.size= unit(0.2, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    ## Strips for facet_wrap
    strip.text=element_text(size=8,family="sans"),
    strip.background=element_rect(fill="#f0f0f0"),
    ## Overall plot
    panel.grid.major=element_line(size=0.25,color="#f0f0f0"),
    plot.tag=element_text(size=12,family="sans",face="bold")
    
  )
cols <- ggsci::pal_npg()(10)

setwd("~/Documents/GitHub/covback/")
devtools::load_all()
savewd <- "~/Google Drive/nCoV/backcalculation_paper/figures_final_adjusted/"
adaptive_period <- 100000
nsamp <- 1000
scale_reporting <- TRUE

n_wuhan <- 9785388

#chain_top_wd <- "~/Documents/GitHub/covback_chains_final/main_results_final_maybe/"
chain_top_wd <- "~/Documents/GitHub/covback_chains_final/final_20200522/"

run_names <- list.files(chain_top_wd)
scenario_key <- read_csv("~/Documents/GitHub/covback/scripts/scenario_key.csv")

## Real export probs
export_probs <- read_csv("data/export_probs_matched.csv")$export_prob
export_probs_lower <- read_csv("data/export_probs_lower.csv")$export_prob

tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")

time_varying_report_pars <- data.frame(date=as.Date(times,origin="2019-11-01"),shape=3.18,scale=1/0.59)
time_varying_report_pars[time_varying_report_pars$date <= as.Date("2020-01-27",origin="2019-11-01"),"shape"] <- 3.72
time_varying_report_pars[time_varying_report_pars$date <= as.Date("2020-01-27",origin="2019-11-01"),"scale"] <- 1/0.42


metrics_23rdJan <- NULL
metrics_30thJan <- NULL
metrics_cumu_incidence <- NULL
metrics_peak_times <- NULL
metrics_prevalence_23rdJan <- NULL
metrics_prevalence_30thJan <- NULL

parameter_estimates_all <- NULL


#run_names <- c("diffuse_r_local","main","t_switch_6")
#run_names <- c("main")
#res <- foreach(i=1:length(run_names),.packages=c("covback","lazymcmc","tidyverse",
#                                                 "data.table","ggpubr","patchwork",
#                                                 "coda","grid")) %dopar% {
for(i in seq_along(run_names)){
  runname <- run_names[i]
  print(runname)
  filename <- runname
  dir.create(paste0(savewd,"/", filename))
  filename <- paste0(filename,"/",filename)

    
  #chain_wd <- paste0("~/Documents/GitHub/covback_chains_final/main_results_final_maybe/",runname)
  chain_wd <- paste0("~/Documents/GitHub/covback_chains_final/final_20200522/",runname)
  
  parTab <- read_csv("pars/partab_logistic_growth.csv")
  parTab[parTab$names == "t0","fixed"] <- 0
  
  parTab[parTab$names == "K","values"] <- log(parTab[parTab$names == "K","values"])
  parTab[parTab$names == "K","upper_bound"] <- log(parTab[parTab$names == "K","upper_bound"])
  parTab[parTab$names == "K","lower_bound"] <- log(parTab[parTab$names == "K","lower_bound"])
  parTab[parTab$names == "K","lower_start"] <- log(parTab[parTab$names == "K","lower_start"])
  parTab[parTab$names == "K","upper_start"] <- log(parTab[parTab$names == "K","upper_start"])
  parTab[parTab$names == "t0","values"] <- scenario_key$t0_val[i]
  parTab[parTab$names == "t0","fixed"] <- scenario_key$t0_fixed[i]
  parTab[parTab$names == "local_r_sd","values"] <- scenario_key$r_local_sd[i]
  parTab[parTab$names == "t_switch","fixed"] <- scenario_key$t_switch_fixed[i]
  parTab[parTab$names == "t_switch","values"] <- scenario_key$t_switch_val[i]
  
  prior_func_rlocal <- function(pars){
    r_local_sd1 <- pars["local_r_sd"]
    r_local_mean1 <- pars["local_r_mean"]
    gamma_pars <- gamma_pars_from_mean_sd(r_local_mean1, r_local_sd1^2)
    r_locals <- pars[which(names(pars) == "local_r")]
    lik <- dgamma(r_locals[r_index],gamma_pars[[1]],scale=gamma_pars[[2]],log=TRUE)
    return(sum(lik))
  }
  
  if(scenario_key$travellers[i] == 1){
    export_probs_use <- export_probs
  } else {
    export_probs_use <- export_probs_lower
  }
  

  chains <- load_mcmc_chains(chain_wd,parTab,unfixed=TRUE,thin=1,burnin=adaptive_period,multi=TRUE)
  chains <- as.list(chains[[1]])
  for(i in 1:length(chains)){
    chains[[i]] <- as.data.frame(chains[[i]])
    chains[[i]]$sampno <- 1:nrow(chains[[i]])
    chains[[i]]$chain <- i
  }
  chains <- do.call("rbind", chains)
  
  melted_chains <- reshape2::melt(chains,id.vars=c("sampno","chain"))
  
  import_probs <- read.csv("data/import_probs_matched.csv")
  provinces <- as.character(import_probs[,1])
  provinces[which(provinces == "Inner Mongolia")] <- "Inner_Mongolia"
  
  labels <- paste0("R[",provinces[2:length(provinces)],"]^{local}")
  names(labels) <- c("local_r",paste0("local_r.",1:(length(provinces)-2)))
  
  labels <- c(labels, c("K"="log(K)","size"="phi", "serial_interval_mean"="SI_mean", 
                        "serial_interval_var"="SI_var"))
  
  
  melted_chains$variable <- as.character(melted_chains$variable)
  melted_chains <- melted_chains %>% mutate(group=ifelse(variable %like% "local_r","R","other"))
  melted_chains$variable <- labels[melted_chains$variable]
  melted_chains$variable <- factor(melted_chains$variable, levels=labels)
  melted_chains <- melted_chains %>% drop_na()
  min_samp <- min(melted_chains$sampno)
  max_samp <- max(melted_chains$sampno)
  p_r <- melted_chains %>% filter(group == "R") %>% 
    ggplot() + geom_line(aes(x=sampno, y=value,col=as.factor(chain))) + 
    facet_wrap(~variable,labeller=label_parsed, ncol=5) +
    xlab("Sample") +
    ylab("")+
    scale_x_continuous(limits=c(min_samp, max_samp),breaks=seq(min_samp, max_samp, by=500)) +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x=element_text(size=6),
          axis.text.y=element_text(size=6)) + labs(tag = "B")
  p_other <- melted_chains %>% filter(group == "other") %>% 
    ggplot() + geom_line(aes(x=sampno, y=value,col=as.factor(chain))) + 
    facet_wrap(~variable,labeller=label_parsed, nrow=1,scales="free_y") +
    xlab("") +
    ylab("") +
    scale_x_continuous(limits=c(min_samp, max_samp),
                       breaks=seq(min_samp, max_samp, by=500)) +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x=element_text(size=6),
          axis.text.y=element_text(size=6)) + labs(tag = "A")
  
  p_all <- (p_other + p_r) + plot_layout(ncol=1,heights=c(0.2,1))
  pdf(paste0(savewd, "/", filename,"_convergence.pdf"),height=8,width=8)
  print(p_all)
  dev.off()
  
  temp_estimates <- melted_chains %>% ungroup() %>%
    group_by(variable) %>%
    summarise(mean=median(value),
           lower=quantile(value, 0.025),
           upper=quantile(value,0.975)) %>%
    mutate(estimate=paste0(signif(mean,3), " (", signif(lower,3),"-", signif(upper,3), ")")) %>%
    select(variable, estimate) %>%
    mutate(Scenario=runname)
  parameter_estimates_all <- rbind(parameter_estimates_all, temp_estimates)
  
  
  png(paste0(savewd, "/", filename,"_convergence.png"),height=8,width=8,units="in",res=300)
  print(p_all)
  dev.off()
  
  chains2 <- as.data.frame(load_mcmc_chains(chain_wd,parTab,unfixed=FALSE,thin=1,burnin=adaptive_period,multi=TRUE)[[2]])
  chains2$sampno <- 1:nrow(chains2)
  
  ## Real import probs
  import_probs <- read_csv("data/import_probs_matched.csv")
  import_probs <- as.matrix(import_probs[,2:ncol(import_probs)])
  
  confirmed_data1 <- as.data.frame(read_csv("data/real/midas_data_final.csv"))
  confirmed_data1 <- confirmed_data1 %>% mutate(n=ifelse(province==1, NA, n))
  confirmed_data1 <- confirmed_data1 %>% select(-province_raw)
  print("Generating prediction intervals 1")
  quants_summary <- generate_prediction_intervals(chains2, parTab, confirmed_data1, 
                                                  daily_import_probs = import_probs, daily_export_probs = export_probs_use,
                                                  time_varying_confirm_delay_pars = time_varying_report_pars,
                                                  nsamp=nsamp,return_draws = FALSE,model_ver="logistic",noise_ver="poisson",
                                                  incubation_ver="lnorm",scale_reporting=FALSE, 
                                                  report_rate_switch=83,report_rate_1=0.14,report_rate_2=0.65)
  print("... done")
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
  
  var_key <- c("presymptomatic_prevalence"="Prevalence of pre-symptomatic infections",
               "cantravel_prevalence"="Prevalence of infections eligible\n for international travel",
               "symptomatic_prevalence"="Prevalence of symptomatic, \nnot yet recovered infections",
               "disease_prevalence"="Prevalence of all infected individuals"
               )
  tmp <- quants_summary_other[quants_summary_other$var %in% names(var_key),]
  tmp$var <- var_key[tmp$var]
  prev_all_p <- ggplot(tmp[tmp$date >= "2020-01-01",]) +
    geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
    geom_line(aes(x=date,y=median,col=var)) +
    #geom_point(data=confirmed_data2[confirmed_data1$date >= "2020-01-01",],aes(x=date,y=n),size=0.5) +
    scale_x_date(breaks="7 days") +
    scale_fill_manual(values=c("#E69F00","#0072B2","#009E73","#CC79A7"))+
    scale_color_manual(values=c("#E69F00","#0072B2","#009E73","#CC79A7"))+
    geom_vline(xintercept=as.Date("2020-01-23",origin="2019-11-01"),linetype="dashed") +
    theme_pubr() + 
    ylab("Daily prevalence (absolute numbers)") +
    xlab("Date") +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
          axis.text.y=element_text(size=7),
          strip.text=element_text(size=8),
          axis.title = element_text(size=8),
          legend.text = element_text(size=7),
          legend.title = element_blank(),
          legend.position=c(0.6,0),
          legend.direction = "horizontal",
          panel.grid.minor=element_blank()) +
    facet_wrap(~province,ncol=4,scales="free_y")
  
  write_csv(tmp,paste0(savewd,"/",filename,"_all_prevalence.csv"))
  
  pdf(paste0(savewd,"/",filename,"_prev_all_p.pdf"),height=8,width=8)
  plot(prev_all_p)
  dev.off()
  png(paste0(savewd,"/",filename,"_prev_all_p.png"),height=8,width=8,res=300,units="in")
  plot(prev_all_p)
  dev.off()
  
  pop_wuhan <- 9785388
  
  tmp_hubei <- quants_summary_hubei %>% filter(date >= "2019-12-08" & var %in% names(var_key))
  tmp_hubei$var <- var_key[tmp_hubei$var]
  hubei_plot <- tmp_hubei %>%
    ggplot() + 
    geom_rect(xmin=as.Date("2020-01-23"),xmax=as.Date("2021-01-23"),ymin=0,ymax=1,fill="grey80",alpha=0.5)+
    geom_vline(xintercept=(as.Date("2020-01-23", format="%Y-%m-%d", tz="LMT", origin="2019-11-01")),
               linetype="dashed") +
    geom_line(aes(x=date,y=median/pop_wuhan,col=var)) + 
    geom_ribbon(aes(x=date,ymin=lower/pop_wuhan,ymax=upper/pop_wuhan,fill=var),alpha=0.25) +
    scale_fill_manual(values=c("#E69F00","#0072B2","#009E73","#CC79A7"))+
    scale_color_manual(values=c("#E69F00","#0072B2","#009E73","#CC79A7"))+
    scale_x_date(breaks="7 days") + 
    scale_y_continuous(expand=c(0,0),limits=c(0,0.02),breaks=seq(0,0.02,by=0.0025)) +
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
  
  tmp_hubei2 <- quants_summary_hubei %>% filter(date >= "2019-12-08" & date <= "2020-01-01" & var %in% names(var_key))
  tmp_hubei2$var <- var_key[tmp_hubei2$var]
  hubei_plot_inset <- tmp_hubei2 %>%
    ggplot() + 
    geom_vline(xintercept=(as.Date("2020-01-23", format="%Y-%m-%d", tz="LMT", origin="2019-11-01")),linetype="dashed") +
    geom_line(aes(x=date,y=median/pop_wuhan,col=var)) + 
    geom_ribbon(aes(x=date,ymin=lower/pop_wuhan,ymax=upper/pop_wuhan,fill=var),alpha=0.25) +
    scale_fill_manual(values=c("#E69F00","#0072B2","#009E73","#CC79A7"))+
    scale_color_manual(values=c("#E69F00","#0072B2","#009E73","#CC79A7"))+
    scale_x_date(breaks="2 days") + 
    scale_y_continuous(expand=c(0,0),limits=c(0,0.0002),breaks=seq(0,0.0002,by=0.000025)) +
    theme_bw() +
    ylab("") +
    xlab("") +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=6),
          panel.background = element_blank(),
          axis.text.y=element_text(size=6),
          legend.position="none",
          panel.grid.minor=element_blank())
  vp <- viewport(width=0.4,height=0.5, x=0.28,y=0.73)
  
  quants_summary_hubei1 <- quants_summary_hubei %>% filter(var %in% names(var_key))
  quants_summary_hubei1$var <- var_key[quants_summary_hubei1$var]
  quants_summary_hubei1$runname <- runname
  
  write_csv(quants_summary_hubei1, paste0(savewd,"/",filename,"_hubei_prevalence.csv"))
  
  pdf(paste0(savewd,"/",filename,"_hubei_prevalence.pdf"),height=6,width=8)
  print(hubei_plot)
  print(hubei_plot_inset, vp=vp)
  dev.off()
  png(paste0(savewd,"/",filename,"_hubei_prevalence.png"),height=6,width=8,res=300,units="in")
  print(hubei_plot)
  print(hubei_plot_inset, vp=vp)
  dev.off()
  
  ## Incidence
  quants1 <- quants_summary
  
  confirmed_data2 <- confirmed_data1
  confirmed_data2$date <- as.Date(confirmed_data2$date, origin="2019-11-01")
  confirmed_data2$province <- as.numeric(confirmed_data2$province)
  confirmed_data2$province <- provinces[confirmed_data2$province]
  confirmed_data2$province <- factor(confirmed_data2$province,levels=factor_levels_names)
  confirmed_data2 <- confirmed_data2 %>% drop_na()
  quants1$province <- factor(quants1$province,levels=factor_levels_names)
  quants1 <- quants1 %>% drop_na()
  
  var_key2 <- c("observations"="Simulated incidence of case confirmations",
                "onsets"="Symptom onset incidence",
                "infections"="Infection incidence")
  quants1 <- quants1 %>% filter(var %in% names(var_key2))
  quants1$var <- var_key2[quants1$var]
  quants1$var <- factor(quants1$var, levels=var_key2)
  
  p <- ggplot(quants1[quants1$date >= "2020-01-01",]) +
    geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
    geom_line(aes(x=date,y=median,col=var)) +
    geom_point(data=confirmed_data2[confirmed_data2$date >= "2020-01-01",],aes(x=date,y=n),size=0.5) +
    scale_x_date(breaks="7 days") +
    scale_fill_manual(values=cols[c(2,3,1)])+
    scale_color_manual(values=cols[c(2,3,1)])+
    geom_vline(xintercept=as.Date("2020-01-23",origin="2019-11-01"),linetype="dashed") +
    theme_pubr() +
    xlab("Date") +
    ylab("Daily incidence (absolute numbers)") +
    export_theme +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
          axis.text.y=element_text(size=7),
          axis.title = element_text(size=8),
          legend.text = element_text(size=7),
          strip.text = element_text(size=8),
          legend.title = element_blank(),
          legend.direction = "horizontal",
          legend.position=c(0.6,0),
          panel.grid.minor=element_blank()) +
    facet_wrap(~province,ncol=4,scales="free_y") 
  p
  write_csv(quants1,paste0(savewd,"/",filename,"_all_incidence.csv"))
  save(confirmed_data2, file=paste0(savewd,"/",filename,"_confirmed_data.RData"))
  
  pdf(paste0(savewd,"/",filename,"_all_incidence.pdf"),height=8,width=8)
  print(p)
  dev.off()
  
  png(paste0(savewd,"/",filename,"_all_incidence.png"),height=8,width=8,res=300,units="in")
  print(p)
  dev.off()
  quants <- generate_prediction_intervals(chains2, parTab, confirmed_data1, 
                                                  daily_import_probs = import_probs, daily_export_probs = export_probs_use,
                                                  time_varying_confirm_delay_pars = time_varying_report_pars,
                                                  nsamp=nsamp,return_draws = TRUE,model_ver="logistic",noise_ver="poisson",
                                                  incubation_ver="lnorm",scale_reporting=scale_reporting, 
                                          report_rate_switch=83,report_rate_1=0.14,report_rate_2=0.65)

  quants <- quants$draws
  quants$province <- provinces[quants$province]
  quants$date <- as.Date(quants$date, origin="2019-11-01")
  
  
  ## Get prevalence on 23rd and 30th
  tmp_dat <- quants %>% filter(province == "Hubei" & var == "cantravel_prevalence" & date == "2020-01-23") %>% 
    summarise(median=median(n/n_wuhan), lower=quantile(n/n_wuhan, c(0.025)), upper=quantile(n/n_wuhan, c(0.975)))
  tmp_dat$runname <- runname
  write_csv(tmp_dat, paste0(savewd,"/", filename,"_23rd_prev.csv"))
  metrics_prevalence_23rdJan <- bind_rows(metrics_prevalence_23rdJan, tmp_dat)
  
  ## Get prevalence on 23rd and 30th
  tmp_dat <- quants %>% filter(province == "Hubei" & var == "cantravel_prevalence" & date == "2020-01-30") %>% 
    summarise(median=median(n/n_wuhan), lower=quantile(n/n_wuhan, c(0.025)), upper=quantile(n/n_wuhan, c(0.975)))
  tmp_dat$runname <- runname
  write_csv(tmp_dat, paste0(savewd,"/", filename,"_30th_prev.csv"))
  metrics_prevalence_30thJan <- bind_rows(metrics_prevalence_30thJan, tmp_dat)
  
  tmp_dat <- quants %>% ungroup() %>% filter(var == "onsets") %>% group_by(province) %>% filter(n == max(n)) %>%
    summarise(median=median(date))
  tmp_dat$runname <- runname
  write_csv(tmp_dat, paste0(savewd,"/", filename,"_peak_times.csv"))

  metrics_peak_times <- bind_rows(metrics_peak_times, tmp_dat)
    
  ## 549 cases confirmed by 23rd Jan
  number_conf <- quants %>% filter(province == "Hubei" & date <= "2020-01-23" & var == "observations") %>%
    group_by(sampno) %>% 
    summarise(x=sum(n)) %>% pull(x)
  
  tmp_dat <- quantile(549/number_conf, c(0.025,0.5,0.975))
  tmp_dat <- data.frame(runname=runname, lower=tmp_dat[1],median=tmp_dat[2],upper=tmp_dat[3], metric="23rd Jan underreporting")
  write_csv(tmp_dat, paste0(savewd,"/", filename,"_23rd_underreporting.csv"))
  metrics_23rdJan <- bind_rows(metrics_23rdJan, tmp_dat)
  
  ## Or 30th, 5806
  number_conf2 <- quants %>% filter(province == "Hubei" & date <= "2020-01-30" & var == "observations") %>%
    group_by(sampno) %>% 
    summarise(x=sum(n)) %>% pull(x)
  tmp_dat <- quantile(5806/number_conf2, c(0.025,0.5,0.975))
  tmp_dat <- data.frame(runname=runname, lower=tmp_dat[1],median=tmp_dat[2],upper=tmp_dat[3], metric="30th Jan underreporting")
  write_csv(tmp_dat, paste0(savewd,"/", filename,"_30th_underreporting.csv"))
  metrics_30thJan <- bind_rows(metrics_30thJan, tmp_dat)
  
  quants_summary %>% filter(province == "Hubei" & 
                                var %in% c("presymptomatic_prevalence","symptomatic_prevalence","disease_prevalence") & 
                                date %in% c(as.Date("2020-01-23"),as.Date("2020-01-30"))) %>% 
    select(date, var, lower, median, upper) %>% ungroup() %>% 
    mutate(lower=lower/n_wuhan, median=median/n_wuhan,upper=upper/n_wuhan)
  
  ## Final cumulative incidence
  tmp_dat <- quantile(exp(chains$K)/n_wuhan,c(0.025,0.5,0.975))
  tmp_dat <- data.frame(runname=runname, lower=tmp_dat[1], median=tmp_dat[2], upper=tmp_dat[3], metric="cumu_incidence")
  write_csv(tmp_dat, paste0(savewd,"/", filename,"_final_cumu_incidence.csv"))
  
  import_probs <- read.csv("data/import_probs_matched.csv")
  provinces <- as.character(import_probs[,1])
  provinces[which(provinces == "Inner Mongolia")] <- "Inner_Mongolia"
  
  tmp <- chains2[,colnames(chains2) %like% "local_r"]
  tmp <- tmp[,!(colnames(tmp) %in% c("local_r_mean","local_r_sd"))]
  tmp <- tmp[,2:ncol(tmp)]
  tmp <- tmp/rowMeans(tmp)
  colnames(tmp) <- unique(confirmed_data2$province)
  tmp1 <- reshape2::melt(tmp)
  import_probs_melt <- reshape2::melt(import_probs)
  import_probs_melt$province_use <- as.character(import_probs_melt$province_use)
  import_probs_melt <- import_probs_melt %>% mutate(province_use = ifelse(province_use=="Inner Mongolia", "Inner_Mongolia",province_use))
  import_probs_melt <- import_probs_melt %>% group_by(province_use) %>% summarise(val=mean(value)) %>% filter(province_use != "Hubei") %>% ungroup()
  colnames(import_probs_melt) <- c("variable", "Mean import probability")
  import_probs_melt$variable <- factor(import_probs_melt$variable,levels=levels(tmp1$variable))
  tmp1 <- tmp1 %>% left_join(import_probs_melt)
  tmp1$variable <- factor(tmp1$variable, levels=levels(confirmed_data2$province))
  local_r_plot <- ggplot(tmp1) + 
    geom_violin(aes(x=variable, y=value, fill=`Mean import probability`),
                draw_quantiles=c(0.025,0.5,0.975),scale="width") + 
    scale_fill_gradient2(low="blue",mid="red",high="orange",midpoint=0.075) +
    geom_hline(yintercept=1,linetype="dashed") +
    guides(fill=guide_colourbar(title.position="top",
                                  direction="horizontal")) +
    theme_bw() + 
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.text = element_text(angle=45,hjust=1,size=7),
          legend.title=element_text(size=7,hjust=0.5),
          legend.position=c(0.8,0.83)) + 
    xlab("Province") + 
    labs(y=expression("Average number of local cases\n per import, R" [local]))+
    scale_y_continuous(expand=c(0,0),breaks=seq(0,6,by=1),limits=c(0,6))
  
  write_csv(tmp1, paste0(savewd,"/",filename,"_local_r_relative.csv"))
  
  pdf(paste0(savewd,"/",filename,"_local_r.pdf"),height=4,width=7)
  local_r_plot
  dev.off()
  
  png(paste0(savewd,"/",filename,"_local_r.png"),height=4,width=7,res=300,units="in")
  local_r_plot
  dev.off()

}

runname_key <- c(
  "main"="Analysis 1 (baseline, intermediate Wuhan prevalence)",
  "diffuse_r_local"="Analysis 2 (increased local transmission, lower Wuhan prevalence)", 
  "t_switch_6"="Analysis 3 (symptom peak 29th Jan, higher Wuhan prevalence)",
  "earliest_seed"="Analysis 4 (17th Nov seed)",
  "early_seed"="Analysis 5 (1st Dec seed)",
  "fewer_travellers"="Analysis 6 (4 million travellers)",
  "t_switch_1"="Analysis 7 (symptom peak 23rd Jan)",
  "t_switch_2"="Analysis 8 (symptom peak 24th Jan)",
  "t_switch_3"="Analysis 9 (symptom peak 26th Jan)",
  "t_switch_4"="Analysis 10 (symptom peak 27th Jan)",
  "t_switch_5"="Analysis 11 (symptom peak 28th Jan)",
  "fixed_serial_int" = "Analysis 12 (fixed serial interval distribution)"
)
metrics_23rdJan_backup <- metrics_23rdJan
metrics_23rdJan <- metrics_23rdJan %>% filter(runname != "r_local_unconstrained")
metrics_23rdJan$runname <- runname_key[metrics_23rdJan$runname]
metrics_23rdJan$estimate <- paste0(signif(metrics_23rdJan$median,3)," (", signif(metrics_23rdJan$lower,3), "-", signif(metrics_23rdJan$upper,3),")")
metrics_23rdJan1 <- metrics_23rdJan[,c("metric","runname","estimate")]
colnames(metrics_23rdJan1) <- c("Parameter","Scenario","Estimate")

metrics_30thJan_backup <- metrics_30thJan
metrics_30thJan <- metrics_30thJan %>% filter(runname != "r_local_unconstrained")
metrics_30thJan$runname <- runname_key[metrics_30thJan$runname]
metrics_30thJan$estimate <- paste0(signif(metrics_30thJan$median,3)," (", signif(metrics_30thJan$lower,3), "-", signif(metrics_30thJan$upper,3),")")
metrics_30thJan1 <- metrics_30thJan[,c("metric","runname","estimate")]
colnames(metrics_30thJan1) <- c("Parameter","Scenario","Estimate")

metrics_prevalence_23rdJan_backup <- metrics_prevalence_23rdJan
metrics_prevalence_23rdJan <- metrics_prevalence_23rdJan %>% filter(runname != "r_local_unconstrained")
metrics_prevalence_23rdJan$runname <- runname_key[metrics_prevalence_23rdJan$runname]
metrics_prevalence_23rdJan$estimate <- paste0(signif(metrics_prevalence_23rdJan$median,3)," (", signif(metrics_prevalence_23rdJan$lower,3), "-", signif(metrics_prevalence_23rdJan$upper,3),")")
metrics_prevalence_23rdJan$parameter <- "23rd January prevalence"
metrics_prevalence_23rdJan1 <- metrics_prevalence_23rdJan[,c("parameter","runname","estimate")]
colnames(metrics_prevalence_23rdJan1) <- c("Parameter","Scenario","Estimate")

metrics_prevalence_30thJan_backup <- metrics_prevalence_30thJan
metrics_prevalence_30thJan <- metrics_prevalence_30thJan %>% filter(runname != "r_local_unconstrained")
metrics_prevalence_30thJan$runname <- runname_key[metrics_prevalence_30thJan$runname]
metrics_prevalence_30thJan$estimate <- paste0(signif(metrics_prevalence_30thJan$median,3)," (", signif(metrics_prevalence_30thJan$lower,3), "-", signif(metrics_prevalence_30thJan$upper,3),")")
metrics_prevalence_30thJan$parameter <- "30th January prevalence"
metrics_prevalence_30thJan1 <- metrics_prevalence_30thJan[,c("parameter","runname","estimate")]
colnames(metrics_prevalence_30thJan1) <- c("Parameter","Scenario","Estimate")

#metrics_peak_times_backup <- metrics_peak_times
#metrics_peak_times$runname <- runname_key[metrics_peak_times$runname]

parameter_estimates_all_backup <- parameter_estimates_all
parameter_estimates_all$Scenario <- runname_key[parameter_estimates_all$Scenario]
colnames(parameter_estimates_all) <- c("Parameter","Estimate","Scenario")
parameter_estimates_all <- parameter_estimates_all[,c("Parameter","Scenario","Estimate")]

all_final <- bind_rows(metrics_23rdJan1,metrics_30thJan1,metrics_prevalence_23rdJan1,metrics_prevalence_30thJan1,parameter_estimates_all)
all_final <- all_final %>% arrange(Scenario, Parameter)

metrics_23rdJan1 %>% drop_na() %>% write_csv(path=paste0(savewd,"/23rd_underreporting.csv"))
metrics_30thJan1 %>% drop_na() %>% write_csv(path=paste0(savewd,"/30th_underreporting.csv"))
metrics_prevalence_23rdJan1 %>% drop_na() %>% write_csv(path=paste0(savewd,"/23rd_prevalence.csv"))
metrics_prevalence_30thJan1 %>% drop_na() %>% write_csv(path=paste0(savewd,"/30th_prevalence.csv"))
parameter_estimates_all %>% drop_na() %>% write_csv(path=paste0(savewd,"/par_estimates.csv"))

write_csv(all_final, path=paste0(savewd,"/all_estimates.csv"))


