#' @export
plot_simulations <- function(sim_dat){
    p <- ggplot(sim_dat) +
        geom_line(aes(x=date,y=n,col=var)) +
        theme_bw() +
        theme(legend.position=c(0.2,0.8))
    p
}

#' @export
plot_model_fit <- function(chain, parTab, data, confirm_delay_pars=NULL,
                           daily_import_probs, daily_export_probs,
                           nsamp=1000,
                           add_noise=TRUE, noise_ver="poisson"){
    imports_stop <- parTab[parTab$names == "import_stop","values"]
    dat_plot <- data %>% filter(var %in% c("date_report_observable", "date_infection_true"))
    dat1 <- data %>% filter(var == "date_report_observable") %>% select(date, n, province)
    quants <- generate_prediction_intervals(chain, parTab, dat1, confirm_delay_pars=confirm_delay_pars,
                                            daily_import_probs = daily_import_probs, daily_export_probs = daily_export_probs,
                                            nsamp, add_noise, noise_ver)
    quants1 <- quants %>% filter(var %in% c("infections","observations",
                                 "date_onset_true","date_infection_true"))
    p <- ggplot(quants1) +
        geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
        geom_line(aes(x=date,y=median,col=var)) +
        geom_point(data=dat_plot,aes(x=date,y=n,col=var),size=0.5) +
        geom_ribbon(data=dat_plot, aes(x=date,ymin=0,ymax=0,fill=var)) +
        geom_vline(xintercept=imports_stop,linetype="dashed") +
        theme_bw() + theme(legend.position=c(0.8, 0.05)) +
        facet_wrap(~province,scales="free_y",ncol=5)
    p
}

#' @export
generate_prediction_intervals <- function(chain, parTab, data, time_varying_confirm_delay_pars=NULL,
                                          daily_import_probs, daily_export_probs,
                                          nsamp=1000,
                                          add_noise=TRUE, noise_ver="poisson",
                                          return_draws=FALSE,
                                          model_ver=1){
    model_func <- create_model_func_provinces(parTab, data, time_varying_confirm_delay_pars=time_varying_confirm_delay_pars,
                                              daily_import_probs = daily_import_probs, daily_export_probs = daily_export_probs,
                                              ver="model",model_ver=model_ver)
    par_names <- parTab$names

    samps <- sample(unique(chain$sampno), nsamp)
    
    store_all <- NULL
    for(i in seq_along(samps)){
        pars <- get_index_par(chain, samps[i])
        names(pars) <- par_names
        
        prob_presymptomatic  <- calculate_probs_presymptomatic(100, pars["weibull_alpha"], pars["weibull_sigma"])
        prob_preconfirmation  <- calculate_probs_preconfirmation(100, pars["shape"], pars["scale"])
        
        res <- model_func(pars)

        infection_prevalence <- res %>% 
          filter(var == "infections") %>% 
          group_by(province) %>%
          mutate(n = calculate_infection_prevalence(n, prob_presymptomatic),
                 var = "infection_prev")
        
        symptomatic_prevalence <- res %>% 
          filter(var == "infections") %>% 
          group_by(province) %>%
          mutate(n = calculate_infection_prevalence(n, 1-prob_presymptomatic),
                 var = "symptomatic_prevalence")
        
        onset_prevalence <- res %>% 
          filter(var == "onsets") %>% 
          group_by(province) %>%
          mutate(n = calculate_infection_prevalence(n, prob_preconfirmation),
                 var = "onset_prev")
        
        total_prev <- bind_rows(onset_prevalence, infection_prevalence) %>% 
          group_by(province, date) %>%
          summarise(n= sum(n)) %>%
          mutate(var="total_prevalence") %>% ungroup()

        #confirmations <- res %>% filter(var == "confirmations") %>% pull(n)
        if (add_noise) {
            subset_confirmations <- res %>% filter(var == "confirmations")
            if (noise_ver == "poisson") {
                subset_confirmations <- subset_confirmations %>% 
                    mutate(n = rpois(n(),n)) %>% 
                    mutate(var = "observations")
            } else if (noise_ver == "nbinom") {
                subset_confirmations <- subset_confirmations %>% 
                    mutate(n = rnbinom(n(),mu=n,size=pars["size"])) %>% 
                    mutate(var = "observations")
            } else {
                subset_confirmations <- subset_confirmations %>% mutate(var = "observations")
            }
        }
        res <- bind_rows(res, subset_confirmations, infection_prevalence, onset_prevalence, total_prev, symptomatic_prevalence)
        res$sampno <- i
        store_all <- bind_rows(store_all, res)
    }
    if(return_draws) {
      return(list(draws=store_all,samp_ids=samps))
    }
    quants <- store_all %>% group_by(province, date, var) %>%
        do(data.frame(t(c(quantile(.$n, probs = c(0.01,0.025,0.25,0.5,0.75,0.975,0.99),na.rm=TRUE),mean(.$n)))))
    colnames(quants) <- c("province","date","var",
                                   "min","lower","midlow","median","midhigh","upper","max","mean")
  
    return(quants)    
}

#' @export
plot_posteriors <- function(chain, parTab){
    parTab_estimates <- parTab[parTab$fixed == 0,]
    par_names <- parTab_estimates$names
    chain1 <- chain[,c("sampno",par_names)]

    real_values <- parTab_estimates[,c("names","values")]
    colnames(real_values) <- c("variable","values")

    melt_chain <- chain1 %>% pivot_longer(-sampno, names_to="variable",values_to="value")
    
    estimate_p <- ggplot(melt_chain) + geom_density(aes(x=value)) +
        geom_vline(data=real_values,aes(xintercept=values),linetype="dashed") +
        theme_bw() +
        ggtitle("Dashed line shows true simulated value") +
        facet_wrap(~variable,scales="free")
    estimate_p
    
}
#' @export
plot_reporting_landscape <- function(shapes, scales){
  tmp <- calculate_reporting_delay_matrix(shapes, scales)
  tmp <- as.data.frame(tmp)
  tmp$time <- 0:(nrow(tmp)-1)
  tmp <- reshape2::melt(tmp,id.vars="time")
  tmp$variable <- as.numeric(as.factor(tmp$variable))
  tmp$variable <- tmp$time - tmp$variable
  colnames(tmp) <- c("Time of report","Delay from onset", "Probability")
  tmax <- max(tmp$`Time of report`)
  tmin <- min(tmp$`Time of report`)
  ggplot(tmp) + geom_raster(aes(x=`Time of report`,y=`Delay from onset`, fill=`Probability`)) + 
    scale_fill_gradient(low="blue",high="red") + coord_cartesian(xlim=c(tmin,tmax),ylim=c(tmin,tmax)) + 
    scale_y_continuous(expand=c(0,0)) + 
    scale_x_continuous(expand=c(0,0)) +
    theme_bw() + 
    theme(legend.position=c(0.2,0.7))
}

#' @export
plot_confirm_delay <- function(chain, nsamp=100,xmax=40){
  samps <- sample(unique(chain$sampno), nsamp)
  res <- matrix(nrow=nsamp, ncol=xmax+1)
  
  for(i in seq_along(samps)){
    pars <- get_index_par(chain, samps[i])
    res[i,] <- dgamma(0:xmax, shape=pars["shape"],scale=pars["scale"])
  }
  quants <- t(apply(res, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
  quants <- data.frame(quants)
  colnames(quants) <- c("lower","median","upper")
  quants$day <- 0:xmax
  p <- ggplot(quants) + geom_ribbon(aes(x=day,ymin=lower,ymax=upper),alpha=0.25) + geom_line(aes(x=day,y=median)) +
    theme_bw()
  return(p)
}

#' @export
plot_incubation_period <- function(chain, nsamp=100,xmax=40,prior_pars=NULL){
  samps <- sample(unique(chain$sampno), nsamp)
  res <- matrix(nrow=nsamp, ncol=xmax+1)
  prior_res <- matrix(nrow=nsamp, ncol=xmax+1)
  
  if(!is.null(prior_pars)){
    samp_prior <- sample(1:nrow(prior_pars), nsamp)
  }
  
  for(i in seq_along(samps)){
    pars <- get_index_par(chain, samps[i])
    res[i,] <- dweibull(0:xmax, shape=pars["weibull_alpha"],scale=pars["weibull_sigma"])
    
    if(!is.null(prior_pars)){
    prior_alpha <- prior_pars$alpha[samp_prior[i]]
    prior_sigma <- prior_pars$sigma[samp_prior[i]]
    
    prior_res[i,] <- dweibull(0:xmax, shape=prior_alpha, scale=prior_sigma)
    }
    }
  quants <- t(apply(res, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
  quants <- data.frame(quants)
  colnames(quants) <- c("lower","median","upper")
  quants$day <- 0:xmax
  quants$var <- "Posterior"
  if(!is.null(prior_pars)){
    quants_prior <- t(apply(prior_res, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
    quants_prior <- data.frame(quants_prior)
    colnames(quants_prior) <- c("lower","median","upper")
    quants_prior$day <- 0:xmax
    quants_prior$var <- "Prior"
    quants <- bind_rows(quants, quants_prior)
  }
  
  p <- ggplot(quants) + 
    geom_ribbon(aes(x=day,ymin=lower,ymax=upper, fill=var),alpha=0.25) + 
    geom_line(aes(x=day,y=median,col=var)) +
    theme_bw() +
    theme(legend.position=c(0.8,0.8))
  return(p)
}

