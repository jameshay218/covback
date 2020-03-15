#' @export
create_model_func_provinces <- function(parTab, data=NULL, PRIOR_FUNC=NULL,
                                        tmax=NULL, confirm_delay_pars=NULL,
                                        daily_import_probs=NULL, daily_export_probs=NULL,
                                        ver="posterior", noise_ver="poisson",
                                        model_ver=1, solve_prior=FALSE){
    par_names <- parTab$names
    par_provinces <- parTab$province
    unique_provinces <- unique(parTab$province)
    
    ## Start from last province
    unique_provinces <- unique_provinces[unique_provinces != "all"]
    n_provinces <- length(unique_provinces)
    
    if (!is.null(data)) {
        data <- data %>% arrange(province, date)
        times <- data %>% filter(province == "1") %>% pull(date)
        cases <- data$n
        bigT <- length(times)
        tmax <- max(times)
        data <- data %>% arrange(province, date)
      } else {
        if (is.null(tmax)) {
            stop("One of data and tmax must be specified")
        }
        if (ver == "posterior") {
            stop("Data must be specified to calculate the posterior")
        }
        data <- expand.grid(date=0:tmax,province=unique_provinces)
        data <- as_tibble(data)
        data <- data %>% arrange(province, date)
        
        bigT <- tmax + 1
    }

    ## Exportation and importation probabilities
    leave_matrix <- prob_leave_on_day(daily_export_probs, tmax)
    ## For each province, what's the daily probability of receiving a person from the seed province?
    arrival_matrices <- NULL
    ## First province isn't meaningful
    for (i in 1:n_provinces){
        arrival_matrices[[i]] <- prob_daily_arrival(daily_export_probs, daily_import_probs[i,], tmax)
    }

    if(!is.null(confirm_delay_pars)){
        report_delay_mat <- calculate_reporting_delay_matrix(confirm_delay_pars$shape, confirm_delay_pars$scale)
    }
    
     model_func <- function(pars_all) {
        names(pars_all) <- par_names
        if(!solve_prior){
          pars_seed <- pars_all[which(par_provinces == "1")]
          
          ## Gamma distribution
          shape <- pars_all["shape"]
          scale <- pars_all["scale"]
  
          ## If time-varying parameters not specified, enumerate out the point
          ## estimates
          ## Move into model func call if need to estimate these...
          if (is.null(confirm_delay_pars)) {
              report_delay_mat <- calculate_reporting_delay_matrix_constant(shape,scale,tmax)
          }
          
          ## Negative binomial size
          size <- pars_all["size"]
  
          ## Incubation period    
          weibull_alpha <- pars_all["weibull_alpha"]
          weibull_sigma <- pars_all["weibull_sigma"]
          #print(weibull_alpha)
          #print(weibull_sigma)
          ## Serial interval
          #serial_lmean <- pars_all["lnorm_mean"]
          #serial_lsd <- pars_all["lnorm_sd"]
          serial_interval_par1 <- pars_all["serial_interval_gamma_par1"]
          serial_interval_par2 <- 1/pars_all["serial_interval_gamma_par2"]
  
          serial_probs <- calculate_serial_interval_probs(tmax, serial_interval_par1, serial_interval_par2)
          
          ## For each day with a potential infection onset, get the probability of leaving at some point in the future before
          ## symptom onset
          ## daily_prob_leaving <- prob_left_pre_sympt(pars_seed["export_prob"], weibull_alpha, weibull_sigma, 100)
          presymptom_probs <- calculate_probs_presymptomatic(tmax, weibull_alpha, weibull_sigma)
          daily_prob_leaving <- prob_leave_pre_symptoms_vector(leave_matrix, presymptom_probs)*pars_all["export_prob"]
          
          onset_probs <- calculate_onset_probs(tmax, weibull_alpha, weibull_sigma)
          
  
          ## Get local growth of seed province. Will add/subtract this from later estimates
          if(model_ver == 1){
            infections_seed <- pars_seed["i0"]*daily_exp_interval_cpp(pars_seed["growth_rate"], tmax, pars_seed["t0"])
          } else {
            t_switch <- pars_all["t_switch"]
            if(weibull_alpha > 1) {
              inc_mode <- weibull_mode(weibull_alpha, weibull_sigma)
            } else {
              inc_mode <- 0
            }
            #print(inc_mode)
            t_switch <- t_switch - inc_mode
            growth_rate <- log(pars_all["K"]-1)/t_switch
            infections_seed <- daily_sigmoid_interval_cpp(growth_rate, pars_all["K"], tmax, pars_seed["t0"])
          }
          
          all_infections <- numeric(bigT*n_provinces)
          all_confirmations <- numeric(bigT*n_provinces)
          all_onsets <- numeric(bigT*n_provinces)
          
          all_dat <- NULL
          ## Need to loop through each province
          index <- 1
          for(province in unique_provinces) {
              pars <- pars_all[which(par_provinces %in% c("all",province))]
              
              ## Proportion of seed cases that will get exported
              ## If province 1, then export_propn cases are lost
              ## Otherwise, some proportion of these lost cases are gained elsewhere
              if(province == "1") {
                  ## Cases that leave seed province
                  import_cases <- -1 * daily_prob_leaving * infections_seed
                  t0 <- pars_seed["t0"]
              } else {
                  ## Otherwise, cases that come from seed province
                  daily_prob_arrival <- prob_arrive_pre_symptoms_vector(arrival_matrices[[index]], presymptom_probs)*pars_all["export_prob"]
                  import_cases <- daily_prob_arrival * infections_seed
  
                  ## t0 is days since t0 in seed province
                  t0_import <- pars_seed["t0"]
                  t0 <- pars["t0"] + t0_import
                  t0 <- min(t0, tmax-1)
              }
              
              ## Growth model
              growth_rate <- pars_seed["growth_rate"]
              i0 <- pars["i0"]
              
              ## Solve model for this province
              import_cases_local <- calculate_local_from_import_infections(import_cases*pars["local_r"], serial_probs, tmax)
              import_cases_local[is.na(import_cases_local)] <- 0
              import_cases <- import_cases + import_cases_local
              
              if(model_ver == 1) {
                res <- calculate_all_incidences(growth_rate, t0, i0, import_cases,
                                                onset_probs, report_delay_mat,
                                                tmax)
              } else {
                t_switch <- pars_all["t_switch"]
                inc_mode <- weibull_mode(weibull_alpha, weibull_sigma)
                t_switch <- t_switch - inc_mode
                growth_rate <- log(pars_all["K"]-1)/t_switch
                res <- calculate_all_incidences_logistic(growth_rate, t0, i0, pars_all["K"],
                                                         import_cases,
                                                onset_probs, report_delay_mat,
                                                tmax)
              }
              
              
              ## Extract the 3 incidence types
              infections <- res$infections
              onsets <- res$onsets
              confirmations <- res$confirmations
              
              indices <- (bigT*(index - 1) + 1):(bigT*index)
              index <- index + 1
              all_infections[indices] <- infections
              all_onsets[indices] <- onsets
              all_confirmations[indices] <- confirmations
          }
           ## Once done, combined all of the provinces into one tibble
          all_dat <- data %>% bind_cols(infections=all_infections, onsets=all_onsets, confirmations=all_confirmations)
                          
          ## If version to solve is model, we're done. Otherwise solve likelihood for confirmations
          if(ver == "model"){
            all_dat <- all_dat %>% pivot_longer(-c("date","province"),names_to="var",values_to="n")
              return(all_dat)
          } else {
              if (noise_ver == "poisson") {
                  lik <- sum(dpois(x=cases, all_confirmations, log=TRUE),na.rm=TRUE)
              } else {
                  lik <- sum(dnbinom(x=cases,mu=all_confirmations,size=size,log=TRUE),na.rm=TRUE)
              }
            
              if (!is.null(PRIOR_FUNC)) {
                  lik <- lik + PRIOR_FUNC(pars_all)
              }
              #print(lik)
              return(lik)
          }
        } else {
          lik <- -100000
          if (!is.null(PRIOR_FUNC)) {
            lik <- lik + PRIOR_FUNC(pars_all)
          }
          return(lik)
        }
    }
    return(model_func)
}
