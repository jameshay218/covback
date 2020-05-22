#' @export
create_model_func_provinces_fixed <- function(parTab, 
                                              data=NULL, 
                                              PRIOR_FUNC=NULL,
                                              tmax=NULL, 
                                              time_varying_confirm_delay_pars=NULL,
                                              daily_import_probs=NULL, 
                                              daily_export_probs=NULL,
                                              ver="posterior", 
                                              noise_ver="poisson",
                                              incubation_ver="lnorm",
                                              model_ver="logistic", 
                                              solve_prior=FALSE, 
                                              calculate_prevalence=FALSE){
  par_names <- parTab$names
  par_provinces <- parTab$province
  unique_provinces <- unique(parTab$province)
  pars_all <- parTab$values
  names(pars_all) <- par_names
  ## Start from last province
  unique_provinces <- unique_provinces[unique_provinces != "all"]
  n_provinces <- length(unique_provinces)
  ##############################################################
  ## i) SETUP DATA
  ##############################################################
  ## If not provided data to fit to, then generate skeleton for solving model
  if (!is.null(data)) {
    data <- data %>% arrange(province, date)
    times <- data %>% filter(province == "1") %>% pull(date)
    times_date <- as.Date(times, origin="2019-11-01")
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
  ##############################################################
  ## ii) SETUP TRAVEL PROBABILITIES
  ##############################################################
  ## Exportation and importation probabilities, if >1 province
  if(n_provinces > 1) {
    leave_matrix <- prob_leave_on_day(daily_export_probs, tmax)
    ## For each province, what's the daily probability of receiving a person from the seed province?
    arrival_matrices <- NULL
    ## First province isn't meaningful
    for (i in 1:n_provinces){
      arrival_matrices[[i]] <- prob_daily_arrival(daily_export_probs, daily_import_probs[i,], tmax)
    }
  }
  ##############################################################
  ## iii) SETUP REPORTING DELAYS
  ##############################################################
  ## If a time-varying confirmation delay distribution is not provided, generate a fixed one
  if(!is.null(time_varying_confirm_delay_pars)){
    report_delay_mat <- calculate_reporting_delay_matrix(time_varying_confirm_delay_pars$shape, time_varying_confirm_delay_pars$scale)
  }
  
  ## Check if will need to recalculate confirmation delay each iteration
  recalc_confirm_delay <- TRUE
  if(all(parTab[parTab$names %in% c("confirm_delay_shape","confirm_delay_scale"),"fixed"] == 1)){
    recalc_confirm_delay <- FALSE
    
    shape <- pars_all["confirm_delay_shape"]
    scale <- pars_all["confirm_delay_scale"]
        
    ## If time-varying parameters not specified, enumerate out the point estimates
    if (is.null(time_varying_confirm_delay_pars)) {
      report_delay_mat <- calculate_reporting_delay_matrix_constant(shape,scale,tmax)
    }
  }
  ##############################################################
  ## iv) SETUP INCUBATION PERIOD
  ##############################################################
  ## Check if need to recalculate incubation period each iteration
  recalc_incubation_period <- TRUE
  if((incubation_ver == "weibull" & all(parTab[parTab$names %in% c("weibull_alpha","weibull_sigma"),"fixed"] == 1)) |
     (incubation_ver == "lnorm" & all(parTab[parTab$names %in% c("lnorm_incu_par1","lnorm_incu_par2"),"fixed"] == 1))) {
    recalc_incubation_period <- FALSE
    if(incubation_ver == "weibull"){
      weibull_alpha <- pars_all["weibull_alpha"]
      weibull_sigma <- pars_all["weibull_sigma"]
      ## For each day with a potential infection onset, get the probability of leaving at some point in the future before
      ## symptom onset
      presymptom_probs <- calculate_probs_presymptomatic_weibull(tmax, weibull_alpha, weibull_sigma)
      onset_probs <- calculate_onset_probs_weibull(tmax, weibull_alpha, weibull_sigma)
    } else {
      incu_par1 <- pars_all["lnorm_incu_par1"]
      incu_par2 <- pars_all["lnorm_incu_par2"]
      ## For each day with a potential infection onset, get the probability of leaving at some point in the future before
      ## symptom onset
      presymptom_probs <- calculate_probs_presymptomatic_lnorm(tmax, incu_par1, incu_par2)
      onset_probs <- calculate_onset_probs_lnorm(tmax, incu_par1, incu_par2)
    }
    if(n_provinces > 1){
      ## This is used to get probability of arriving pre-symptomatic, used for calculation of local cases
      daily_prob_leaving_seed <- prob_leave_pre_symptoms_vector(leave_matrix, presymptom_probs)
      daily_prob_arrival_toa_precalc <- NULL
      daily_prob_arrival_toi_precalc <- NULL
      for (i in 1:n_provinces){
        prob_arrival_matrix <- prob_arrive_pre_symptoms(arrival_matrices[[i]],presymptom_probs)
        
        ## Get daily probability of someone arriving in this province
        daily_prob_arrival_toa_precalc[[i]] <- local_travel_matrix_precalc(prob_arrival_matrix)
        ## Daily probability of arriving at some point before becoming symptomatic
        daily_prob_arrival_toi_precalc[[i]] <- prob_arrive_pre_symptoms_vector(arrival_matrices[[i]], presymptom_probs)
      }
    }
  }
  ##############################################################
  ## v) SETUP SERIAL INTERVAL
  ##############################################################
  ## If not re-estimating serial interval parameters, calculate the serial interval here
  recalc_serial_interval <- TRUE
  if(all(parTab[parTab$names %in% c("serial_interval_mean","serial_interval_var"),"fixed"] == 1)){
    recalc_serial_interval <- FALSE
    ## Serial interval
    serial_pars <- gamma_pars_from_mean_sd(pars_all["serial_interval_mean"], pars_all["serial_interval_var"])
    
    serial_interval_alpha <- serial_pars[[1]]
    serial_interval_scale <- serial_pars[[2]]
    
    serial_probs <- calculate_serial_interval_probs(tmax, serial_interval_alpha, serial_interval_scale)
   
  }

  model_func <- function(pars_all) {
    names(pars_all) <- par_names
    browser()
    if(!solve_prior){
      ###########################################################
      ## STEP A - EXTRACT PARAMETERS
      ###########################################################
      pars_seed <- pars_all[which(par_provinces == "1")]
      
      ## Negative binomial size
      size <- pars_all["size"]
      
      ###########################################################
      ## i) get confirmation delay matrix, used for model fitting
      ###########################################################
      if(recalc_confirm_delay){
      ## Gamma distribution
        
        shape <- pars_all["confirm_delay_shape"]
        scale <- pars_all["confirm_delay_scale"]
        ## If time-varying parameters not specified, enumerate out the point
        ## estimates
        ## Move into model func call if need to estimate these...
        if (is.null(time_varying_confirm_delay_pars)) {
          report_delay_mat <- calculate_reporting_delay_matrix_constant(shape,scale,tmax)
        }
      }
      ###########################################################
      ## ii) get incubation period distribution
      ###########################################################
      if(recalc_incubation_period){
        ## Incubation period
        if(incubation_ver == "weibull"){
          weibull_alpha <- pars_all["weibull_alpha"]
          weibull_sigma <- pars_all["weibull_sigma"]
          ## For each day with a potential infection onset, get the probability of leaving at some point in the future before
          ## symptom onset
          presymptom_probs <- calculate_probs_presymptomatic_weibull(tmax, weibull_alpha, weibull_sigma)
          onset_probs <- calculate_onset_probs_weibull(tmax, weibull_alpha, weibull_sigma)
        } else {
          incu_par1 <- pars_all["lnorm_incu_par1"]
          incu_par2 <- pars_all["lnorm_incu_par2"]
          ## For each day with a potential infection onset, get the probability of leaving at some point in the future before
          ## symptom onset
          presymptom_probs <- calculate_probs_presymptomatic_lnorm(tmax, incu_par1, incu_par2)
          onset_probs <- calculate_onset_probs_lnorm(tmax, incu_par1, incu_par2)
        }
      }
      
      ###########################################################
      ## iii) get serial interval distribution to generate secondary cases over
      ###########################################################
      if(recalc_serial_interval){
        ## Serial interval
        serial_pars <- gamma_pars_from_mean_sd(pars_all["serial_interval_mean"], pars_all["serial_interval_var"])
        
        serial_interval_alpha <- serial_pars[[1]]
        serial_interval_scale <- serial_pars[[2]]
        
        serial_probs <- calculate_serial_interval_probs(tmax, serial_interval_alpha, serial_interval_scale)
      }
      ###########################################################
      ## iv) get probability of leaving pre-symptoms, given travel data
      ###########################################################
      ## If only one province, then you don't have anywhere to leave to
      if(n_provinces > 1) {
        if(recalc_serial_interval) {
          daily_prob_leaving <- prob_leave_pre_symptoms_vector(leave_matrix, presymptom_probs)
        } else {
          daily_prob_leaving <- daily_prob_leaving_seed
        }
      } else {
        daily_prob_leaving <- numeric(tmax + 1)
      }
      
      ###########################################################
      ## v) find shift to infection incidence inflection point 
      ##    needed to get onset inflection point on known date
      ###########################################################
      if(model_ver == "logistic" & 
         parTab[parTab$names == "growth_rate" & 
                parTab$province == "1","fixed"] == 1){
        ## Transform back to linear scale
        K <- exp(pars_all["K"])
        
        ## Peak symptom onsets t_switch_onsets days after seeding
        ## pars_all["t_switch"] is the actual date of peak symptom onset
        t_switch_onsets <- pars_all["t_switch"] - pars_seed["t0"]
        
        ## When do infections need to have peaked to give symptom onset peak
        ## on the specified day?
        find_tswitch_offset <- function(){
           f <- function(t_unknown){
             ## Infections peak t_unknown days before onsets
               t_switch_prime <- t_switch_onsets - t_unknown
               ## Growth rate in terms of K and time of peak infections
               growth <- log(K-1)/t_switch_prime
               ## Do the convolution with the current incubation period
               y <- daily_sigmoid_interval_cpp(growth, K, tmax, 0)
               onsets <- calculate_onset_incidence(y, onset_probs,tmax)
         
               ## What's the difference between the peak symptom onset this gives
               ## and desired peak symptom onset?
               t_switch_unknown <- which.max(onsets)
               (t_switch_unknown - t_switch_onsets)^2
        }
       
      ## Increase delay between peak infections and peak onsets until time of peak onsets
      ## matches pars_all["t_switch"]-pars_seed["t0"]
       t_unknown <- 0
       test1 <- f(t_unknown)
       t_unknown <- t_unknown + 1
       test2 <- f(t_unknown)
       while(test2 < test1){
         test1 <- test2
         t_unknown <- t_unknown + 1
         test2 <- f(t_unknown)
       }
       t_unknown <- t_unknown - 1
       return(t_unknown)
     }
      # ## T-switch is based on symptom onsets. So t-switch (peak) for infections is
      # ## peak of symptom onsets minus mode of incubation period distribution
       t_offset <- find_tswitch_offset()
       #t_offset <- 5
       #t_offset <- pars_all["incu_mode"]
      }
      ###########################################################
      ## STEP 1 - GROWTH OF INFECTIONS IN SEED PROVINCE
      ###########################################################          
      ## Get local growth of seed province. Will add/subtract this from later estimates
      ## Exponential growth for model 1, logistic for model 2
      if(model_ver == "exp"){
        infections_seed <- pars_seed["i0"]*daily_exp_interval_cpp(pars_seed["growth_rate"], tmax, pars_seed["t0"])
      } else {
        ## Can define logistic growth rate in terms of K and inflection point
        K <- exp(pars_all["K"])
        if(parTab[parTab$names == "growth_rate" & parTab$province == "1","fixed"] == 1){
          t_switch <- pars_all["t_switch"] - pars_seed["t0"] - t_offset
          growth_rate <- log(K-1)/t_switch
        } else {
          growth_rate <- pars_seed["growth_rate"]
          t_switch <- log(K-1)/growth_rate
        }
        infections_seed <- daily_sigmoid_interval_cpp(growth_rate, K, tmax, pars_seed["t0"])
      }
      
      ## Storage for numbers about to be solved
      all_infections <- numeric(bigT*n_provinces)
      all_confirmations <- numeric(bigT*n_provinces)
      all_onsets <- numeric(bigT*n_provinces)
      all_local_infections <- numeric(bigT*n_provinces)
      all_imported_infections <- numeric(bigT*n_provinces)
      
      all_presymptomatic_prevalence <- numeric(bigT*n_provinces)
      all_symptomatic_prevalence <- numeric(bigT*n_provinces)
      all_preconfirmation_prevalence <- numeric(bigT*n_provinces)
      all_cantravel_prevalence <- numeric(bigT*n_provinces)
      all_disease_prevalence <- numeric(bigT*n_provinces)
      
      all_dat <- NULL
      ## Need to loop through each province
      index <- 1
      for(province in unique_provinces) {
        
        ## Get parameters that apply to this province
        pars <- pars_all[which(par_provinces %in% c("all",province))]
        
        ## Proportion of seed cases that will get exported
        ## If province 1, then export_propn cases are lost
        ## Otherwise, some proportion of these lost cases are gained elsewhere
        if(province == "1") {
          ## Cases that leave seed province
          import_cases_time_of_infection <- -1 * daily_prob_leaving * infections_seed
          t0 <- pars_seed["t0"]
          import_cases_local <- numeric(length(import_cases_time_of_infection))
        } else {
          ## If we re-calcualted the incubation period, need to recalculate the probabilities of arriving
          ## presymptomatic. Otherwise, can use pre-computed values
          if(recalc_incubation_period){
            ## Otherwise, cases that come from seed province
            prob_arrival_matrix <- prob_arrive_pre_symptoms(arrival_matrices[[index]],presymptom_probs)
            
            daily_prob_arrival_toi <- prob_arrive_pre_symptoms_vector(arrival_matrices[[index]], presymptom_probs)
            daily_prob_arrival_toa <- local_travel_matrix_precalc(prob_arrival_matrix)
          } else {
            daily_prob_arrival_toi <- daily_prob_arrival_toi_precalc[[index]]
            daily_prob_arrival_toa <- daily_prob_arrival_toa_precalc[[index]]
          }
          
          ## This gives number of imported case **at the time that they got infected**
          ## For the backcalculation, we need time of infection. But for calculating
          ## local cases and prevalence we need time of arrival in destination.
          import_cases_time_of_infection <- daily_prob_arrival_toi * infections_seed
          
          ## t0 is days since t0 in seed province
          t0_import <- pars_seed["t0"]
          t0 <- pars["t0"] + t0_import
          t0 <- min(t0, tmax-1)
          if(is.na(t0)){
            t0 <- t0_import
          }
          
          ## Solve model for this province
          local_r <- pars["local_r"]
          ## This finds the number of locally generated cases. We find the time that infections in the seed location arrive,
          ## and then generate local_r cases over the remainder of their serial interval spent in the new province
          import_cases_local <- calculate_local_cases(daily_prob_arrival_toa, infections_seed, serial_probs, local_r)
          import_cases_local[is.na(import_cases_local)] <- 0
        }
        
        ## Growth model
        if(model_ver == "exp"){
          growth_rate <- pars_seed["growth_rate"]
        } else {
          ## Can define logistic growth rate in terms of K and inflection point
          K <- exp(pars_all["K"])
          if(parTab[parTab$names == "growth_rate" & parTab$province == "1","fixed"] == 1){
            t_switch <- pars_all["t_switch"] - pars_seed["t0"] - t_offset
            growth_rate <- log(K-1)/t_switch
          } else {
            growth_rate <- pars_seed["growth_rate"]
            t_switch <- log(K-1)/growth_rate
          }
        }
        
        i0 <- pars["i0"]
        if(is.na(i0)) {
          i0 <- 0
        }
        
        
        ## Total number of cases is number of locally generated cases plus imported cases
        import_cases <- import_cases_time_of_infection + import_cases_local
        ## Exponential growth model
        if(model_ver == "exp") {
          res <- calculate_all_incidences(growth_rate, t0, i0, import_cases,
                                          onset_probs, report_delay_mat,
                                          tmax)
          ## Logistic growth model
        } else {          
          res <- calculate_all_incidences_logistic(growth_rate, t0, i0, exp(pars_all["K"]),
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
        all_local_infections[indices] <- import_cases_local
        all_imported_infections[indices] <- import_cases_time_of_infection
        
        if(calculate_prevalence){
          if(n_provinces > 1){
            if(province == "1") {
              prob_not_left_matrix <- probs_not_left_by_day(daily_export_probs, tmax)
              tmp_presymptomatic_prevalence_local <- calculate_infection_prevalence_hubei(infections, presymptom_probs,
                                                                                          prob_not_left_matrix)
              tmp_presymptomatic_prevalence_imported <- numeric(tmax+1)
            } else {
              tmp_presymptomatic_prevalence_local <- calculate_infection_prevalence_local(import_cases_local, presymptom_probs)
              tmp_presymptomatic_prevalence_imported <- calculate_infection_prevalence_imported(infections_seed, 
                                                                                                presymptom_probs, daily_prob_arrival_toa)
            }
          }
          
          #confirm_pars <- gamma_pars_from_mean_sd(pars["confirm_delay_mean"], pars["confirm_delay_sd"]^2)
          #shape <- confirm_pars[[1]]
          #scale <- confirm_pars[[2]]
          
          shape <- pars_all["confirm_delay_shape"]
          scale <- pars_all["confirm_delay_scale"]
          
          prob_preconfirmation  <- calculate_probs_preconfirmation(tmax, shape, 
                                                                   scale)
          prob_prerecovery <- calculate_probs_notrecovered(tmax, pars["recovery_shape"], 
                                                           pars["recovery_scale"])
          
          if(n_provinces > 1){
            all_presymptomatic_prevalence[indices] <- tmp_presymptomatic_prevalence_local + tmp_presymptomatic_prevalence_imported
          } else {
            all_presymptomatic_prevalence[indices] <- 0
          }
          all_symptomatic_prevalence[indices] <- calculate_unrecovered_prevalence(onsets, prob_prerecovery)
          
          if(!is.null(time_varying_confirm_delay_pars)){
            all_preconfirmation_prevalence[indices] <- calculate_preconfirmation_prevalence(onsets, prob_preconfirmation)
          } else {
            browser()
            all_preconfirmation_prevalence[indices] <- calculate_preconfirmation_prevalence_vector(onsets, time_varying_confirm_delay_pars$shape, time_varying_confirm_delay_pars$scale)
          }
          all_cantravel_prevalence[indices] <- all_presymptomatic_prevalence[indices] + all_preconfirmation_prevalence[indices]
          all_disease_prevalence[indices] <- all_presymptomatic_prevalence[indices] + all_symptomatic_prevalence[indices]
        }
        
      }
      ## Once done, combined all of the provinces into one tibble
      all_dat <- data %>% bind_cols(infections=all_infections, onsets=all_onsets, confirmations=all_confirmations,
                                    local_infections=all_local_infections, imported_infections=all_imported_infections)
      if(calculate_prevalence) {
        all_dat <- bind_cols(all_dat, presymptomatic_prevalence=all_presymptomatic_prevalence,
                             symptomatic_prevalence=all_symptomatic_prevalence,
                             preconfirmation_prevalence=all_preconfirmation_prevalence,
                             cantravel_prevalence=all_cantravel_prevalence,
                             disease_prevalence=all_disease_prevalence)
      }
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
