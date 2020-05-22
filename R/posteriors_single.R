#' @export
create_model_func <- function(parTab, data=NULL, PRIOR_FUNC=NULL,
                              tmax=NULL, ver="posterior", noise_ver="poisson",
                              model_ver="logistic",incubation_ver="lnorm",
                              time_varying_confirm_delay_pars=NULL){
  par_names <- parTab$names
  if (!is.null(data)) {
    data <- data %>% arrange(date)
    times <- data$date
    cases <- data$n
    bigT <- nrow(data)
    tmax <- max(times)
  } else {
    if (is.null(tmax)) {
      stop("One of data and tmax must be specified")
    }
    if (ver == "posterior") {
      stop("Data must be specified to calculate the posterior")
    }
    data <- expand.grid(date=0:tmax)
    data <- as_tibble(data)
    data <- data %>% arrange(date)
    bigT <- tmax + 1
  }
  
  ## i) setup confirmation delay bit
  ## If a time-varying confirmation delay distribution is not provided, generate a fixed one
  if(!is.null(time_varying_confirm_delay_pars)){
    report_delay_mat <- calculate_reporting_delay_matrix(time_varying_confirm_delay_pars$shape, time_varying_confirm_delay_pars$scale)
  }
  
  ## Check if will need to recalculate confirmation delay each iteration
  recalc_confirm_delay <- TRUE
  if(all(parTab[parTab$names %in% c("confirm_delay_shape","confirm_delay_scale"),"fixed"] == 1)){
    recalc_confirm_delay <- FALSE
    #confirm_pars <- gamma_pars_from_mean_sd(pars["confirm_delay_shape"], pars["confirm_delay_scale"]^2)
    #shape <- confirm_pars[[1]]
    #scale <- confirm_pars[[2]]
    shape <- pars["confirm_delay_shape"]
    scale <- pars["confirm_delay_scale"]
    
    ## If time-varying parameters not specified, enumerate out the point estimates
    if (is.null(time_varying_confirm_delay_pars)) {
      report_delay_mat <- calculate_reporting_delay_matrix_constant(shape,scale,tmax)
    }
  }
  
  ## ii) setup incubation period
  ##############################################################
  ## Check if need to recalculate incubation period each iteration
  recalc_incubation_period <- TRUE
  if((incubation_ver == "weibull" & all(parTab[parTab$names %in% c("weibull_alpha","weibull_sigma"),"fixed"] == 1)) |
     (incubation_ver == "lnorm" & all(parTab[parTab$names %in% c("lnorm_incu_par1","lnorm_incu_par2"),"fixed"] == 1))) {
    recalc_incubation_period <- FALSE
    if(incubation_ver == "weibull"){
      weibull_alpha <- pars["weibull_alpha"]
      weibull_sigma <- pars["weibull_sigma"]
      ## For each day with a potential infection onset, get the probability of leaving at some point in the future before
      ## symptom onset
      presymptom_probs <- calculate_probs_presymptomatic_weibull(tmax, weibull_alpha, weibull_sigma)
      onset_probs <- calculate_onset_probs_weibull(tmax, weibull_alpha, weibull_sigma)
    } else {
      incu_par1 <- pars["lnorm_incu_par1"]
      incu_par2 <- pars["lnorm_incu_par2"]
      ## For each day with a potential infection onset, get the probability of leaving at some point in the future before
      ## symptom onset
      presymptom_probs <- calculate_probs_presymptomatic_lnorm(tmax, incu_par1, incu_par2)
      onset_probs <- calculate_onset_probs_lnorm(tmax, incu_par1, incu_par2)
    }
  }
  
    
  import_cases <- rep(0, length(times))
    
  model_func <- function(pars){
    names(pars) <- par_names
    
    ## Negative binomial size
    size <- pars["size"]
    
    ###########################################################
    ## i) get confirmation delay matrix, used for model fitting
    ###########################################################
    if(recalc_confirm_delay){
      ## Gamma distribution
      #confirm_pars <- gamma_pars_from_mean_sd(pars["confirm_delay_shape"], pars["confirm_delay_scale"]^2)
      #shape <- confirm_pars[[1]]
      #scale <- confirm_pars[[2]]
      shape <- pars["confirm_delay_shape"]
      scale <- pars["confirm_delay_scale"]
      
      
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
        weibull_alpha <- pars["weibull_alpha"]
        weibull_sigma <- pars["weibull_sigma"]
        ## For each day with a potential infection onset, get the probability of leaving at some point in the future before
        ## symptom onset
        presymptom_probs <- calculate_probs_presymptomatic_weibull(tmax, weibull_alpha, weibull_sigma)
        onset_probs <- calculate_onset_probs_weibull(tmax, weibull_alpha, weibull_sigma)
      } else {
        incu_par1 <- pars["lnorm_incu_par1"]
        incu_par2 <- pars["lnorm_incu_par2"]
        ## For each day with a potential infection onset, get the probability of leaving at some point in the future before
        ## symptom onset
        presymptom_probs <- calculate_probs_presymptomatic_lnorm(tmax, incu_par1, incu_par2)
        onset_probs <- calculate_onset_probs_lnorm(tmax, incu_par1, incu_par2)
      }
    }
    
    ###########################################################
    ## GROWTH OF INFECTIONS
    ###########################################################          
    ## Get local growth of seed province. Will add/subtract this from later estimates
    ## Exponential growth for model 1, logistic for model 2
    t0 <- pars["t0"]
    i0 <- pars["i0"]
    
    if(model_ver == "exp"){
      growth_rate <- pars["growth_rate"]
      infections <- i0*daily_exp_interval_cpp(growth_rate, tmax, t0)
    } else if(model_ver=="logistic"){
      ## Can define logistic growth rate in terms of K and inflection point
      K <- pars["K"]
      if(parTab[parTab$names == "growth_rate","fixed"] == 1){
        t_switch <- pars["t_switch"] - t0 - pars["incu_mode"]
        growth_rate <- log(pars["K"]-1)/t_switch
      } else {
        growth_rate <- pars["growth_rate"]
        t_switch <- log(pars["K"]-1)/growth_rate
      }
      infections <- daily_sigmoid_interval_cpp(growth_rate, K, tmax, t0)
    } else {
      infections <- solve_seir_model(pars,tmax, pars["N"],tstep=1)$new_infections
    }
    onsets <- calculate_onset_incidence(infections*(1-pars["asymptomatic_propn"]), onset_probs, tmax)
    confirmations <- calculate_confirmation_incidence(onsets, tmax, report_delay_mat)*pars["reporting_propn"]
  
    if(ver == "model"){
      infections_dat <- tibble(n=infections,date=times,var="infections")
      onsets_dat <- tibble(n=onsets, date=times,var="onsets")
      confirmations_dat <- tibble(n=confirmations,date=times,var="confirmations")
      final <- bind_rows(infections_dat, onsets_dat, confirmations_dat)
      final$province <- "1"
      return(final)
    } else {
      if (noise_ver == "poisson") {
        lik <- sum(dpois(x=cases, confirmations, log=TRUE),na.rm=TRUE)
      } else {
        lik <- sum(dnbinom(x=cases,mu=confirmations,size=size,log=TRUE),na.rm=TRUE)
      }
      if (!is.null(PRIOR_FUNC)) {
        lik <- lik + PRIOR_FUNC(pars)
      }
      return(lik)
    }
  }
  return(model_func)
}
