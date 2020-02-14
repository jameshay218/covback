#' @export
create_model_func_provinces <- function(parTab, data=NULL, PRIOR_FUNC=NULL,
                                        tmax=NULL, confirm_delay_pars=NULL,
                                        daily_import_probs=NULL, daily_export_probs=NULL,
                                        ver="posterior", noise_ver="poisson"){
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

    ## If time-varying parameters not specified, enumerate out the point
    ## estimates
    ## Move into model func call if need to estimate these...
    if (is.null(confirm_delay_pars)) {
      gamma_shape <- shape
      gamma_scale <- scale
      confirm_delay_pars <- tibble(date_onset=0:tmax, shape=gamma_shape, scale=gamma_scale)
    }
    report_delay_mat <- calculate_reporting_delay_matrix(confirm_delay_pars$shape, confirm_delay_pars$scale)
    
    model_func <- function(pars_all) {
        names(pars_all) <- par_names
        pars_seed <- pars_all[which(par_provinces == "1")]
        
        ## Gamma distribution
        shape <- pars_all["shape"]
        scale <- pars_all["scale"]

        ## Negative binomial size
        size <- pars_all["size"]

        ## Incubation period    
        weibull_alpha <- pars_all["weibull_alpha"]
        weibull_sigma <- pars_all["weibull_sigma"]
        
        
        daily_prob_leaving <- prob_left_pre_sympt(pars_seed["export_prob"], weibull_alpha, weibull_sigma, 100)

        onset_probs <- calculate_onset_probs(tmax, weibull_alpha, weibull_sigma)
 
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
                export_propn <- -daily_prob_leaving
                t0 <- pars_seed["t0"]
                t0_import <- pars_seed["t0"]
            } else {
                t0_import <- pars_seed["t0"]
                t0 <- pars["t0"] + t0_import
                t0 <- min(t0, tmax-1)
                export_propn <- daily_prob_leaving * (pars["import_propn"])#/propn_imports
            }
            
            ## Growth model
            growth_rate <- pars["growth_rate"]
            i0 <- pars["i0"]
            
            ## Import model
            growth_rate_imports <- pars_seed["growth_rate"] ## Important, growth rate from province 1
            imports_stop <- pars["imports_stop"]
            ## Solve model for this province
            #res <- calculate_all_incidences(growth_rate, growth_rate_imports, t0, t0_import, i0, export_propn, imports_stop,
            #                                weibull_alpha, weibull_sigma, confirm_delay_pars$shape, confirm_delay_pars$scale,
            #                                tmax)
            res <- calculate_all_incidences_new(growth_rate, growth_rate_imports, t0, t0_import, i0, export_propn, imports_stop,
                                            onset_probs, report_delay_mat,
                                            tmax)
            ## Extract the 3 incidence types
            infections <- res$infections
            onsets <- res$onsets
            confirmations <- res$confirmations
            
            indices <- (bigT*(index - 1) + 1):(bigT*index)
            index <- index + 1
            all_infections[indices] <- infections
            all_onsets[indices] <- onsets
            all_confirmations[indices] <- confirmations
            
            #infections <- tibble(date=times, n=infections,province=province,var="infections")
            #onsets <- tibble(date=times,n=onsets,province=province,var="onsets")
            #confirmations <- tibble(date=times,n=confirmations,province=province,var="confirmations")
            
            ## Combine and save in list
            #all_dat[[province]] <- bind_rows(infections, onsets, confirmations)
        }
         ## Once done, combined all of the provinces into one tibble
        all_dat <- data %>% bind_cols(infections=all_infections, onsets=all_onsets, confirmations=all_confirmations)
                        
        ## If version to solve is model, we're done. Otherwise solve likelihood for confirmations
        if(ver == "model"){
          all_dat <- all_dat %>% pivot_longer(-c("date","province"),names_to="var",values_to="n")
            return(all_dat)
        } else {
            if (noise_ver == "poisson") {
                lik <- sum(dpois(x=cases, all_confirmations, log=TRUE))
            } else {
                lik <- sum(dnbinom(x=cases,mu=all_confirmations,size=size,log=TRUE))
            }
            if (!is.null(PRIOR_FUNC)) {
                lik <- lik + PRIOR_FUNC(pars)
            }
            #print(lik)
            return(lik)
        }
    }
    return(model_func)
}
