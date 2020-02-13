#' @export
create_model_func_provinces <- function(parTab, data=NULL, PRIOR_FUNC=NULL,
                                        tmax=NULL, confirm_delay_pars=NULL,
                                        ver="posterior", noise_ver="poisson"){
    par_names <- parTab$names
    par_provinces <- parTab$province
    unique_provinces <- unique(parTab$province)
    ## Start from last province
    unique_provinces <- rev(unique_provinces[unique_provinces != "all"])
    
    
    if (!is.null(data)) {
        data <- data %>% arrange(province, date)
        times <- data %>% filter(province == "1") %>% pull(date)
        cases <- data$n
        bigT <- length(times)
        tmax <- max(times)
    } else {
        if (is.null(tmax)) {
            stop("One of data and tmax must be specified")
        }
        if (ver == "posterior") {
            stop("Data must be specified to calculate the posterior")
        }
        bigT <- tmax + 1
    }

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
        
        #propn_imports <- sum(exp(pars_all[which(par_names == "import_propn")])) - exp(pars_seed["import_propn"])
        
        daily_prob_leaving <- prob_left_pre_sympt(pars_seed["export_prob"], weibull_alpha, weibull_sigma, 100)

        #print("")
        #print(paste0("Sum propn imports: ", propn_imports))
        ## If time-varying parameters not specified, enumerate out the point
        ## estimates
        if (is.null(confirm_delay_pars)) {
            gamma_shape <- shape
            gamma_scale <- scale
            confirm_delay_pars <- tibble(date_onset=0:tmax, shape=gamma_shape, scale=gamma_scale)
        }
        all_dat <- NULL
        ## Need to loop through each province
        for(province in unique_provinces) {
            pars <- pars_all[which(par_provinces %in% c("all",province))]
            
            #print(paste0("Province: ",province))
            #print(paste0("Unmodified propn import: ", pars["import_propn"]))
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
            #print(paste0("Import proportion: ",export_propn))
            
            ## Growth model
            growth_rate <- pars["growth_rate"]
            i0 <- pars["i0"]
            
            ## Import model
            growth_rate_imports <- pars_seed["growth_rate"] ## Important, growth rate from province 1
            imports_stop <- pars["imports_stop"]
            ## Solve model for this province
            res <- calculate_all_incidences(growth_rate, growth_rate_imports, t0, t0_import, i0, export_propn, imports_stop,
                                            weibull_alpha, weibull_sigma, confirm_delay_pars$shape, confirm_delay_pars$scale,
                                            tmax)

            ## Extract the 3 incidence types
            infections <- res$infections
            onsets <- res$onsets
            confirmations <- res$confirmations

            infections <- tibble(date=times, n=infections,province=province,var="infections")
            onsets <- tibble(date=times,n=onsets,province=province,var="onsets")
            confirmations <- tibble(date=times,n=confirmations,province=province,var="confirmations")
            #print(paste0("Province confirmations: ", sum(confirmations$n)))
            ## Combine and save in list
            all_dat[[province]] <- bind_rows(infections, onsets, confirmations)
        }
        ## Once done, combined all of the provinces into one tibble
        all_dat <- do.call("bind_rows", all_dat)
        all_dat <- all_dat %>% arrange(province, date)
                
        ## If version to solve is model, we're done. Otherwise solve likelihood for confirmations
        if(ver == "model"){
            return(all_dat)
        } else {
            confirmed <- all_dat %>% filter(all_dat$var == "confirmations") %>% pull(n)
            if (noise_ver == "poisson") {
                lik <- sum(dpois(x=cases, confirmed, log=TRUE))
            } else {
                lik <- sum(dnbinom(x=cases,mu=confirmed,size=size,log=TRUE))
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
