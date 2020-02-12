#' @export
create_model_func <- function(parTab, data=NULL, PRIOR_FUNC=NULL,
                              tmax=NULL, confirm_delay_pars=NULL,
                              ver="posterior", noise_ver="poisson"){
  par_names <- parTab$names
  
  if (!is.null(data)) {
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
    
    bigT <- tmax + 1
  }
  
  model_func <- function(pars){
    names(pars) <- par_names
    
    ## Growth model
    growth_rate <- pars["growth_rate"]
    i0 <- pars["i0"]
    t0 <- pars["t0"]
    
    ## Import model
    growth_rate_imports <- pars["growth_rate_imports"]
    imports_stop <- pars["imports_stop"]
    import_propn <- pars["import_propn"]
    
    ## Gamma distribution
    shape <- pars["shape"]
    scale <- pars["scale"]
    
    ## If time-varying parameters not specified, enumerate out the point
    ## estimates
    if (is.null(confirm_delay_pars)) {
      gamma_shape <- shape
      gamma_scale <- scale
      confirm_delay_pars <- tibble(date_onset=0:tmax, shape=gamma_shape, scale=gamma_scale)
    }
    
    weibull_alpha <- pars["weibull_alpha"]
    weibull_sigma <- pars["weibull_sigma"]
    
    ## Negative binomial size
    size <- pars["size"]
    
    res <- calculate_all_incidences(growth_rate, growth_rate_imports, t0, 0, i0, import_propn, imports_stop,
                                    weibull_alpha, weibull_sigma, confirm_delay_pars$shape, confirm_delay_pars$scale,
                                    tmax)
    if(ver == "model"){
      return(res)
    } else {
      confirmed <- res$confirmations
      if (noise_ver == "poisson") {
        lik <- sum(dpois(x=cases, confirmed, log=TRUE))
      } else {
        lik <- sum(dnbinom(x=cases,mu=confirmed,size=size,log=TRUE))
      }
      if (!is.null(PRIOR_FUNC)) {
        lik <- lik + PRIOR_FUNC(pars)
      }
      return(lik)
    }
  }
  return(model_func)
}
