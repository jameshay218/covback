fit_weibull_normal_prior <- function(inc_period_draws){
  
  to_fit_alpha <- density(inc_period_draws$alpha)
  to_fit_sigma <- density(inc_period_draws$sigma)
    f <- function(pars) {
        shape1 <- pars[1]
        shape2 <- pars[2]
        out <- dnorm(to_fit_alpha$x, mean = shape1, sd = shape2)
        return(sum((out - to_fit_alpha$y)^2))
    }
    alpha_pars <- optim(c(mean(inc_period_draws$alpha), sd(inc_period_draws$alpha)), f, method = "Nelder-Mead",
                        control = list(abstol = 1e-8, reltol = 1e-8))
    alpha_pars <- alpha_pars$par

    f <- function(pars) {
        shape1 <- pars[1]
        shape2 <- pars[2]
        out <- dnorm(to_fit_sigma$x, mean = shape1, sd = shape2)
        return(sum((out - to_fit_sigma$y)^2))
    }
    sigma_pars <- optim(c(mean(inc_period_draws$sigma), sd(inc_period_draws$sigma)), f, method = "Nelder-Mead",
                        control = list(abstol = 1e-8, reltol = 1e-8))
    sigma_pars <- sigma_pars$par
    return(list(alpha_pars=alpha_pars, sigma_pars=sigma_pars))
}

fit_weibull_gamma_prior <- function(inc_period_draws){
  to_fit_alpha <- density(inc_period_draws$alpha)
  to_fit_sigma <- density(inc_period_draws$sigma)
  f <- function(pars) {
    shape1 <- pars[1]
    shape2 <- pars[2]
    out <- dgamma(to_fit_alpha$x, shape1, shape2)
    return(sum((out - to_fit_alpha$y)^2))
  }
  alpha_pars <- optim(c(10,2), f, method = "Nelder-Mead",
                      control = list(abstol = 1e-8, reltol = 1e-8))
  alpha_pars <- alpha_pars$par
  
  f <- function(pars) {
    shape1 <- pars[1]
    shape2 <- pars[2]
    out <- dgamma(to_fit_sigma$x, shape1, shape2)
    return(sum((out - to_fit_sigma$y)^2))
  }
  sigma_pars <- optim(c(10, 2), f, method = "Nelder-Mead",
                      control = list(abstol = 1e-8, reltol = 1e-8))
  sigma_pars <- sigma_pars$par
  return(list(alpha_pars=alpha_pars, sigma_pars=sigma_pars))
}

fit_lnorm_normal_prior <- function(lnorm_draws){
  f <- function(pars) {
    shape1 <- pars[1]
    shape2 <- pars[2]
    out <- dgamma(serial_interval_draws$param1, shape1, shape2)
    return(sum((out - serial_interval_draws$param1)^2))
  }
  mean_pars <- optim(c(1,1), f, method = "Nelder-Mead",
                      control = list(abstol = 1e-8, reltol = 1e-8))
  mean_pars <- mean_pars$par
  
  f <- function(pars) {
    shape1 <- pars[1]
    shape2 <- pars[2]
    out <- dgamma(serial_interval_draws$param2, shape1, shape2)
    return(sum((out - serial_interval_draws$param2)^2))
  }
  sd_pars <- optim(c(mean(serial_interval_draws$param2), sd(serial_interval_draws$param2)), f, method = "Nelder-Mead",
                      control = list(abstol = 1e-8, reltol = 1e-8))
  sd_pars <- sd_pars$par
  return(list(mean_pars=mean_pars, sd_pars=sd_pars))
}


#' @export
create_incubation_prior <- function(inc_period_draws){
    incu_pars <- fit_weibull_normal_prior(inc_period_draws)
    alpha_pars <- incu_pars$alpha_pars
    sigma_pars <- incu_pars$sigma_pars

    
    prior_func <- function(pars){
        a <- dnorm(pars["weibull_alpha"],alpha_pars[1],alpha_pars[2],TRUE)
        b <- dnorm(pars["weibull_sigma"],sigma_pars[1],sigma_pars[2],TRUE)
        return(a+b)
    }
    prior_func
}


#' @export
create_prior_startdate <- function(parTab, inc_period_draws, t0_mean, t0_sd){
  incu_pars <- fit_weibull_normal_prior(inc_period_draws)
  alpha_pars <- incu_pars$alpha_pars
  sigma_pars <- incu_pars$sigma_pars
  
  par_names <- parTab$names
  par_provinces <- parTab$province
  unique_provinces <- unique(parTab$province)
  ## Start from last province
  unique_provinces <- unique_provinces[unique_provinces != "all"]
  n_provinces <- length(unique_provinces)
  
  prior_func <- function(pars_all){
    names(pars_all) <- par_names
    pars_seed <- pars_all[which(par_provinces == "1")]
    a <- dnorm(pars_all["weibull_alpha"],alpha_pars[1],alpha_pars[2],TRUE)
    if(a > 0){
      print(alpha_pars[1])
      print(alpha_pars[2])
    print(paste0("Weibull alpha: ", pars_all["weibull_alpha"]))
    print(paste0("Val: ", a))
    }
    b <- dnorm(pars_all["weibull_sigma"],sigma_pars[1],sigma_pars[2],TRUE)
    c <- dnorm(pars_seed["t0"], t0_mean, t0_sd, TRUE)
    #d <- dgamma(pars_all["lnorm_mean"],mean_pars[1],mean_pars[2],log=TRUE)
    #e <- dnorm(pars_all["lnorm_sd"],sd_pars[1],sd_pars[2],log=TRUE)
    
    return(a+b+c)
  }
  prior_func
}
