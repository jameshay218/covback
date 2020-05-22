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
create_random_effects_rlocal <- function(parTab, province_prior="1",
                                         r_local_mean=0.4, r_local_sd=0.05,
                                         func_ver="norm",
                                         tswitchmean=87,tswitchsd=1){
  subset_parTab <- parTab[parTab$names == "local_r",]
  r_index <- which(subset_parTab$province != province_prior)
  par_names <- parTab$names
  if(func_ver == "norm"){
    prior_func <- function(pars){
      r_local_sd1 <- pars["local_r_sd"]
      r_local_mean1 <- pars["local_r_mean"]
      r_locals <- pars[which(names(pars) == "local_r")]
      lik <- dnorm(r_locals[r_index],r_local_mean1,r_local_sd1,log=TRUE)

      
      return(sum(lik))
    }
  } else {
    prior_func <- function(pars){
      r_local_sd1 <- pars["local_r_sd"]
      r_local_mean1 <- pars["local_r_mean"]
      gamma_pars <- gamma_pars_from_mean_sd(r_local_mean1, r_local_sd1^2)
      r_locals <- pars[which(names(pars) == "local_r")]
      lik <- dgamma(r_locals[r_index],gamma_pars[[1]],scale=gamma_pars[[2]],log=TRUE)
      
      #print(lik)
      ## T switch bit
      growth_rate <- pars[which(names(pars) == "growth_rate")[1]]
      t_switch <- log(pars["K"]-1)/growth_rate
      t_switch <- t_switch + pars[which(names(pars) == "t0")[1]] + 5
      lik2 <- dnorm(t_switch, tswitchmean,tswitchsd,1)
      
      return(sum(lik) + lik2)
    }
  }
  prior_func
}

#' @export
create_random_effects_rlocal2 <- function(parTab, r_local_mean=0.4, r_local_sd=0.05, prior_ver="gamma"){
  gamma_pars <- gamma_pars_from_mean_sd(r_local_mean, r_local_sd^2)
  if(prior_ver == "gamma"){
    prior_func <- function(pars){
      r_locals <- pars[which(names(pars) == "local_r")]
      r_locals <- r_locals[2:length(r_locals)]
      lik <- dgamma(r_locals,gamma_pars[[1]],scale=gamma_pars[[2]],log=TRUE)
      return(sum(lik))
    }
  } else {
    prior_func <- function(pars){
      r_locals <- pars[which(names(pars) == "local_r")]
      r_locals <- r_locals[2:length(r_locals)]
      lik <- dexp(r_locals, r_local_mean,log=TRUE)
      return(sum(lik))
    }
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

#' @export
find_prior_sd <- function(prior_mean, prior_model="norm",
                          prior_lower_quantile, prior_upper_quantile,
                          upper_bound=5){
  
  diff_in_quants <- prior_upper_quantile - prior_lower_quantile
  if(prior_model == "norm"){
    f <- function(prior_sd){
      res <- qnorm(c(0.025,0.975),prior_mean, prior_sd)
      sum((res - c(prior_lower_quantile, prior_upper_quantile))^2)
    }  
  } else {
    f <- function(prior_sd){
      gamma_pars <- gamma_pars_from_mean_sd(prior_mean, prior_sd^2)
      res <- qgamma(c(0.025,0.975),gamma_pars[[1]], scale=gamma_pars[[2]])
      sum((res - c(prior_lower_quantile, prior_upper_quantile))^2)
    } 
  }
  use_par <- optim(1, f, method="Brent",lower=0,upper=upper_bound)$par
  return(use_par)
}

#' @export
find_all_prior_sds <- function(prior_table){
  prior_sds <- numeric(nrow(prior_table))
  for(i in 1:nrow(prior_table)){
    prior_sds[i] <- find_prior_sd(prior_table$mean[i], prior_table$model[i],
                                  prior_table$lower_quantile[i], prior_table$upper_quantile[i],
                                  prior_table$upper_bound_sd[i])  
  }
  return(prior_sds)
}
