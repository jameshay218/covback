fit_weibull_normal_prior <- function(inc_period_draws){
    f <- function(pars) {
        shape1 <- pars[1]
        shape2 <- pars[2]
        out <- dnorm(inc_period_draws$alpha, mean = shape1, sd = shape2)
        return(sum((out - inc_period_draws$alpha)^2))
    }
    alpha_pars <- optim(c(mean(inc_period_draws$alpha), sd(inc_period_draws$alpha)), f, method = "Nelder-Mead",
                        control = list(abstol = 1e-8, reltol = 1e-8))
    alpha_pars <- alpha_pars$par

    f <- function(pars) {
        shape1 <- pars[1]
        shape2 <- pars[2]
        out <- dnorm(inc_period_draws$sigma, mean = shape1, sd = shape2)
        return(sum((out - inc_period_draws$sigma)^2))
    }
    sigma_pars <- optim(c(mean(inc_period_draws$sigma), sd(inc_period_draws$sigma)), f, method = "Nelder-Mead",
                        control = list(abstol = 1e-8, reltol = 1e-8))
    sigma_pars <- sigma_pars$par
    return(list(alpha_pars=alpha_pars, sigma_pars=sigma_pars))
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
