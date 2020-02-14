#' @export
simulate_observed_data_single <- function(parTab, tmax, treport=tmax, confirm_delay_pars=NULL, add_noise=TRUE, noise_ver="poisson"){
    times <- 0:tmax
    pars <- parTab$values
    names(pars) <- parTab$names
    ## If time-varying parameters not specified, enumerate out the point
    ## estimates
    if (is.null(confirm_delay_pars)) {
        gamma_shape <- pars["shape"]
        gamma_scale <- pars["scale"]
        confirm_delay_pars <- tibble(date_onset=0:(tmax*5), shape=gamma_shape, scale=gamma_scale)
    }
    
    ## Simulate the deterministic version
    f <- create_model_func(parTab, tmax=tmax,confirm_delay_pars=confirm_delay_pars,ver="model")
    res <- f(parTab$values)
    infections <- res$infections
    
    ## Put into tibbles
    sim_dat <- tibble(infections=infections,date_infection=times)

    ## Enumerate out each individual
    sim_dat <- sim_dat %>% uncount(infections)

    ## Give each indidividual a random onset date
    sim_dat$delay_sympt <- floor(rweibull(nrow(sim_dat), pars["weibull_alpha"], pars["weibull_sigma"]))
    sim_dat$date_onset <- sim_dat$date_infection + sim_dat$delay_sympt
    confirm_delay_pars <- confirm_delay_pars %>% filter(date_onset <= max(sim_dat$date_onset))
    
    ## Give each individual a confirmation date
    sim_dat <- sim_dat %>% left_join(confirm_delay_pars)
    sim_dat <- sim_dat %>% mutate(delay_confirm = floor(rgamma(n(), shape=shape,scale=scale)))
    sim_dat$date_report <- sim_dat$date_onset + sim_dat$delay_confirm
    ## We now have the full simulated line list
    
    ## Group back into population level for all simulated data
    sim_dat_true <- sim_dat %>% select(-c("delay_sympt","delay_confirm")) %>%
        pivot_longer(c("date_onset","date_report","date_infection"),names_to="var",values_to="date") %>%
        group_by(date, var) %>%
        tally()
    sim_dat_true$var <- paste0(sim_dat_true$var,"_true")
    
    ## Group back into population level for observable data
    sim_dat_observable <- sim_dat %>% filter(date_report <= treport)
    sim_dat_observable <- sim_dat_observable %>% select(-c("delay_sympt","delay_confirm")) %>%
        pivot_longer(c("date_onset","date_report","date_infection"),names_to="var",values_to="date") %>%
        group_by(date, var) %>%
        tally()
    sim_dat_observable$var <- paste0(sim_dat_observable$var,"_observable")

    ## Flesh out times with no counts
    full_observable <- expand.grid(var=unique(sim_dat_observable$var), date=times)
    sim_dat_observable <- sim_dat_observable %>%
        right_join(full_observable) %>%
        mutate(n=ifelse(is.na(n),0,n))

    full_true <- expand.grid(var=unique(sim_dat_true$var), date=times)
    sim_dat_true <- sim_dat_true %>%
        right_join(full_true) %>%
        mutate(n=ifelse(is.na(n),0,n))    

    ## Add noise to grouped observations
    if (add_noise) {
        if (noise_ver == "poisson") {
            sim_dat_observable <- sim_dat_observable %>% mutate(n = ifelse(var=="date_report_observable", rpois(n(), n),n))
        } else {
            sim_dat_observable <- sim_dat_observable %>% mutate(n = ifelse(var=="date_report_observable", rnbinom(n(), mu=n,size=pars["size"]),n))
        }
    }

    all_dat <- bind_rows(sim_dat_observable, sim_dat_true)

    return(list(linelist=sim_dat, aggregated=all_dat))
}



#' @export
generate_confirmation_delays <- function(final_mean, final_var, start_mean, start_var, mean_rate, var_rate, tmax=100){
    gamma_means <- start_mean*exp(-(0:tmax)*mean_rate) + final_mean
    gamma_vars <- start_var*exp(-(0:tmax)*var_rate) + final_var
    
    to_plot <- tibble(times=0:tmax, mean=gamma_means, variance=gamma_vars)
    to_plot <- to_plot %>% pivot_longer(c("mean","variance"), names_to="variable",values_to="value")
    p <- ggplot(to_plot) +
        geom_line(aes(x=times,y=value,col=variable)) +
        facet_wrap(~variable) +
        theme_bw() +
        theme(legend.position="none")
    
    scales <- gamma_vars/gamma_means
    shapes <- gamma_means/scales
    
    gamma_changes <- tibble(date_onset=0:tmax,shape=shapes,scale=scales)
    return(list(pars=gamma_changes,plot=p))
}