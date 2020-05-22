SEIR_odes <- function(t, x, params) {
    S <- x[1]
    E <- x[2]
    I <- x[3]
    R <- x[4]
    inc <- x[5]
    
    beta <- params[1]
    gamma <- params[2]
    sigma <- params[3]
    seed_time <- params[4]
    switch_time <- params[5]
    
    if(t > switch_time){
        beta <- params[6]
    }
    
    if(t < seed_time){
        I <- 0
        dS <- -beta*S*I
        dE <- beta*S*I - sigma*E
        dI <- sigma*E - gamma*I
        dR <- gamma*I
        dInc <- beta*S*I
    } else {
        dS <- -beta*S*I
        dE <- beta*S*I - sigma*E
        dI <- sigma*E - gamma*I
        dR <- gamma*I
        dInc <- beta*S*I
    }
    
    
    list(c(dS,dE,dI,dR, dInc))
}

#' @export
solve_seir_model <- function(pars, tmax, N,tstep=0.1){
    ## Times to solve model over
    t <- seq(0,tmax,by=tstep)
    
    i0 <- pars["i0"]
    S <- (N - i0)/N
    E <- 0
    I <- i0/N
    R <- 0
    inc <- 0
    
    beta <- pars["R0"]/pars["gamma"]
    beta_less <- pars["R0_less"]/pars["gamma"]
    gamma <- 1/pars["gamma"]
    sigma <- 1/pars["sigma"]
    t0 <- pars["t0"]
    switch_time <- pars["switch_time"]
    
    ## Note starting conditions for population size - we're working per capita here
    results <- as.data.frame(deSolve::ode(y=c(S=S,E=E,I=I,R=R,inc=0),
                                          times=t, func=SEIR_odes,
                                          parms=c(beta,gamma,sigma, t0,switch_time,beta_less)))
    results$new_infections <- diff(c(0,results$inc))*N
    return(results)
    
}


#' @export
simulate_observed_data_single <- function(pars, tmax, treport=tmax, add_noise=FALSE, noise_ver="poisson",N, n_indivs=1000, model_ver="logistic"){
    times <- 0:tmax
    ## If time-varying parameters not specified, enumerate out the point
    ## estimates
    gamma_shape <- pars["confirm_delay_shape"]
    gamma_scale <- pars["confirm_delay_scale"]
    confirm_delay_pars <- tibble(date_onset_floor=0:(tmax*5), shape=gamma_shape, scale=gamma_scale)
    if(model_ver == "exp"){
        prob_infection <- pars["i0"]*daily_exp_interval_cpp(pars["growth_rate"],tmax, pars["t0"])   
    } else if(model_ver == "logistic"){
        prob_infection <- daily_sigmoid_interval_cpp(pars["growth_rate"],pars["K"],tmax,pars[which(names(pars)=="t0")[1]])/N
    } else {
        prob_infection <- solve_seir_model(pars, tmax, N,tstep=1)$new_infections/N
    }
    ## How many individuals are going to get infected at all?
    total_infections <- sum(prob_infection)*n_indivs
    
    ## Simulate infection time for each individual, starting at t0 = 0
    t_infs <- sample(seq_along(prob_infection), total_infections, prob=prob_infection/sum(prob_infection),replace=TRUE) - 1
    
    ## Enumerate out each individual
    sim_dat <- tibble(date_infection=t_infs)
    sim_dat <- sample_frac(sim_dat, size=1-pars["asymptomatic_propn"])
    ## Give each indidividual a random onset date
    sim_dat$delay_sympt <- rlnorm(nrow(sim_dat), pars["lnorm_incu_par1"], pars["lnorm_incu_par2"])
    sim_dat$date_onset <- sim_dat$date_infection + sim_dat$delay_sympt
    sim_dat$date_onset_floor <- floor(sim_dat$date_onset)
    confirm_delay_pars <- confirm_delay_pars %>% filter(date_onset_floor <= max(sim_dat$date_onset_floor))
    
    ## Give each individual a confirmation date
    sim_dat <- sim_dat %>% left_join(confirm_delay_pars)
    sim_dat <- sim_dat %>% mutate(delay_confirm = rgamma(n(), shape=shape,scale=scale)) %>% select(-date_onset_floor)
    sim_dat$date_confirm <- sim_dat$date_onset + sim_dat$delay_confirm
    sim_dat <- sample_frac(sim_dat, size=pars["reporting_propn"])
    
    ## We now have the full simulated line list
    
    ## Group back into population level for all simulated data
    sim_dat_true <- sim_dat %>% select(-c("delay_sympt","delay_confirm")) %>%
        pivot_longer(c("date_onset","date_confirm","date_infection"),names_to="var",values_to="date")
    
    sim_dat_true$date <- floor(sim_dat_true$date)
    sim_dat_true <- sim_dat_true %>%
        group_by(date, var) %>%
        tally()
    sim_dat_true$var <- paste0(sim_dat_true$var,"_true")
    
    ## Group back into population level for observable data
    sim_dat_observable <- sim_dat %>% filter(date_confirm <= treport)
    sim_dat_observable <- sim_dat_observable %>% select(-c("delay_sympt","delay_confirm")) %>%
        pivot_longer(c("date_onset","date_confirm","date_infection"),names_to="var",values_to="date")
    sim_dat_observable$date <- floor(sim_dat_observable$date)
    sim_dat_observable <- sim_dat_observable %>%
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
            sim_dat_observable <- sim_dat_observable %>% mutate(n = ifelse(var=="date_confirm_observable", rpois(n(), n),n))
        } else {
            sim_dat_observable <- sim_dat_observable %>% mutate(n = ifelse(var=="date_confirm_observable", rnbinom(n(), mu=n,size=pars["size"]),n))
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
