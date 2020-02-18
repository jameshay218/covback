#' @export
simulate_observed_data_provinces <- function(parTab, tmax, treport=tmax, confirm_delay_pars=NULL,
                                             daily_import_probs, daily_export_probs,
                                             add_noise=TRUE, noise_ver="poisson"){
    times <- 0:tmax
    pars_all <- parTab$values
    names(pars_all) <- parTab$names
    ## If time-varying parameters not specified, enumerate out the point
    ## estimates
    if (is.null(confirm_delay_pars)) {
        gamma_shape <- pars_all["shape"]
        gamma_scale <- pars_all["scale"]
        confirm_delay_pars <- tibble(date_onset=0:(tmax*5), shape=gamma_shape, scale=gamma_scale)
    }
    
    ## Simulate the deterministic version
    f <- create_model_func_provinces(parTab, tmax=tmax,confirm_delay_pars=confirm_delay_pars,
                                     daily_import_probs = daily_import_probs, daily_export_probs = daily_export_probs,
                                     ver="model")
    res <- f(parTab$values)

    ## Put into tibbles
    sim_dat <- res %>% filter(var == "infections")
    colnames(sim_dat)[which(colnames(sim_dat) == "date")] <- "date_infection"
    
    ## Enumerate out each individual
    sim_dat <- sim_dat %>% select(-var) %>% uncount(n)

    ## Give each indidividual a random onset date
    sim_dat$delay_sympt <- floor(rweibull(nrow(sim_dat), pars_all["weibull_alpha"], pars_all["weibull_sigma"]))
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
        group_by(date, var, province) %>%
        tally()
    sim_dat_true$var <- paste0(sim_dat_true$var,"_true")
    
    ## Group back into population level for observable data
    sim_dat_observable <- sim_dat %>% filter(date_report <= treport)
    sim_dat_observable <- sim_dat_observable %>% select(-c("delay_sympt","delay_confirm")) %>%
        pivot_longer(c("date_onset","date_report","date_infection"),names_to="var",values_to="date") %>%
        group_by(date, var, province) %>%
        tally()
    sim_dat_observable$var <- paste0(sim_dat_observable$var,"_observable")

    ## Flesh out times with no counts
    full_observable <- expand.grid(var=unique(sim_dat_observable$var), date=times, province=unique(sim_dat_observable$province),
                                   stringsAsFactors=FALSE)
    sim_dat_observable <- sim_dat_observable %>%
        right_join(full_observable) %>%
        mutate(n=ifelse(is.na(n),0,n))

    full_true <- expand.grid(var=unique(sim_dat_true$var), date=times, province=unique(sim_dat_observable$province),
                             stringsAsFactors=FALSE)
    sim_dat_true <- sim_dat_true %>%
        right_join(full_true) %>%
        mutate(n=ifelse(is.na(n),0,n))    

    ## Add noise to grouped observations
    if (add_noise) {
        if (noise_ver == "poisson") {
            sim_dat_observable <- sim_dat_observable %>% mutate(n = ifelse(var=="date_report_observable", rpois(n(), n),n))
        } else {
            sim_dat_observable <- sim_dat_observable %>% mutate(n = ifelse(var=="date_report_observable", rnbinom(n(), mu=n,size=pars_all["size"]),n))
        }
    }

    all_dat <- bind_rows(sim_dat_observable, sim_dat_true) %>% ungroup()

    return(list(linelist=sim_dat, aggregated=all_dat))
}
