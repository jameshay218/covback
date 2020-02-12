#' @export
plot_simulations <- function(sim_dat){
    p <- ggplot(sim_dat) +
        geom_line(aes(x=date,y=n,col=var)) +
        theme_bw() +
        theme(legend.position=c(0.2,0.8))
    p
}

#' @export
plot_model_fit <- function(chain, parTab, data, confirm_delay_pars=NULL,nsamp=1000,
                           add_noise=TRUE, noise_ver="poisson"){
    imports_stop <- parTab[parTab$names == "import_stop","values"]
    dat1 <- data %>% filter(var == "date_report_observable") %>% select(-var)
    quants <- generate_prediction_intervals(chain, parTab, dat1, confirm_delay_pars, nsamp, add_noise, noise_ver)

    p <- ggplot(quants) +
        geom_ribbon(aes(x=x,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
        geom_line(aes(x=x,y=median,col=var)) +
        geom_point(data=data,aes(x=date,y=n,col=var)) +
        geom_vline(xintercept=imports_stop,linetype="dashed") +
        theme_bw() +
        theme(legend.position=c(0.2,0.7))
    p
}

#' @export
generate_prediction_intervals <- function(chain, parTab, data, confirm_delay_pars=NULL,nsamp=1000,
                                          add_noise=TRUE, noise_ver="poisson"){
    model_func <- create_model_func(parTab, data, confirm_delay_pars=confirm_delay_pars, ver="model")

    samps <- sample(unique(chain$sampno), nsamp)
    
    observations_store <- confirmations_store <- onsets_store <- infections_store <- matrix(0, nrow=nsamp, ncol=nrow(data))
    for(i in seq_along(samps)){
        pars <- get_index_par(chain, samps[i])
        res <- model_func(pars)
        infections <- res$infections
        onsets <- res$onsets
        confirmations <- res$confirmations
        if (add_noise) {
            if (noise_ver == "poisson") {
                observations <- rpois(length(confirmations),confirmations)
            } else if (noise_ver == "nbinom") {
                observations <- rnbinom(length(confirmations), mu=confirmations, size=pars["size"])
            } else {
                observations <- observations
            }
        }

        observations_store[i,] <- observations
        infections_store[i,] <- infections
        onsets_store[i,] <- onsets
        confirmations_store[i,] <- confirmations
    }
    
    quants_obs <- get_quantiles(observations_store, "observations")
    quants_infections <- get_quantiles(infections_store, "infections")
    quants_onsets <- get_quantiles(onsets_store, "onsets")
    quants_confirmations <- get_quantiles(confirmations_store, "confirmations")

    all_dat <- bind_rows(quants_obs, quants_infections, quants_onsets, quants_confirmations)
    return(all_dat)    
}

#' @export
plot_posteriors <- function(chain, parTab){
    parTab_estimates <- parTab[parTab$fixed == 0,]
    par_names <- parTab_estimates$names
    chain1 <- chain[,c("sampno",par_names)]

    real_values <- parTab_estimates[,c("names","values")]
    colnames(real_values) <- c("variable","values")

    melt_chain <- chain1 %>% pivot_longer(-sampno, names_to="variable",values_to="value")
    
    estimate_p <- ggplot(melt_chain) + geom_density(aes(x=value)) +
        geom_vline(data=real_values,aes(xintercept=values),linetype="dashed") +
        theme_bw() +
        ggtitle("Dashed line shows true simulated value") +
        facet_wrap(~variable,scales="free")
    estimate_p
    
}
