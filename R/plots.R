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
    dat1 <- data %>% filter(var %in% c("date_report_observable", "date_infection_true"))
    quants <- generate_prediction_intervals(chain, parTab, dat1, confirm_delay_pars, nsamp, add_noise, noise_ver)
    quants1 <- quants %>% filter(var %in% c("infections","observations",
                                 "date_onset_true","date_infection_true"))
    p <- ggplot(quants1) +
        geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
        geom_line(aes(x=date,y=median,col=var)) +
        geom_point(data=dat1,aes(x=date,y=n,col=var)) +
        geom_ribbon(data=dat1, aes(x=date,ymin=0,ymax=0,fill=var)) +
        geom_vline(xintercept=imports_stop,linetype="dashed") +
        theme_bw() +
        facet_wrap(~province,scales="free_y")
    p
}

#' @export
generate_prediction_intervals <- function(chain, parTab, data, confirm_delay_pars=NULL,nsamp=1000,
                                          add_noise=TRUE, noise_ver="poisson"){
    model_func <- create_model_func_provinces(parTab, data, confirm_delay_pars=confirm_delay_pars, ver="model")
    par_names <- parTab$names
    
    samps <- sample(unique(chain$sampno), nsamp)
    
    store_all <- NULL
    for(i in seq_along(samps)){
        pars <- get_index_par(chain, samps[i])
        names(pars) <- par_names
        res <- model_func(pars)
        #infections <- res %>% filter(var == "infections") %>% pull(n)
        #onsets <- res %>% filter(var == "onsets") %>% pull(n)
        #confirmations <- res %>% filter(var == "confirmations") %>% pull(n)
        if (add_noise) {
            subset_confirmations <- res %>% filter(var == "confirmations")
            if (noise_ver == "poisson") {
                subset_confirmations <- subset_confirmations %>% 
                    mutate(n = rpois(n(),n)) %>% 
                    mutate(var = "observations")
            } else if (noise_ver == "nbinom") {
                subset_confirmations <- subset_confirmations %>% 
                    mutate(n = rnbinom(n(),mu=n,size=pars["size"])) %>% 
                    mutate(var = "observations")
            } else {
                subset_confirmations <- subset_confirmations %>% mutate(var = "observations")
            }
        }
        res <- bind_rows(res, subset_confirmations)
        res$sampno <- i
        store_all <- bind_rows(store_all, res)
    }
    
    quants <- store_all %>% group_by(province, date, var) %>%
        do(data.frame(t(c(quantile(.$n, probs = c(0.01,0.025,0.25,0.5,0.75,0.975,0.99),na.rm=TRUE),mean(.$n)))))
    colnames(quants) <- c("province","date","var",
                                   "min","lower","midlow","median","midhigh","upper","max","mean")
  
    return(quants)    
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
