#' @export
logistic_growth <- function(K, i0, r, t0, ts){
  y <- K/(1 + ((K-i0)/i0)*exp(-r*(ts-t0)))
  y
}

#' @export
convert_date_to_integer <- function(date, times){
  which(times == as.POSIXct(date,format="%Y-%m-%d",tz="UTC"))
}

#' @export
get_index_par <- function(chain, index){
    par_names <- colnames(chain)[2:(ncol(chain)-1)]
    par <- as.numeric(chain[chain$sampno == index, 2:(ncol(chain)-1)])
    names(par) <- par_names
    return(par)
}

#' @export
get_quantiles <- function(outputs, var_name="none"){
    quants <- as.data.frame(t(apply(outputs, 2, function(x) quantile(x,c(0.025,0.5,0.975)))))
    colnames(quants) <- c("lower","median","upper")
    quants$x <- 0:(nrow(quants)-1)
    quants$var <- var_name
    return(quants)
}

#' @export
generate_start_tab <- function(par_tab){
    for(i in 1:nrow(par_tab)){
        if(par_tab[i,"fixed"] == 0){
            par_tab[i, "values"] <- runif(1,par_tab[i,"lower_start"], par_tab[i, "upper_start"])
        }
    }
    return(par_tab)        
}
#' @export
prob_left_pre_sympt <- function(geom_prob, weibull_alpha, weibull_sigma, tmax){
  p_tot <- 0
  for(i in seq_len(tmax+1)){
    ## Prob that you leave on day i-1 * prob you didn't leave up to day i-1 * prob you didn't develop symptoms up to day i-1
    p_tot <- p_tot + (geom_prob*(1-geom_prob)^(i-1))*(1-pweibull(i-1, weibull_alpha, weibull_sigma))
  }
  p_tot
}
#' @export
create_many_province_partab <- function(parTab_start, n_provinces,tmax=100){
  parTab_top <- parTab_start[parTab_start$province == "all",]
  provinces <- as.character(seq_len(n_provinces))
  parTab <- NULL
  parTab_province1 <- data.frame(values=c(0,1,0,0),
                                 names=c("growth_rate","i0","t0","local_r"),
                                 province="1",
                                 fixed=c(0,1,0,1),
                                 lower_bound=c(0,0,0,0),
                                 upper_bound=c(0.2,100,31,1),
                                 steps=0.1,
                                 lower_start=c(0.1,1,0,0),
                                 upper_start=c(0.12,1,10,1),stringsAsFactors = FALSE)
  parTab_province1 <- generate_start_tab(parTab_province1)
  for(province in provinces[2:length(provinces)]){
    parTab_province <- data.frame(values=c(0,1,0,0.2),
                                  names=c("growth_rate","i0","t0","local_r"),
                                  province=province,
                                  fixed=c(0,1,0,0),
                                  lower_bound=c(0,0,0,0),
                                  upper_bound=c(0.2,1,40,10),
                                  steps=0.1,
                                  lower_start=c(0.05,1,10,0),
                                  upper_start=c(0.15,1,40,1),stringsAsFactors = FALSE)
    parTab_province <- generate_start_tab(parTab_province)
    parTab <- bind_rows(parTab, parTab_province)
  }
  
  parTab <- bind_rows(parTab_top, parTab_province1, parTab)
  return(parTab)
}


#' @export
calculate_reporting_delay_matrix_constant  <- function(shape, scale, tmax){
    full_dist <- pgamma(tmax+1 - (0:tmax), shape=shape, scale=scale) - pgamma(tmax - (0:tmax), shape=shape, scale=scale)

    use <- matrix(0,nrow=tmax+1, ncol=tmax+1)
    for(t in 0:tmax){
        use[t+1,1:(t+1)] <- full_dist[(tmax-t+1):(tmax+1)]
    }
    use
}

#' @export
weibull_mode <- function(weibull_alpha, weibull_sigma){
  weibull_sigma * ((weibull_alpha-1)/weibull_alpha) ^ (1/weibull_alpha)
}

#' @export
create_export_prob_matrix <- function(total_travellers, wuhan_pop_ini,
                                      export_dat,
                                      date_start, date_end,
                                      index_date_start=as.POSIXct("2020-01-10", format="%Y-%m-%d",tz="UTC"), 
                                      index_date_end=as.POSIXct("2020-01-25", format="%Y-%m-%d",tz="UTC")){
  times <- seq(date_start, date_end,by="1 day")
  ## Total emigration index
  emigration_indices <- export_dat %>% 
    filter(Date >= index_date_start & Date <= index_date_end) %>%
    select(emigration) %>% 
    summarise(total_indices=sum(emigration)) %>% 
    pull(total_indices)
  
  ## How many travellers per index?
  per_index <- total_travellers/emigration_indices
  
  ## Convert indices to actual numbers of travellers,
  ## then get the number that leave Hubei
  export_dat <- export_dat %>%
    mutate(total_left=emigration*per_index,
           total_arrived=immigration*per_index,
           prop_export=total_left*frac_leave_hubei)
  
  ## Now go through and ask, for a given initial Wuhan population size,
  ## how many individuals are in Wuhan at the start of each day?
  total_pop <- numeric(nrow(export_dat))
  total_pop[1] <- wuhan_pop_ini
  for(i in 2:nrow(export_dat)) {
    ## -1 because it's number at the start of yesterday, plus those that arrived yesterday minus those who left yesterday
    total_pop[i] <- total_pop[i-1] + export_dat$total_arrived[i-1] - export_dat$total_left[i-1]
  }
  
  ## Convert number that left Wuhan to proportion that left Wuhan
  export_probs <- export_dat$total_left/total_pop
  export_probs_china <- export_dat$prop_export/total_pop
  
  probs <- numeric(length(times))
  probs[which(times < min(export_dat$Date))] <- export_probs[1]
  probs[match(export_dat$Date, times)] <- export_probs
  
  probs_china <- numeric(length(times))
  probs_china[which(times < min(export_dat$Date))] <- export_probs_china[1]
  probs_china[match(export_dat$Date, times)] <- export_probs_china
  
  return(list(probs=probs,total_pop=total_pop, probs_within_china=probs_china))
}

#' @export
gamma_pars_from_mean_sd <-  function(gamma_mean, gamma_var){
  scale <- gamma_var/gamma_mean
  shape <- gamma_mean/scale
  return(list("shape"=shape,"scale"=scale))
}
