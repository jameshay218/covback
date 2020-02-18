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
  parTab_province1 <- data.frame(values=c(0,1,0),
                                 names=c("growth_rate","i0","t0"),
                                 province="1",
                                 fixed=c(0,1,0),
                                 lower_bound=c(0,0,0),
                                 upper_bound=c(0.2,100,31),
                                 steps=0.1,
                                 lower_start=c(0.1,1,0),
                                 upper_start=c(0.12,1,10),stringsAsFactors = FALSE)
  parTab_province1 <- generate_start_tab(parTab_province1)
  for(province in provinces[2:length(provinces)]){
    parTab_province <- data.frame(values=c(0,1,0),
                                  names=c("growth_rate","i0","t0"),
                                  province=province,
                                  fixed=c(0,1,0),
                                  lower_bound=c(0,0,0),
                                  upper_bound=c(0.2,1,40),
                                  steps=0.1,
                                  lower_start=c(0.05,1,10),
                                  upper_start=c(0.15,1,40),stringsAsFactors = FALSE)
    parTab_province <- generate_start_tab(parTab_province)
    parTab <- bind_rows(parTab, parTab_province)
  }
  
  parTab <- bind_rows(parTab_top, parTab_province1, parTab)
  return(parTab)
}
