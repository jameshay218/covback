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
