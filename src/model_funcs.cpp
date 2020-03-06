#include <Rcpp.h>
using namespace Rcpp;


//' @export
//[[Rcpp::export]]
NumericVector daily_sigmoid_interval_cpp(double growth, double i_max, int tmax, double t0){
  double t0_diff = t0 - floor(t0);
  int t0_shift = floor(t0);
  
  NumericVector overall_y(tmax+1, 0.0);
  
  IntegerVector times_seed = seq(0, tmax-t0_shift);
  NumericVector times_seed_use = as<NumericVector>(times_seed);
  
  NumericVector times_seed_upper = i_max/(1+(i_max-1)*Rcpp::exp((times_seed_use + 1.0 - t0_diff)*-growth));
  NumericVector times_seed_lower = i_max/(1+(i_max-1)*Rcpp::exp((times_seed_use - t0_diff)*-growth));
  
  NumericVector y = (times_seed_upper - times_seed_lower);
  times_seed = times_seed + t0_shift;
  overall_y[times_seed] = y;
  return(overall_y);
}

//' @export
//[[Rcpp::export]]
NumericVector daily_exp_interval_cpp(double growth, int tmax, double t0){
  double t0_diff = t0 - floor(t0);
  int t0_shift = floor(t0);
  
  NumericVector overall_y(tmax+1, 0.0);
  
  IntegerVector times_seed = seq(0, tmax-t0_shift);
  NumericVector times_seed_use = as<NumericVector>(times_seed);
  
  NumericVector times_seed_upper = Rcpp::exp((times_seed_use + 1.0 - t0_diff)*growth);
  NumericVector times_seed_lower = Rcpp::exp((times_seed_use - t0_diff)*growth);
  
  //NumericVector y = (1.0/growth) * (times_seed_upper - times_seed_lower);
  NumericVector y = (times_seed_upper - times_seed_lower);
  times_seed = times_seed + t0_shift;
  overall_y[times_seed] = y;
  return(overall_y);
}

//[[Rcpp::export]]
NumericVector calculate_infection_incidence_time(double growth_rate, double growth_rate_imports, 
                                            int tmax, double t0, double t0_import, double i0, NumericVector import_propns){
  NumericVector local = daily_exp_interval_cpp(growth_rate,tmax, t0)*i0;
  NumericVector imports = daily_exp_interval_cpp(growth_rate_imports, tmax, t0_import)*import_propns;
  NumericVector total = local + imports;
  return(total);
}


//[[Rcpp::export]]
NumericVector calculate_local_from_import_infections(NumericVector imported_infections, NumericVector serial_probs, int tmax){
  NumericVector local(tmax+1);
  NumericVector infections_subset;
  IntegerVector seqs;
  NumericVector probs;
  
  for(int t = 0; t <= tmax; ++t){
    seqs = seq(0, t);
    infections_subset = imported_infections[seqs];
    probs = serial_probs[t-seqs];
    probs = probs*infections_subset;

    local[t] = sum(probs);
  }
  return(local);
}



//[[Rcpp::export]]
NumericVector calculate_onset_incidence(NumericVector infections, NumericVector onset_probs, int tmax){
  NumericVector onsets(tmax+1);
  NumericVector infections_subset;
  IntegerVector seqs;
  NumericVector probs;
  
  for(int t = 0; t <= tmax; ++t){
    seqs = seq(0, t);
    infections_subset = infections[seqs];
    probs = onset_probs[t-seqs];
    probs = probs*infections_subset;
    
    onsets[t] = sum(probs);
  }
  return(onsets);
}

//[[Rcpp::export]]
NumericMatrix calculate_reporting_delay_matrix(NumericVector shapes, NumericVector scales){
  double shape1;
  double scale1;
  int tmax = shapes.size();
  NumericMatrix probs(tmax, tmax);
  double tmp;
  // Time of observation
  for(int t = 0; t < tmax; ++t){
    // Going back in time from now to 0
    for(int j = 0; j <= t; ++j){
      // Use the reporting delay distribution as of time j
      shape1 = shapes[j];
      scale1 = scales[j];
      
      // Probability of reporting between t-j+1 and t-j days after onset
      tmp = R::pgamma(t-j+1, shape1, scale1, true, false) - R::pgamma(t-j, shape1, scale1, true,false);
      probs(t,j) = tmp;
    }
  }
  return(probs);
}


//[[Rcpp::export]]
NumericVector calculate_confirmation_incidence(NumericVector onsets, int tmax,NumericMatrix report_delay_mat){
  NumericVector confirmations(tmax+1);
  NumericVector onsets_subset;
  IntegerVector seqs;
  
  NumericVector probs;
  
  double tmp;
  double prob;
  
  // For each observation time
  for(int t = 0; t <= tmax; ++t){
    seqs = seq(0, t);
    tmp = 0;
    
    // Get all onsets on all preceding days
    for(int j = 0; j <= t; ++j){
      // Probability that an onset t-j days prior to t is reported today
      prob = report_delay_mat(t, j);
      tmp += prob*onsets[j];
    }
    confirmations[t] = tmp;
  }
  return(confirmations);
}
//[[Rcpp::export]]
List calculate_all_incidences(double growth_rate, double t0, double i0, NumericVector import_cases,
                              NumericVector onset_probs, NumericMatrix report_delay_mat,
                              int tmax){
  NumericVector infections;
  
  if(i0 > 0){
    infections = i0*daily_exp_interval_cpp(growth_rate, tmax, t0);
    infections = infections + import_cases;
  } else {
    infections = import_cases;
  }
  
  NumericVector onsets = calculate_onset_incidence(infections, onset_probs, tmax);
  NumericVector confirmations = calculate_confirmation_incidence(onsets, tmax, report_delay_mat);
  
  List res = List::create(Named("infections")=infections,Named("onsets")=onsets,Named("confirmations")=confirmations);
  return(res);
}


//[[Rcpp::export]]
List calculate_all_incidences_logistic(double growth_rate, double t0, double i0,
                                       double K,
                                       NumericVector import_cases,
			      NumericVector onset_probs, NumericMatrix report_delay_mat,
                              int tmax){
  NumericVector infections;
 
  if(i0 > 0){
    infections = daily_sigmoid_interval_cpp(growth_rate, K, tmax, t0);
    infections = infections + import_cases;
  } else {
    infections = import_cases;
  }
  
  NumericVector onsets = calculate_onset_incidence(infections, onset_probs, tmax);
  NumericVector confirmations = calculate_confirmation_incidence(onsets, tmax, report_delay_mat);
  
  List res = List::create(Named("infections")=infections,Named("onsets")=onsets,Named("confirmations")=confirmations);
  return(res);
}




//[[Rcpp::export]]
double prob_not_symptomatic(double weibull_alpha, double weibull_sigma, double t){
  // Prob that you beame symptomatic in or before this time interval
  double prob = R::pweibull(t+1, weibull_alpha, weibull_sigma,true,false);
  return(1.0-prob);
}

