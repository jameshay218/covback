#include <Rcpp.h>
using namespace Rcpp;

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
  
  NumericVector y = (1.0/growth) * (times_seed_upper - times_seed_lower);
  times_seed = times_seed + t0_shift;
  overall_y[times_seed] = y;
  return(overall_y);
}

//[[Rcpp::export]]
NumericVector calculate_infection_incidence(double growth_rate, double growth_rate_imports, 
                                            int tmax, double t0, double t0_import, double i0, double import_propn,
                                            int imports_stop){
  NumericVector local = daily_exp_interval_cpp(growth_rate,tmax, t0)*i0;
  NumericVector imports = daily_exp_interval_cpp(growth_rate_imports, tmax, t0_import)*import_propn;
  
  if (imports_stop <= tmax){
    IntegerVector blanks = seq(imports_stop, tmax);
    imports[blanks] = 0;
  }
  
  NumericVector total = local + imports;
  return(total);
}
//[[Rcpp::export]]
NumericVector calculate_onset_incidence(NumericVector infections, double weibull_alpha, double weibull_sigma, int tmax){
  NumericVector onsets(tmax+1);
  NumericVector infections_subset;
  IntegerVector seqs;
  
  NumericVector probs;
  
  for(int t = 0; t <= tmax; ++t){
    seqs = seq(0, t);
    infections_subset = infections[seqs];
    probs = Rcpp::pweibull(t-seqs+1, weibull_alpha, weibull_sigma,true,false) - Rcpp::pweibull(t-seqs, weibull_alpha, weibull_sigma,true,false);
    probs = probs*infections_subset;
      
    onsets[t] = sum(probs);
  }
  return(onsets);
}
//[[Rcpp::export]]
double prob_not_symptomatic(double weibull_alpha, double weibull_sigma, double t){
  // Prob that you beame symptomatic in or before this time interval
  double prob = R::pweibull(t+1, weibull_alpha, weibull_sigma,true,false);
  return(1.0-prob);
}

//[[Rcpp::export]]
NumericVector calculate_infection_prevalence(double growth_rate, double growth_rate_imports, 
                                             int tmax, double t0, double t0_import, double i0, double import_propn,
                                             int imports_stop, double weibull_alpha, double weibull_sigma){
  NumericVector incidence = calculate_infection_incidence(growth_rate, growth_rate_imports,
                                                          tmax, t0, t0_import, i0, import_propn, imports_stop);
  
  NumericVector prevalence(tmax+1, 0.0);
  
  for(int t = 0; t <= tmax; ++t){
    for(int x = 0; x <= t; ++x){
      prevalence[t] = prevalence[t] + incidence[x]*prob_not_symptomatic(weibull_alpha, weibull_sigma, t-x);
    }
  }
  return(prevalence);
}


//[[Rcpp::export]]
NumericVector calculate_confirmation_incidence(NumericVector onsets, NumericVector shapes, NumericVector scales, int tmax){
  NumericVector confirmations(tmax+1);
  NumericVector onsets_subset;
  IntegerVector seqs;
  
  NumericVector probs;
  
  double shape1;
  double scale1;
  
  double tmp;
  double prob;
  
  for(int t = 0; t <= tmax; ++t){
    seqs = seq(0, t);
    tmp = 0;
    
    for(int j = 0; j <= t; ++j){
      shape1 = shapes[j];
      scale1 = scales[j];
      prob = R::pgamma(t-j + 1, shape1, scale1, true, false) - R::pgamma(t-j, shape1, scale1, true, false);
      tmp += prob*onsets[j];
    }
    confirmations[t] = tmp;
  }
  return(confirmations);
}

//[[Rcpp::export]]
List calculate_all_incidences(double growth_rate, double growth_rate_imports, 
			      double t0, double t0_import, double i0, double import_propn,
                              int imports_stop, double weibull_alpha, double weibull_sigma, 
                              NumericVector shapes, NumericVector scales,
                              int tmax){
  NumericVector infections = calculate_infection_incidence(growth_rate, growth_rate_imports,
							   tmax, t0, t0_import, i0, import_propn, imports_stop);
  NumericVector onsets = calculate_onset_incidence(infections, weibull_alpha, weibull_sigma, tmax);
  NumericVector confirmations = calculate_confirmation_incidence(onsets, shapes, scales, tmax);
  
  List res = List::create(Named("infections")=infections,Named("onsets")=onsets,Named("confirmations")=confirmations);
  return(res);
}
                          
