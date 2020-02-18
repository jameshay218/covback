#include <Rcpp.h>
using namespace Rcpp;

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
NumericVector calculate_confirmation_incidence_old(NumericVector onsets, NumericVector shapes, NumericVector scales, int tmax){
  NumericVector confirmations(tmax+1);
  NumericVector onsets_subset;
  IntegerVector seqs;
  
  NumericVector probs;
  
  double shape1;
  double scale1;
  
  double tmp;
  double prob;
  
  // For each observation time
  for(int t = 0; t <= tmax; ++t){
    seqs = seq(0, t);
    tmp = 0;
    
    // Get all onsets on all preceding days
    for(int j = 0; j <= t; ++j){
      shape1 = shapes[j];
      scale1 = scales[j];
      
      // Probability that an onset t-j days prior to t is reported today
      prob = R::pgamma(t-j + 1, shape1, scale1, true, false) - R::pgamma(t-j, shape1, scale1, true, false);
      tmp += prob*onsets[j];
    }
    confirmations[t] = tmp;
  }
  return(confirmations);
}

//[[Rcpp::export]]
NumericVector calculate_onset_incidence_old(NumericVector infections, double weibull_alpha, double weibull_sigma, int tmax){
  NumericVector onsets(tmax+1);
  NumericVector infections_subset;
  IntegerVector seqs;
  NumericVector probs;
  
  for(int t = 0; t <= tmax; ++t){
    seqs = seq(0, t);
    infections_subset = infections[seqs];
    //probs = onset_probs[t-seqs+1];
    probs = Rcpp::pweibull(t-seqs+1, weibull_alpha, weibull_sigma,true,false) - Rcpp::pweibull(t-seqs, weibull_alpha, weibull_sigma,true,false);
    probs = probs*infections_subset;
      
    onsets[t] = sum(probs);
  }
  return(onsets);
}



//[[Rcpp::export]]
List calculate_all_incidences_old(double growth_rate, double growth_rate_imports, 
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
