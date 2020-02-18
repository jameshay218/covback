#include <Rcpp.h>
using namespace Rcpp;

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
