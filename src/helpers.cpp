#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector prob_not_left_per_day(NumericVector probs){
  int tmax=  probs.size();
  double tmp;
  NumericVector res(tmax+1);
  // The probability that you haven't left by each day post onset
  res[0] = 1;
  for(int t = 1; t < tmax; ++t){
    tmp = 1;
    // Is the probability to you haven't left on all previous days
    for(int j = 0; j < t; ++j){
      tmp *= (1-probs[j]);
    }
    res[t] = tmp;
  }
  return(res);
}

//[[Rcpp::export]]
NumericMatrix prob_leave_on_day(NumericVector probs, int tmax){
  NumericMatrix res(tmax+1, tmax+1);
  double tmp = 1;
  // For each day of potential onset
  for(int t = 0; t <= tmax; ++t){
    // Cannot have left before day 0
    res(t,0) = probs[t];
    
    // Get prob that you haven't left by each day in the future
    for(int j = 1; j <= tmax; ++j){
      tmp = 1;
      // Is the probability that you haven't left on all previous days
      // using probs from today onward
      for(int i = 0; i < j ; ++i){
          tmp *= (1-probs[i+t]);
      }
      res(t,j) = tmp*probs[j+t];
    }
  }
  return(res);
}

//[[Rcpp::export]]
NumericMatrix prob_leave_pre_symptoms(NumericMatrix leave_matrix, double weibull_alpha, double weibull_sigma){
  int n_col = leave_matrix.ncol();
  int n_row = leave_matrix.nrow();
  NumericMatrix probs(n_row, n_col);
  
  // For each day of potential infection
  for(int t = 0; t < n_row; ++t){
    // Go forward from this day
    for(int j = t; j < n_col; ++j){
      // What proportion of infections on this day leave by each day in the future given that they haven't
      // get experienced symptoms?
      probs(t, j) = leave_matrix(t, j-t)*(1-R::pweibull(j-t, weibull_alpha, weibull_sigma, true, false));
    }
  }
  return(probs);
  
}