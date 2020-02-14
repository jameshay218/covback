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
  int n = probs.size();
  NumericMatrix res(tmax+1, n);
  double tmp = 1;
  // For each day of potential onset
  for(int t = 0; t <= tmax; ++t){
    // Cannot have left before day 0
    res(t,0) = probs[t];
    
    // Get prob that you haven't left by each day in the future
    for(int j = 1; j < n; ++j){
      tmp = 1;
      // Is the probability that you haven't left on all previous days
      // using probs from today onward
      for(int i = 0; i <= j; ++i){
          tmp *= (1-probs[i+t]);
      }
      res(t,j) = tmp*probs[j];
    }
  }
  return(res);
}