#include <Rcpp.h>
using namespace Rcpp;

/////////////////////////////////////
// SYMPTOM ONSET PROBABILITIES
/////////////////////////////////////
//[[Rcpp::export]]
NumericVector calculate_onset_probs(int tmax, double weibull_alpha, double weibull_sigma){
  NumericVector probs(tmax+1);
  for(int t = 0; t <= tmax; ++t){
    probs[t] = R::pweibull(t+1, weibull_alpha, weibull_sigma, true, false) - R::pweibull(t, weibull_alpha, weibull_sigma, true, false);
  }
  return(probs);
}

//[[Rcpp::export]]
NumericVector calculate_probs_presymptomatic(int tmax, double weibull_alpha, double weibull_sigma){
  NumericVector probs(tmax+1);
  for(int t = 0; t <= tmax; ++t){
    probs[t] = 1 - R::pweibull(t, weibull_alpha, weibull_sigma, true, false);
  }
  return(probs);
}


/////////////////////////////////////
// EXPORTATION PROBABILITIES
/////////////////////////////////////
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
    res(t,t) = probs[t];
    
    // Get prob that you haven't left by each day in the future
    for(int j = t; j <= tmax; ++j){
      tmp = 1;
      // Is the probability that you haven't left on all previous days
      // using probs from today onward
      for(int i = t; i < j ; ++i){
          tmp *= (1-probs[i]);
      }
      res(t,j) = tmp*probs[j];
    }
  }
  return(res);
}

//[[Rcpp::export]]
NumericMatrix prob_leave_pre_symptoms(NumericMatrix leave_matrix, NumericVector presymptom_probs){
  int n_col = leave_matrix.ncol();
  int n_row = leave_matrix.nrow();
  NumericMatrix probs(n_row, n_col);
  
  // For each day of potential infection
  for(int t = 0; t < n_row; ++t){
    // Go forward from this day
    for(int j = t; j < n_col; ++j){
      // What proportion of infections on this day leave by each day in the future given that they haven't
      // yet experienced symptoms?
      probs(t, j) = leave_matrix(t, j)*presymptom_probs[j-t];//(1-R::pweibull(j-t, weibull_alpha, weibull_sigma, true, false));
    }
  }
  return(probs);  
}

//[[Rcpp::export]]
NumericVector prob_leave_pre_symptoms_vector(NumericMatrix leave_matrix,NumericVector presymptom_probs){ //double weibull_alpha, double weibull_sigma){
  NumericMatrix res = prob_leave_pre_symptoms(leave_matrix, presymptom_probs);//weibull_alpha, weibull_sigma);
  return(rowSums(res));
}

/////////////////////////////////////
// IMPORTATION PROBABILITIES
/////////////////////////////////////
// For a given day (row), gives the probability of importing a case from the seed province at each day in the future (column)
//[[Rcpp::export]]
NumericVector prob_daily_arrival(NumericVector export_probs, NumericVector import_probs, int tmax){
  NumericMatrix res(tmax+1, tmax+1);

  double tmp = 1;
  // For each day of potential onset
  for(int t = 0; t <= tmax; ++t){
    // Cannot have left before day 0
    res(t,t) = export_probs[t];
    
    // Get prob that you haven't left by each day in the future
    for(int j = t; j <= tmax; ++j){
      tmp = 1;
      // Is the probability that you haven't left on all previous days
      // using probs from today onward
      for(int i = t; i < j ; ++i){
          tmp *= (1-export_probs[i]);
      }
      res(t,j) = tmp*export_probs[j]*import_probs[j];
    }
  }
  return(res);
}

//[[Rcpp::export]]
NumericMatrix prob_arrive_pre_symptoms(NumericMatrix arrive_matrix, NumericVector presymptom_probs){//double weibull_alpha, double weibull_sigma){
  int n_col = arrive_matrix.ncol();
  int n_row = arrive_matrix.nrow();
  NumericMatrix probs(n_row, n_col);
  
  // For each day of potential infection
  for(int t = 0; t < n_row; ++t){
    // Go forward from this day
    for(int j = t; j < n_col; ++j){
      // What proportion of infections on this day leave by each day in the future given that they haven't
      // yet experienced symptoms?
      probs(t, j) = arrive_matrix(t, j)*presymptom_probs[j-t];//(1-R::pweibull(j-t, weibull_alpha, weibull_sigma, true, false));
    }
  }
  return(probs);  
}

//[[Rcpp::export]]
NumericVector prob_arrive_pre_symptoms_vector(NumericMatrix arrive_matrix, NumericVector presymptom_probs){//double weibull_alpha, double weibull_sigma){
  NumericMatrix res = prob_arrive_pre_symptoms(arrive_matrix, presymptom_probs);//weibull_alpha, weibull_sigma);
  return(rowSums(res));
}
