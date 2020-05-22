#include <Rcpp.h>
using namespace Rcpp;

/////////////////////////////////////
// SYMPTOM ONSET PROBABILITIES
/////////////////////////////////////
//[[Rcpp::export]]
NumericVector calculate_onset_probs_weibull(int tmax, double weibull_alpha, double weibull_sigma){
  NumericVector probs(tmax+1);
  for(int t = 0; t <= tmax; ++t){
    probs[t] = R::pweibull(t+1, weibull_alpha, weibull_sigma, true, false) - R::pweibull(t, weibull_alpha, weibull_sigma, true, false);
  }
  return(probs);
}

//[[Rcpp::export]]
NumericVector calculate_probs_presymptomatic_weibull(int tmax, double weibull_alpha, double weibull_sigma){
  NumericVector probs(tmax+1);
  for(int t = 0; t <= tmax; ++t){
    probs[t] = 1 - R::pweibull(t, weibull_alpha, weibull_sigma, true, false);
  }
  return(probs);
}

//[[Rcpp::export]]
NumericVector calculate_onset_probs_lnorm(int tmax, double par1, double par2){
  NumericVector probs(tmax+1);
  for(int t = 0; t <= tmax; ++t){
    probs[t] = R::plnorm(t+1, par1, par2, true, false) - R::plnorm(t, par1, par2, true, false);
  }
  return(probs);
}

//[[Rcpp::export]]
NumericVector calculate_probs_presymptomatic_lnorm(int tmax, double par1, double par2){
  NumericVector probs(tmax+1);
  for(int t = 0; t <= tmax; ++t){
    probs[t] = 1 - R::plnorm(t, par1, par2, true, false);
  }
  return(probs);
}

//[[Rcpp::export]]
NumericVector calculate_probs_preconfirmation(int tmax, double alpha, double scale){
  NumericVector probs(tmax + 1);
  for(int t = 0; t <= tmax; ++t){
    probs[t] = 1 - R::pgamma(t, alpha, scale, true, false);
  }
  return(probs);
}


/////////////////////////////////////
// RECOVERY PROBABILITIES
/////////////////////////////////////
//[[Rcpp::export]]
NumericVector calculate_probs_notrecovered(int tmax, double alpha, double scale){
  NumericVector probs(tmax + 1);
  for(int t = 0; t <= tmax; ++t){
    probs[t] = 1 - R::pgamma(t, alpha, scale, true, false);
  }
  return(probs);
}


/////////////////////////////////////
// LOCAL TRANSMISSION PROBABILITIES
/////////////////////////////////////
//[[Rcpp::export]]
NumericVector calculate_serial_interval_probs(int tmax, double alpha, double scale){
  NumericVector probs(tmax+1);
  for(int t = 0; t <= tmax; ++t){
//probs[t] = R::plnorm(t+1, lnorm_mean, lnorm_sd, true, false) - R::plnorm(t, lnorm_mean, lnorm_sd, true, false);
    probs[t] = R::pgamma(t+1, alpha, scale, true, false) - R::pgamma(t, alpha, scale, true, false);
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
NumericMatrix probs_not_left_by_day(NumericVector probs, int tmax){
  NumericMatrix res(tmax+1, tmax+1);
  
  double tmp = 1;
  // For each day of potential onset
  for(int t = 0; t <= tmax; ++t){
    // Cannot have left before day 0
    res(t,t) = 1;
    
    // Get prob that you haven't left by each day in the future
    for(int j = t; j <= tmax; ++j){
      tmp = 1;
      // Is the probability that you haven't left on all previous days
      // using probs from today onward
      for(int i = t; i < j ; ++i){
        tmp *= (1-probs[i]);
      }
      res(t,j) = tmp;
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
NumericMatrix prob_daily_arrival(NumericVector export_probs, NumericVector import_probs, int tmax){
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

//[[Rcpp::export]]
NumericMatrix local_travel_matrix_precalc(NumericMatrix prob_arrival_mat){
  int n_col = prob_arrival_mat.ncol();
  int n_row = prob_arrival_mat.nrow();
  NumericMatrix precalc(n_row, n_col);
  
  for(int t = 0; t < n_row; ++t){
    for(int i = 0; i <= t; ++i){
      for(int j = i; j <= t; ++j){
        precalc(t,i) += prob_arrival_mat(i, j);
      }
    }
  }
  return(precalc);  
}

/////////////////////////////////////
// LOCAL TRANSMISSION
/////////////////////////////////////
//[[Rcpp::export]]
NumericVector calculate_local_cases(NumericMatrix precalc_local, NumericVector infections, 
                                    NumericVector serial_probs, double r_local){
  int tmax = infections.size();
  NumericVector local_cases(tmax);
  // For all times
  for(int t = 0; t < tmax; ++t){
    // Sum up infections that were generated on each day in the past,
    // find the expected number of secondary cases that they generate. 
    // This is the serial interval distribution multiplied by the probability that
    // they are in that location on that day post infection
    for(int i = 0; i <= t; ++i){
      local_cases[t] += infections[i]*serial_probs[t-i]*precalc_local(t,i);
    }
    local_cases[t] *= r_local;
  }
  return(local_cases);
}


/////////////////////////////////////
// PREVALENCE
/////////////////////////////////////

//[[Rcpp::export]]
NumericVector calculate_infection_prevalence_local(NumericVector incidence, NumericVector prob_presymptomatic){
  int tmax = incidence.size() - 1;
  NumericVector prevalence(tmax+1);
  for(int t = 0; t <= tmax; ++t){
    for(int i = 0; i <= t; ++i) {
      prevalence[t] += incidence[i]*prob_presymptomatic[t-i];
    }
  }
  return(prevalence);
}

//[[Rcpp::export]]
NumericVector calculate_infection_prevalence_imported(NumericVector incidence, NumericVector prob_presymptomatic,
                                                      NumericMatrix prob_have_arrived){
  int tmax = incidence.size() - 1;
  NumericVector prevalence(tmax+1);
  for(int t = 0; t <= tmax; ++t){
    for(int i = 0; i <= t; ++i) {
      prevalence[t] += incidence[i]*prob_presymptomatic[t-i]*prob_have_arrived(t,i);
    }
  }
  return(prevalence);
}

//[[Rcpp::export]]
NumericVector calculate_infection_prevalence_hubei(NumericVector incidence, NumericVector prob_presymptomatic, NumericMatrix probs_not_left_by_day){
  int tmax = incidence.size() - 1;
  NumericVector prevalence(tmax+1);
  for(int t = 0; t <= tmax; ++t){
    for(int i = 0; i <= t; ++i) {
      prevalence[t] += incidence[i]*prob_presymptomatic[t-i]*probs_not_left_by_day(i,t);
    }
  }
  return(prevalence);
}

//[[Rcpp::export]]
NumericVector calculate_preconfirmation_prevalence_vector(NumericVector onsets, NumericVector alphas, NumericVector scales){
  int tmax = onsets.size() - 1;
  NumericVector prevalence(tmax+1);
  for(int t = 0; t <= tmax; ++t){
    for(int i = 0; i <= t; ++i) {
      // Confirmation delay uses the parameters from the day of onset
      prevalence[t] += onsets[i]*(1.0 - R::pgamma(t-i, alphas[i], scales[i], true, false));
    }
  }
  return(prevalence);
}

//[[Rcpp::export]]
NumericVector calculate_preconfirmation_prevalence(NumericVector onsets, NumericVector prob_not_confirmed){
  int tmax = onsets.size() - 1;
  NumericVector prevalence(tmax+1);
  for(int t = 0; t <= tmax; ++t){
    for(int i = 0; i <= t; ++i) {
      prevalence[t] += onsets[i]*prob_not_confirmed[t-i];
    }
  }
  return(prevalence);
}

//[[Rcpp::export]]
NumericVector calculate_unrecovered_prevalence(NumericVector onsets, NumericVector prob_not_recovered){
  int tmax = onsets.size() - 1;
  NumericVector prevalence(tmax+1);
  for(int t = 0; t <= tmax; ++t){
    for(int i = 0; i <= t; ++i) {
      prevalence[t] += onsets[i]*prob_not_recovered[t-i];
    }
  }
  return(prevalence);
}
