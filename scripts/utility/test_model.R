lnorm_par1 <- 1.57
lnorm_par2 <- 0.65
report_shape <- 2
report_scale <- 6

## Infection incidence
K <- 100000
## We know when symptom onsets peaked
t_switch_onsets <- 50
ts <- 0:2000

## Real export probs
export_probs <- read_csv("data/export_probs_matched.csv")$export_prob
## Real import probs
import_probs <- read_csv("data/import_probs_matched.csv")
import_probs <- as.matrix(import_probs[,2:ncol(import_probs)])

tmax <- length(export_probs) -1
leave_matrix <- prob_leave_on_day(export_probs, tmax)

## For each province, what's the daily probability of receiving a person from the seed province?
arrival_matrices <- NULL
## First province isn't meaningful
for (i in 1:2){
  arrival_matrices[[i]] <- prob_daily_arrival(export_probs, import_probs[i,], tmax)
}

########################################
## 1. Infection incidence
########################################
## Calculate daily probability of onset
onset_probs <- calculate_onset_probs_lnorm(max(ts), lnorm_par1,lnorm_par2)

## We need to find the inflection point of infection incidence
## that gives an inflection point of symptom onset at the given time
## when combined with the given incubation period distribution
find_tswitch_offset <- function(){
  f <- function(t_unknown){
    ## Infection incidence inflection t_unknown days before symptom onset
    t_switch_prime <- t_switch_onsets - t_unknown
    ## What would growth rate need to be to peak t_switch_prime days from start?
    growth <- log(K-1)/t_switch_prime
    
    ## Daily infection incidence
    y <- daily_sigmoid_interval_cpp(growth, K, tmax, 0)
    ## Daily onset incidence
    onsets <- calculate_onset_incidence(y, onset_probs,tmax)
    
    ## What day do onsets peak?
    t_switch_unknown <- which.max(onsets)
    ## Difference from set switch point
    (t_switch_unknown - t_switch_onsets)^2
  }
  ## Try infection incidence peaking at X days back, starting at 0
  ## Keep incrementing days back by 1 until we find the correct 
  ## infection incidence inflection time
  t_unknown <- 0
  test1 <- f(t_unknown)
  t_unknown <- t_unknown + 1
  test2 <- f(t_unknown)
  ## While we keep getting closer...
  while(test2 < test1){
    test1 <- test2
    t_unknown <- t_unknown + 1
    test2 <- f(t_unknown)
  }
  t_unknown <- t_unknown - 1
  return(t_unknown)
}
t_offset <- find_tswitch_offset()
t_switch <- t_switch_onsets - t_offset

growth <- log(K-1)/(t_switch)

## Daily incidence of new infections
y <- daily_sigmoid_interval_cpp(growth, K, tmax, 0)
## Cumulative incidence of infections
cumu_y <- K/(1+(K-1)*exp(-growth*0:tmax))

plot(cumu_y,type='l',col="blue")
lines(y,col="red")
abline(v=log(K-1)/growth)

## So we're getting incidence at the end of each day
## ie. I(0) = 1, and I(1) = 1.162, so dI(0) = 0.162

########################################
## 2. Checking forward convolution to get symptom onsets
########################################
onsets <- calculate_onset_incidence(y, onset_probs,tmax)


print(paste0("Time between peak infections and peak onsets: ", 
             which.max(onsets) - which.max(y)))

plot(y, type='l',col="red")
lines(onsets,col="blue")

########################################
## 3. Checking forward convolution to get confirmations
########################################
report_delay_univariate <- pgamma(1:(tmax+1), report_shape, scale=report_scale) - 
  pgamma(0:tmax, report_shape, scale=report_scale)

confirmations <- numeric(tmax+1)
for(i in 1:(tmax+1)){
  for(j in i:(tmax+1)){
    confirmations[j] <- confirmations[j] + onsets[i]*report_delay_univariate[j-i+1]
  }
}

report_delay_mat <- calculate_reporting_delay_matrix_constant(report_shape, report_scale, tmax)
#image(report_delay_mat)
confirmations <- calculate_confirmation_incidence(onsets, tmax, report_delay_mat)
lines(confirmations, col="purple")
print(paste0("Time between peak onset and peak confirmations: ", 
             which.max(confirmations) - which.max(onsets)))

########################################
## 4. Checking importation
########################################
presymptom_probs <- calculate_probs_presymptomatic_lnorm(tmax, lnorm_par1, lnorm_par2)
daily_prob_leaving <- prob_leave_pre_symptoms_vector(leave_matrix, presymptom_probs)
daily_prob_arrival <- prob_arrive_pre_symptoms_vector(arrival_matrices[[2]], presymptom_probs)

arrival_matrices_tmp <- arrival_matrices[[2]]
for(i in 1:nrow(arrival_matrices_tmp)){
  arrival_matrices_tmp[i,i:ncol(arrival_matrices_tmp)] <- arrival_matrices_tmp[i,i:ncol(arrival_matrices_tmp)] *
    presymptom_probs[1:(length(presymptom_probs)-i+1)]
}

import_cases <- daily_prob_arrival * y

serial_probs <- calculate_serial_interval_probs(tmax, 2.29, 1/0.36)

import_cases_local <- calculate_local_from_import_infections(import_cases*3, serial_probs, tmax)
import_cases_local[is.na(import_cases_local)] <- 0
plot(import_cases,type='l',col="blue",ylim=c(0,100))
import_cases_overall <- import_cases + import_cases_local
lines(import_cases_local,col='red')
abline(v=which.max(import_cases_local))


res <- calculate_all_incidences(growth, 0, 0, import_cases_overall,
                                onset_probs, report_delay_mat,
                                tmax)
par(mfrow=c(2,1))
plot(res$infections,type='l')
lines(res$onsets,col="red")
lines(res$confirmations,col="blue")

plot(y, type='l')
lines(onsets,col="red")
lines(confirmations,col="blue")
