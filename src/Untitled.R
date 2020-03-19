prob_arrival_matrix <- prob_arrive_pre_symptoms(arrival_matrices[[2]],presymptom_probs)
serial_probs <- calculate_serial_interval_probs(tmax, 2.29, 1/0.36)

locals <- numeric(nrow(prob_arrival_matrix))
r_local <- 3
for(t in seq_along(locals)){
  outside_tmp <- 0
  for(i in 1:t){
    i_tmp <- y[i]
    travel_tmp <- 0
    for(j in i:t){
      tmp_serial <- 0
      for(k in j:t){
        tmp_serial <- tmp_serial + serial_probs[k-i+1]
      }
      travel_tmp <- travel_tmp + prob_arrival_matrix[i,j]*tmp_serial
    }
    outside_tmp <- outside_tmp + i_tmp*travel_tmp
  }
  locals[t] <- outside_tmp*r_local
}
#plot(locals,type='l')


import_cases2 <- numeric(nrow(prob_arrival_matrix))
for(t in seq_along(import_cases)){
  for(i in 1:t){
    import_cases2[t] <- import_cases2[t] + prob_arrival_matrix[i,t]*y[i]
  }
}

## Get duration that these infected inviduals stay pre-symptomatic




#lines(import_cases)


locals1 <- numeric(nrow(prob_arrival_matrix))
for(t in seq_along(locals1)){
  outside_tmp <- 0
  for(i in 1:t){
    i_tmp <- y[i]
    travel_tmp <- 0
    for(j in i:t){
      travel_tmp <- travel_tmp + prob_arrival_matrix[i,j]*serial_probs[t-i+1]
    }
    outside_tmp <- outside_tmp + i_tmp*travel_tmp
  }
  locals1[t] <- outside_tmp*r_local
}

locals2 <- numeric(nrow(prob_arrival_matrix))
for(t in seq_along(locals2)){
  outside_tmp <- 0
  for(i in 1:t){
    outside_tmp <- outside_tmp + y[i]*serial_probs[t-i+1]*sum(prob_arrival_matrix[i,i:t])
  }
  locals2[t] <- outside_tmp*r_local
}

solve_tmp <- function(){
  precalc_stuff <- matrix(0, nrow=length(locals2),ncol=length(locals2))
  for(t in seq_along(locals2)){
    for(i in 1:t){
      precalc_stuff[t,i] <-  sum(prob_arrival_matrix[i,i:t])
    }
  }
  locals2 <- numeric(nrow(prob_arrival_matrix))
  for(t in seq_along(locals2)){
    outside_tmp <- 0
    for(i in 1:t){
      outside_tmp <- outside_tmp + y[i]*serial_probs[t-i+1]*precalc_stuff[t,i]
    }
    locals2[t] <- outside_tmp*r_local
  }
}

solve_tmp_cpp <- function(){
  x <- local_travel_matrix_precalc(prob_arrival_matrix)
  omg <- calculate_local_cases(x, y, serial_probs, r_local)
}

plot(locals1,type='l',col="red",ylim=c(0,100))
lines(import_cases2,col="blue")
y_tot <- import_cases2 + locals1
lines(y_tot,col="green")
abline(v=which.max(locals1))
#locals <- c(locals)
inc <- c(0,diff(locals))
#lines(inc,col="blue")

