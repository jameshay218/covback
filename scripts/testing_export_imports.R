library(ggplot2)

## Daily probability of leaving the seed province
## Assume exportations stop halfway through
tmax <- 100
export_probs <- c(runif(tmax/2,0.001,0.01),rep(0,tmax/2 +1))

## Let's say 5 provinces - exportations get spread across these
import_probs <- matrix(0, nrow=5, ncol=tmax+1)
import_prob_averages <- runif(5, 0.1,0.8)
for(i in 1:nrow(import_probs)){
  import_probs[i,] <- runif(tmax+1, import_prob_averages[i]-0.05, import_prob_averages[i]+0.05)
}
## Normalize
import_probs_sums <- colSums(import_probs)
import_probs <- t(apply(import_probs, 1, function(x) x/import_probs_sums))
import_probs_tmp <- reshape2::melt(import_probs)
colnames(import_probs_tmp) <- c("location","time","import_prop")
ggplot(import_probs_tmp) + geom_line(aes(x=time,y=import_prop,col=as.factor(location)))

## For each day of possible infection onset, get daily probability of leaving 
## on each day in the future and not before
y_export <- prob_leave_on_day(export_probs, tmax)
image(t(y_export))
## So probability of leaving at all is just rowSums of this
plot(rowSums(y_export), type='l', xlab="Days since t0", ylab="Probability of leaving before event")

## But if can only leave pre symptoms, then need to take into account
## probability that you haven't developed symptoms
## This gives the proportion of infections from each day that end up in other 
x_export <- prob_leave_pre_symptoms_vector(y_export, 2.5, 6)
lines(x_export, col="red")
legend(20, 0.4, legend=c("Probability of leaving at some point in the future",
                         "Probability of leaving before symptom onset"),
       col=c("black","red"),lty=c(1,1),cex=0.8)

## Now get probability of each exportation going to each province
arrival_probs <- matrix(0, nrow=5, ncol=tmax+1)
for(i in 1:nrow(import_probs)){
  tmp <- prob_daily_arrival(export_probs, import_probs[i,], tmax)
  arrival_probs[i,] <- prob_arrive_pre_symptoms_vector(tmp, 2.5, 6)
}
lines(colSums(arrival_probs),col="green")

plot(x_export, type='l', xlab="Day of infection", ylab="Probability of leaving before symptom onset", col="purple")
lines(colSums(arrival_probs),col="purple")
for(i in 1:nrow(arrival_probs)) {
  lines(arrival_probs[i,], col=i)
}

legend(60, 0.03, legend=c("Exportation probability",
                         "Importation probability (one province)"),
       col=c("purple","black"),lty=c(1,1),cex=0.8)

## How long each province's importation prob calculation takes per iteration
times <- microbenchmark::microbenchmark(prob_arrive_pre_symptoms_vector(tmp, 2.5, 6))
## 400 microseconds per solve, for 26 provinces and ~100000 iterations will take about 18 minutes extra
