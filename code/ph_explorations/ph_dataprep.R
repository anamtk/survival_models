install.packages('KMsurv')
library(KMsurv)
data('larynx')
#stage:disease stage
#time: time in months from first treatment
#until death or end of study
#age: age at diagnosis
#diagyr: year of diagnosis
#delta: 1 = died, 0 = otherwise

#model variables:
#time: is = time if delta = 1 (death), otherwise NA

#cens[i,]: two-column matrix wtih rows of number 
##of individuals
#uncensored: cens[i,1] = 0
#right-censored: cens[i,1] = censoring time
#left-censored: cens[i,2] = censoring time
#interval-censored: cens[i,] = (cens.low, cens.up)
#cens.low = lower limit; cens.up = 
#upper limit of the observed interval

#is.censored: a help ordinal variable to indicate 
#censoring status
#uncensored: is.censored = 0
#right-censored: is.censored = 1
#left-censored: is.censored = 0
#interval censored: is.censored = 1

K <- 3 # = number of intervals
a <- seq(0, max(larynx$time) + 0.001, length.out = K + 1)
#int.obs: vector that tells us which interval each
#observation is
int.obs <- matrix(data = NA, 
                  nrow = nrow(larynx),
                  ncol = length(a)-1)
d <- matrix(data = NA,
            nrow = nrow(larynx),
            ncol = length(a) - 1)

for(i in 1:nrow(larynx)){
  for(k in 1:(length(a)-1)){
    d[i,k] <- ifelse(time[i] - a[k] > 0, 1, 0) *
      ifelse(a[k+1] - time[i] > 0, 1, 0)
    int.obs[i,k] <- d[i,k]*k
  }
}

int.obs <- rowSums(int.obs)




