#Maximum Likelihood Estimation**
source("1_simulation.R")
source("2_mle.R")

#Params
r0 = 10.2
num_days = 45
#Data
data = simulate_branching(num_days, r0, shape_gamma, scale_gamma)

#log likelihood 
poi_log_like <- function(r0_opt, y){
  
  #Params
  num_days = length(y)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  logl = 0
  
  for (t in 2:num_days) {
    
    lambda = r0_opt*sum(y[1:t-1]*rev(prob_infect[1:t-1]))
    logl = logl + y[t]*log(lambda) - lambda
    
  }
  
  return(-logl)
  
}

#*********
#Optim
optim(1, poi_log_like, y = data)

result = optim(1, poi_log_like, y = data)

#********************
#R0 - A range of values 
plot_mle <- function(vec_r0){

  vec_mle_r0 = vector('numeric', length(list_r0))
  
  for (i in 1:length(vec_r0)){
    
    data = simulate_branching(num_days, vec_r0[i], shape_gamma, scale_gamma)
    vec_mle_r0[i] = optim(1, poi_log_like, y = data)$par
  }
  
  #Ticks
  ticks = seq(from = 0, to = 20)
  plot(vec_r0, vec_mle_r0, at=ticks, xlab = 'True R0', ylab = 'MLE R0', cex.lab=1.2, col = 'red', main = 'True R0 vs MLE') #, pch = 8)
  #vec_mle_r0

}

#Apply
vec_r0 = seq(from=0.1, to = 20, by = 0.1)
plot_mle(vec_r0)

#**************************************************
#Optim Test
fr <- function(x) {   ## Rosenbrock Banana function
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}

optim(c(-1.2,1), fr)
