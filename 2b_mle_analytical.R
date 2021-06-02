#Maximum Likelihood Estimation**
source("simulation_branching.R")
source("2b_mle_analytical.R")

par(mar=c(1,1,1,1))

#Params
r0 = 10.2
num_days = 30
#Data
data = simulate_branching(num_days, r0, shape_gamma, scale_gamma)

#log likelihood 
mle_analytical <- function(y){
  
  #Params
  num_days = length(y)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  sum_lambda = 0
  
  for (t in 2:num_days) {
    
    lambda_t = sum(y[1:t-1]*rev(prob_infect[1:t-1]))
    sum_lambda = sum_lambda + lambda_t
    
  }
  
  r0_mle = sum(y)/sum(sum_lambda)
  
  r0_mle
  
}


#********************
#**Plot true vs inferred 
#R0 - A range of values 
plot_mle_analytical <- function(vec_r0){
  
  vec_mle_r0 = vector('numeric', length(vec_r0))
  
  for (i in 1:length(vec_r0)){
    
    data = simulate_branching(num_days, vec_r0[i], shape_gamma, scale_gamma)
    vec_mle_r0[i] = mle_analytical(data)
  }
  
  #Ticks
  ticks2 = seq(from = 0, to = 20)
  plot(vec_r0, vec_mle_r0, xlab = 'True R0', ylab = 'MLE R0', col = 'blue', main = 'True R0 vs Analytical MLE') #, pch = 8)
  axis(1, at =ticks2, labels = ticks2)
  axis(2, at =ticks2, labels = ticks2)
  #vec_mle_r0 at=ticks2 cex.lab=1.2
  
}

#Apply
vec_r0 = seq(from=0.1, to = 20, by = 0.1)
plot_mle_analytical(vec_r0)
