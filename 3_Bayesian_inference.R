source("1_simulation.R")
source("3_bayesian_inference.R")

#Data Params
r0 = 4.0
num_days = 20
shape_gamma = 6
scale_gamma = 1

#Data
x = simulate_branching(num_days, r0, shape_gamma, scale_gamma)

#Priors
prior_r0_k = 1
prior_r0_theta = 1

#Posterior
plot_posterior <- function(x, prior_r0_k, prior_r0_theta){
  
  #Params
  num_days = length(y)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  sum_lambda = 0
  
  for (t in 2:num_days) {
    
    lambda_t = sum(x[1:t-1]*rev(prob_infect[1:t-1]))
    sum_lambda = sum_lambda + lambda_t
    
  }
  
  #Posterior
  sum_x = sum(x)
  vec_r0 = seq(from=0.1, to = 20, by = 0.1)
  post_gamma_shape = sum_x + prior_r0_k
  post_gamma_scale = 1/(sum_lambda + (1/prior_r0_theta))
  density_gamma = dgamma(vec_r0, shape = post_gamma_shape, scale = post_gamma_scale)
  
  #Plot
  plot(vec_r0, density_gamma)
  
}

#Apply
plot_posterior(x, prior_r0_k, prior_r0_theta)

