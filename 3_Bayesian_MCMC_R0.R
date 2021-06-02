#Description
#MCMC for r0

#Setup
source("simulate_branching.R")
#source("3_Bayesian_MCMC_R0.R")
#par(mar=c(1,1,1,1))

#Params
num_days = 30 #100
r0 = 3.1
shape_gamma = 6
scale_gamma = 1

#Priors
prior_r0_k = 1
prior_r0_theta = 1

#Data
x = simulate_branching(num_days, r0, shape_gamma, scale_gamma)

get_posterior_params <- function(){
  
  #Data + params
  sum_lambda = 0
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  for (t in 2:num_days) {
    
    lambda_t = sum(x[1:t-1]*rev(prob_infect[1:t-1]), na.rm = TRUE)
    sum_lambda = sum_lambda + lambda_t
    
  }
  
  #Posterior
  sum_x = sum(x)
  post_gamma_shape = sum_x + prior_r0_k
  #print(post_gamma_shape)
  post_gamma_scale = 1/(sum_lambda + (1/prior_r0_theta))
  
  
  return(c(post_gamma_shape, post_gamma_scale))
  
}
  
#Apply
post_params = get_posterior_params()
post_gamma_shape = post_params[1]
post_gamma_scale = post_params[2]
  
#MCMC
MetropolisHastings_r0 <- function(n, sigma_opt) {
  
  #Posterior Params
  post_params = get_posterior_params()
  post_gamma_shape = post_params[1]
  post_gamma_scale = post_params[2]
  
  #Set up
  r0_vec <- vector('numeric', n)
  mu_r0 = post_gamma_scale*post_gamma_shape
  r0_vec[1] <- mu_r0 #rnorm(1, mean = mu_r0, sd = Sigma)
  U <- runif(n)
  
  #MCMC chain
  for(i in 2:n) {
    Y <- r0_vec[i-1] + rnorm(1, sd = sigma_opt) #, mean = 0, sd = sigma_opt)
    #if(Y < 0){
      #Y = abs(Y)
    #}
    alpha <- (dgamma(Y, post_gamma_shape, post_gamma_scale))/ (dgamma(r0_vec[i-1], post_gamma_shape, post_gamma_scale))
    
    if(!(is.na(alpha)) && U[i] < alpha) {
      r0_vec[i] <- Y
    } else {
      r0_vec[i] <- r0_vec[i-1]
    }
  }
  r0_vec
}

#Target
n = 10000
sigma_opt = 1 #2.38 #optimal
r0_mcmc = MetropolisHastings_r0(n, sigma_opt)

#Plot
ts.plot(r0_mcmc, ylab = 'R0', main = 'MCMC chain of R0')

#Plot mean
r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
plot(seq_along(r0_mean), r0_mean, xlab = 'Time', ylab = 'R0', main = 'Mean of R0 from MCMC chain')

