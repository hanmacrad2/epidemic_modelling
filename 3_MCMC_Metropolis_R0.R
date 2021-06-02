#Description
#Metropolis algorithm with target r0

#Setup
source("simulate_branching.R")
par(mar=c(1,1,1,1))

#Params
num_days = 30 #100
r0 = 3.1
shape_gamma = 6
scale_gamma = 1

#Priors
prior_r0_k = 1
prior_r0_theta = 1


#***********************************
#Likelihood (log)
log_like <- function(y, r0_dash){
  
  #Params
  num_days = length(y)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  logl = 0
  
  for (t in 2:num_days) {
    
    lambda = r0_dash*sum(y[1:t-1]*rev(prob_infect[1:t-1]))
    logl = logl + y[t]*log(lambda) - lambda
    
  }
  
  logl
  
}

#***********************************
#MCMC
MetropolisHastings_r0 <- function(data, n, x0 = 1, sigma_opt, burn_in = 1000) {
  
  #Set up
  r0_vec <- vector('numeric', n)
  r0_vec[1] <- x0
  U <- runif(n)
  
  #MCMC chain
  for(i in 2:n) {
    Y <- r0_vec[i-1] + rnorm(1, sd = sigma_opt) #, mean = 0, sd = sigma_opt)
    if(Y < 0){
      Y = abs(Y)
    }
    
    log_alpha = log_like(data, Y) - log_like(data, r0_vec[i-1]) #Priors cancel:log_prior(theta_dash) - log_prior(theta) = 1 - 1 
    
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      r0_vec[i] <- Y
    } else {
      r0_vec[i] <- r0_vec[i-1]
    }
  }
  
  r0_vec = r0_vec[burn_in:n]
  r0_vec
}

#Apply MCMC
n = 20000
sigma_opt = 1 #2.38 #optimal
data = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
r0_mcmc = MetropolisHastings_r0(data, n, sigma_opt)

#Plots
ts.plot(r0_mcmc, ylab = 'R0', main = 'MCMC chain of R0')

#Plot mean
r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
plot(seq_along(r0_mean), r0_mean, xlab = 'Time', ylab = 'R0', main = 'Mean of r0 - MCMC chain')

#Histogram
hist(r0_mcmc, prob = TRUE)

#Hist
hist1 <- hist(r0_mcmc, breaks = 80)
hist1$counts <- hist1$counts/sum(hist1$counts)
plot(hist1, xlab = 'r0', ylab = 'Density', 
     main = 'Empirical density of r0 - MCMC chain')

#Notes: MCMC
#Actual analytical (Not possible with real data - parameters unknown) alpha <- (dgamma(Y, post_gamma_shape, post_gamma_scale))/ (dgamma(r0_vec[i-1], post_gamma_shape, post_gamma_scale))


#***********************
#Bayesian Analytical plot
#Data
x = simulate_branching(num_days, r0, shape_gamma, scale_gamma)

#Posterior
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

#Plot
a = seq(0.0, 20, by = 0.01)
b = rgamma(n, post_gamma_shape, post_gamma_scale)
hist2 = hist(b, breaks = 100)
hist2$counts <- hist2$counts/sum(hist2$counts)
plot(hist2, xlab = 'r0', ylab = 'Density of gamma(2, 0.5)', #shape, scale
     main = 'Density of r0 - Analytical Bayesian model')


