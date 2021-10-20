#**********************************************************************************
#Super spreading model- Monte Carlo
# - Compare dbinom, dpois function to lse 

#Setup
library(MASS)
setwd("~/GitHub/epidemic_modelling")
source("functions.R")

#Params
n = 100 # 10000 #50000
burn_in = 2
shape_gamma = 6
scale_gamma = 1
#Priors
prior_alpha_k = 1
prior_alpha_theta = 1

#*********************************************************************************************
#* Log Likelihood - Using package distribution functions
log_like_ss_R_funcs <- function(x, alphaX, betaX, gammaX){
  
  'Calculate the log likelihood on the simulated data using package distribution functions in R'
  
  #Params
  num_days = length(x)
  shape_gamma = 6
  scale_gamma = 1
  logl = 0
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Over all timepoints
  for (t in 2:num_days) {
    
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    inner_sum = 0
    
    #Sum for all values of y_t
    for (y_t in 0:x[t]){ 
      
      #Log likelihood
      prob_yt = dpois(y_t, alphaX*lambda_t)
      #cat("prob_yt, dpois : ", prob_yt, "\n")
      
      prob_zt = s(x[t] - y_t, betaX*lambda_t, 1/(1 + gammaX))
      #cat("prob_zt, dnbinom : ", prob_zt, "\n")
      
      prob_xt = prob_zt*prob_yt
      inner_sum = inner_sum + prob_xt
      #cat("prob_xt, : ", prob_xt, "\n")
      
    } 
    
    #cat("log(prob_xt), : ", log(prob_xt), "\n")
    logl = logl + log(inner_sum) 
    #cat("logl, : ", logl, "\n")
    
  }
  
  logl
}


#*MCMC - Investigate one parameter at a time
mcmc_ss_compare_functions <- function(data, n, sigma, alphaX, betaX, gammaX, x0 = 1) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Setup
  #MCMC params
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n)
  alpha_vec[1] <- x0; beta_vec[1] <- x0; gamma_vec[1] <- x0
  
  #Results params
  count_accept1 = 0; count_reject1 = 0; count_accept2 = 0
  count_reject2 = 0; count_accept3 = 0; count_reject3 = 0
  vecloglikeI <- vector('numeric', n)
  vecloglikeII <- vector('numeric', n)
  
  # #MCMC chain
  for(i in 2:n) {
    
    #******************************************************
    #alpha
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma) 
    #cat("alpha dash: ", alpha_dash, "\n")
    
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    
    #log alpha
    logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev #+ prior1 - prior2
    vecloglikeI[i] = log_accept_prob
    
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
      count_reject1 = count_reject1 + 1
    }
    
    #II. Using package distribution functions
    loglike_new = log_like_ss_R_funcs(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    loglike_prev = log_like_ss_R_funcs(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    log_accept_prob2 = loglike_new - loglike_prev
    #Print
    #cat("New log accept prob : ", log_accept_prob2, "\n")
    #cat("\n")
    vecloglikeII[i] = log_accept_prob2
    
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
      count_reject1 = count_reject1 + 1
    }
    
    #************************************************************************
    #beta
    beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma)
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }
    
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    prior1 = dgamma(beta_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(beta_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev #+ prior1 - prior2
    
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      beta_vec[i] <- beta_dash
      count_accept2 = count_accept2 + 1
    } else {
      beta_vec[i] <- beta_vec[i-1]
      count_reject2 = count_reject2 + 1
    }
    
    #************************************************************************
    #gamma
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma)

    if(gamma_dash < 0){
      gamma_dash = abs(gamma_dash)
    }

    #Acceptance Probability
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash)
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
    prior1 = dgamma(gamma_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(gamma_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev #+ prior1 - prior2

    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      gamma_vec[i] <- gamma_dash
      count_accept3 = count_accept3 + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
      count_reject3 = count_reject3 + 1
    }
    
  }
  #Final stats
  #alpha
  total_iters1 = count_accept1 + count_reject1
  accept_rate1 = 100*(count_accept1/(count_accept1+count_reject1))
  num_samples1 = count_accept1
  print("Acceptance rate1 = ")
  print(accept_rate1)
  
  #beta
  total_iters2 = count_accept2 + count_reject2
  accept_rate2 = 100*(count_accept2/(count_accept2+count_reject2))
  num_samples2 = count_accept2
  print("Acceptance rate2 = ")
  print(accept_rate2)
  
  #gamma
  total_iters3 = count_accept3 + count_reject3
  accept_rate3 = 100*(count_accept3/(count_accept3+count_reject3))
  num_samples3 = count_accept3
  print("Acceptance rate3 = ")
  print(accept_rate3)
  
  #Burn-in 
  #alpha_vec = alpha_vec[burn_in:n]
  #beta_vec = beta_vec[burn_in:n]
  #gamma_vec = gamma_vec[burn_in:n]
  
  #Plot log likelihoods Comparison (Original + Function)
  plot_cmp = plot_diff_points_comparison(vecloglikeI, vecloglikeII)
  print(plot_cmp)
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, 
              accept_rate1, num_samples1, sigma, 
              accept_rate2, num_samples2, sigma, 
              accept_rate3, num_samples3, sigma))
}

#********
#*Implement
num_days = 50
#lambda params
shape_gamma = 6
scale_gamma = 1
#params
alphaX = 2 #Without ss event, ~r0. 
betaX = 0.05
gammaX = 10
#Epidemic data
sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#Time 
n = 10000
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
mcmc_params_ad = mcmc_ss_compare_functions(sim_data, n, sigma, alphaX, betaX,gammaX)

end_time = Sys.time()
time_elap = end_time - start_time
print('Time elapsed:')
print(time_elap)

#Extract params
alpha_mcmc = mcmc_params_ad[1]
alpha_mcmc = unlist(alpha_mcmc)

beta_mcmc = mcmc_params_ad[2]
beta_mcmc = unlist(beta_mcmc)

gamma_mcmc = mcmc_params_ad[3]
gamma_mcmc = unlist(gamma_mcmc)

#Plot
plot.ts(alpha_mcmc, ylab = 'alpha', main = paste("MCMC Super spreading model, simulated alpha = ", alphaX))
plot.ts(beta_mcmc, ylab = 'beta', main = paste("MCMC Super spreading model, simulated beta = ", betaX))
plot.ts(gamma_mcmc,  ylab = 'gamma', main = paste("MCMC Super spreading model, simulated gamma = ", gammaX))

#alpha mean
alpha_mean = cumsum(alpha_mcmc)/seq_along(alpha_mcmc)
plot2 = plot(seq_along(alpha_mean), alpha_mean, xlab = 'Time', ylab = 'alpha', main = paste("Mean of alpha MCMC chain, True alpha = ",alphaX))
print(plot2)

#beta mean
beta_mean = cumsum(beta_mcmc)/seq_along(beta_mcmc)
plot2 = plot(seq_along(beta_mean), beta_mean, xlab = 'Time', ylab = 'beta', main = paste("Mean of beta MCMC chain, True beta = ",betaX))
print(plot2)

#gamma Mean
gamma_mean = cumsum(gamma_mcmc)/seq_along(gamma_mcmc)
plot2 = plot(seq_along(gamma_mean), gamma_mean, xlab = 'Time', ylab = 'gamma', main = paste("Mean of gamma MCMC chain, True gamma = ",gammaX))
print(plot2)
