#Contents:Super spreading Events model- MCMC Inference

#Setup
setwd("~/GitHub/epidemic_modelling/Model_super_spreading")
library(MASS)
library(pracma)
source("functions.R")
source("plotting_functions.R")

#Epidemic params
num_days = 50
#Gamma params for infectiousness curve (lambda) distribution
shape_gamma = 6
scale_gamma = 1 

seed_count = 1
#par(mar=c(1,1,1,1))

################################################################################
# MCMC SSE - ALL THREE PARAMETERS AT ONCE
################################################################################

mcmc_super_spreading <- function(data, n, sigma,  sigma_b, x0 = 1) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Set up
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  
  alpha_vec[1] <- x0; beta_vec[1] <- x0;
  gamma_vec[1] <- x0; r0_vec[1] <- x0;
  
  count_accept1 = 0; count_accept2 = 0; count_accept3 = 0;
  
  #MCMC chain
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
    log_accept_prob = logl_new - logl_prev  
    
    #Metropolis Acceptance Step
    #if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
    if(log(runif(1)) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
    }
    
    #************************************************************************
    #beta
    beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma_b) 
    #cat("Beta dash: ", beta_dash, "\n")
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }
    
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    log_accept_prob = logl_new - logl_prev
    
    #Metropolis Acceptance Step
    if(log(runif(1)) < log_accept_prob) {
      beta_vec[i] <- beta_dash
      count_accept2 = count_accept2 + 1
    } else {
      beta_vec[i] <- beta_vec[i-1]
    }
    
    #************************************************************************
    #gamma
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma) 
    
    if(gamma_dash < 1){
      gamma_dash = 2 - gamma_dash #abs(gamma_dash)
    }
    
    #Acceptance Probability
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash) 
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
    log_accept_prob = logl_new - logl_prev 
    
    #Metropolis Acceptance Step
    if(log(runif(1)) < log_accept_prob) {
      gamma_vec[i] <- gamma_dash
      count_accept3 = count_accept3 + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
    }
    
    r0_vec[i] = alpha_vec[i] + beta_vec[i]*gamma_vec[i]
    
  }
  #Final stats
  #alpha
  accept_rate1 = 100*count_accept1/n
  cat("Acceptance rate1 = ",accept_rate1)
  
  #beta
  accept_rate2 = 100*count_accept2/n
  cat("Acceptance rate2 = ", accept_rate2)
  
  #gamma
  accept_rate3 = 100*count_accept3/n
  cat("Acceptance rate3 = ", accept_rate3)
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3))
}

################################################################################
# MCMC SSE ALL THREE PARAMETERS + PRIORS
################################################################################

mcmc_ss_prior <- function(data, n, sigma,  sigma_b, x0 = 1, prior = TRUE) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Set up
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  alpha_vec[1] <- x0; beta_vec[1] <- x0;
  gamma_vec[1] <- x0; r0_vec[1] <- x0;
  
  #Priors + Likelihoods
  like_a <- vector('numeric', n); like_b <- vector('numeric', n); 
  like_g = vector('numeric', n); prior_gamma <- vector('numeric', n); 
  prior_alpha <- vector('numeric', n); prior_beta <- vector('numeric', n); 
  
  count_accept1 = 0; count_accept2 = 0; count_accept3 = 0;
  
  #MCMC chain
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
    log_accept_prob = logl_new - logl_prev  #+ prior1 - prior2
    
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - alpha_dash + alpha_vec[i-1]
    }
    
    #Metropolis Acceptance Step
    #if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
    if(log(runif(1)) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
    }
    #Store
    like_a[i] = exp(log_like_ss_lse(data,  alpha_vec[i], beta_vec[i-1], gamma_vec[i-1]))
    prior_alpha[i] = exp(-alpha_vec[i])
    
    #************************************************************************
    #beta
    beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma_b) 
    #cat("Beta dash: ", beta_dash, "\n")
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }
    
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    log_accept_prob = logl_new - logl_prev
    
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - beta_dash + beta_vec[i-1]
    }
    
    #Metropolis Acceptance Step
    if(log(runif(1)) < log_accept_prob) {
      beta_vec[i] <- beta_dash
      count_accept2 = count_accept2 + 1
    } else {
      beta_vec[i] <- beta_vec[i-1]
    }
    #Store
    like_b[i] = exp(log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1]))
    prior_beta[i] = exp(-beta_vec[i])
    
    #************************************************************************
    #gamma
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma) 
    
    if(gamma_dash < 1){
      gamma_dash = 2 - gamma_dash #abs(gamma_dash)
    }
    
    #Acceptance Probability
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash) 
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
    log_accept_prob = logl_new - logl_prev 
    
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - gamma_dash + gamma_vec[i-1]
    }
    
    #Metropolis Acceptance Step
    #if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
    if(log(runif(1)) < log_accept_prob) {
      gamma_vec[i] <- gamma_dash
      count_accept3 = count_accept3 + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
    }
    #Store priors 
    like_g[i] = exp(log_like_ss_lse(data,  alpha_vec[i], beta_vec[i], gamma_vec[i]))
    prior_gamma[i] = 1 + exp(-gamma_vec[i])
    
    #R0
    r0_vec[i] = alpha_vec[i] + beta_vec[i]*gamma_vec[i]
    
  }
  #Final stats
  #alpha
  accept_rate1 = 100*count_accept1/n
  cat("Acceptance rate1 = ",accept_rate1)
  
  #beta
  accept_rate2 = 100*count_accept2/n
  cat("Acceptance rate2 = ", accept_rate2)
  
  #gamma
  accept_rate3 = 100*count_accept3/n
  cat("Acceptance rate3 = ", accept_rate3)
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3,
              like_a, like_b, like_g,
              prior_alpha, prior_beta, prior_gamma))
}

############# --- INSERT PARAMETERS! --- ######################################
alphaX = 0.8 #1.1 #0.8 #1.1 # 0.8 #2 #0.9 #2 #2 #Without ss event, ~r0.
betaX = 0.2 #0.05 #0.2 #0.05 #0.2 #0.2 #0.05 #0.2 #0.05 #0.05
gammaX = 10
true_r0 = alphaX + betaX*gammaX
true_r0
##!!---##############################################################!!---##

set.seed(seed_count)

#Epidemic data - Neg Bin
sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')
sim_data = sim_data

#Epidemic data - Poisson 
#sim_data = simulate_ss_poisson(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
#plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#MCMC 
n = 30000
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
sigma_b = 0.05
prior = FALSE
#mcmc_params = mcmc_super_spreading(sim_data, n, sigma, sigma_b, x0 = 1)
prior = TRUE
mcmc_params = mcmc_ss_prior(sim_data, n, sigma, sigma_b, x0 = 1)
end_time = Sys.time()
time_elap = round(end_time - start_time, 2)
print('Time elapsed:')
print(time_elap)

#Plotting 
dist_type = 'Neg Bin,'
plot_mcmc_results_total(sim_data, mcmc_params, true_r0, dist_type, time_elap, seed_count)
plot_mcmc_results(sim_data, mcmc_params, true_r0, dist_type, time_elap, seed_count, prior)
plot_mcmc_x4_priors(sim_data, mcmc_params, true_r0, dist_type, time_elap, seed_count, prior)

#Seed
#Data: Set seed 
# for (i in 1:seed_count) {
#   cat('Seed count; ', i)
#   set.seed(i)
# }

seed_count = seed_count + 1
seed_count
