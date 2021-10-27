#Description
#Super spreading model- Monte Carlo
library(MASS)

#Setup
setwd("~/GitHub/epidemic_modelling/Model_super_spreading")
source("functions.R")
#par(mar=c(1,1,1,1))

#*MCMC - Investigate one parameter at a time
mcmc_ss_investigate_r0_total <- function(data, n, sigma, alphaX, betaX, gammaX, flag_infer, x0 = 1) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Set up
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_total_vec <- vector('numeric', n)
  count_accept1 = 0; count_accept2 = 0; count_accept3 = 0;
  
  #Initialise parameters
  r0_total_vec[1] <- x0
  
  #alpha
  if (flag_infer == 'alpha') {
    alpha_vec[1] <- x0
    beta_vec[1] <- betaX
    gamma_vec[1] <- gammaX
  }
  
  #beta
  if (flag_infer == 'beta') {
    alpha_vec[1] <- alphaX
    beta_vec[1] <- x0
    gamma_vec[1] <- gammaX
  }
  
  #gamma
  if (flag_infer == 'gamma') {
    alpha_vec[1] <- alphaX
    beta_vec[1] <- betaX
    gamma_vec[1] <- x0
  }
  
  # #MCMC chain
  for(i in 2:n) {
    
    #******************************************************
    #alpha
    if (flag_infer == 'alpha') {
      #alpha_vec[1] <- alphaX; beta_vec[1] <- betaX; gamma_vec[1] <- x0
      
      alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma)
      #cat("alpha dash: ", alpha_dash, "\n")
      
      if(alpha_dash < 0){
        alpha_dash = abs(alpha_dash)
      }
      
      #log alpha
      logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
      logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
      #prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
      #prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
      log_accept_prob = logl_new - logl_prev #+ prior1 - prior2
      
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        alpha_vec[i] <- alpha_dash
        count_accept1 = count_accept1 + 1
      } else {
        alpha_vec[i] <- alpha_vec[i-1]
      }
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
    }

    #************************************************************************
    #beta
    if (flag_infer == 'beta') {
    beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma)
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }

    logl_new = log_like_ss_lse(data, alpha_vec[i-1], beta_dash, gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    #prior1 = dgamma(beta_dash, shape = 1, scale = 1, log = TRUE)
    #prior2 = dgamma(beta_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev #+ prior1 - prior2

    if (!(is.na(log_accept_prob)) &&
        log(runif(1)) < log_accept_prob) {
      beta_vec[i] <- beta_dash
      count_accept2 = count_accept2 + 1
    } else {
      beta_vec[i] <- beta_vec[i - 1]
    }
    } else {
      beta_vec[i] <- beta_vec[i - 1]
    }
    
    #************************************************************************
    #gamma
    if (flag_infer == 'gamma') {
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma)
    
    if(gamma_dash < 0){
      gamma_dash = abs(gamma_dash)
    }
    
    #Acceptance Probability
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash)
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
    #prior1 = dgamma(gamma_dash, shape = 1, scale = 1, log = TRUE)
    #prior2 = dgamma(gamma_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev #+ prior1 - prior2
    
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      gamma_vec[i] <- gamma_dash
      count_accept3 = count_accept3 + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
    }
    } else {
      gamma_vec[i] <- gamma_vec[i - 1]
    }
    
    #Calculate r0 value
    r0_total_vec[i] = alpha_vec[i] + beta_vec[i]*gamma_vec[i]
  }
  #Final stats
  #alpha
  accept_rate1 = 100*count_accept1/n
  #beta
  accept_rate2 = 100*count_accept2/n
  
  #gamma
  accept_rate3 = 100*count_accept3/n
  
  #Results - Create dataframe
  df_results <- data.frame(
    accept_rate_alpha = accept_rate1, #unlist(list_accept_rate1),
    n_accepted_a = count_accept,
    accept_rate_beta = accept_rate2, #unlist(list_accept_rate2),
    n_accepted_b = count_accept2, 
    accept_rate_gamma = accept_rate3,
    n_accepted_g = count_accept3) #unlist(list_accept_rate3))
  
    #time_taken = unlist(list_time_taken))
  
  print(df_results)
  
  #Burn-in 
  #alpha_vec = alpha_vec[burn_in:n]
  #beta_vec = beta_vec[burn_in:n]
  #gamma_vec = gamma_vec[burn_in:n]
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_total_vec))
}

#********
#Implement
num_days = 50

#Gamma params for infectiousness curve (lambda) distribution
shape_gamma = 6
scale_gamma = 1 

#INSERT PARAMETERS
flag_infer = 'alpha'
alphaX = 0.8 #2 #0.9 #2 #2 #Without ss event, ~r0.
betaX = 0.2 #0.05 #0.2 #0.05 #0.05
gammaX = 10
true_tot_r0 = alphaX + betaX*gammaX

#Epidemic data
sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#Time 
n = 50000
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
mcmc_params_ad = mcmc_ss_investigate_r0_total(sim_data, n, sigma, alphaX, betaX, gammaX, flag_infer)
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

r0_total_mcmc = mcmc_params_ad[4]
r0_total_mcmc = unlist(r0_total_mcmc)

#Plots
#alpha
plot.ts(alpha_mcmc, ylab = 'alpha', main = paste("MCMC Super spreading model, simulated alpha = ", alphaX))
#alpha mean
alpha_mean = cumsum(alpha_mcmc)/seq_along(alpha_mcmc)
plot2 = plot(seq_along(alpha_mean), alpha_mean, xlab = 'Time', ylab = 'alpha', main = paste("Mean of alpha MCMC chain, True alpha = ",alphaX))
print(plot2)

#beta
plot.ts(beta_mcmc, ylab = 'beta', main = paste("MCMC Super spreading model, simulated beta = ", betaX))
#beta mean
beta_mean = cumsum(beta_mcmc)/seq_along(beta_mcmc)
plot2 = plot(seq_along(beta_mean), beta_mean, xlab = 'Time', ylab = 'beta', main = paste("Mean of beta MCMC chain, True beta = ",betaX))
print(plot2)

#gamma
plot.ts(gamma_mcmc,  ylab = 'gamma', main = paste("MCMC Super spreading model, simulated gamma = ", gammaX))
#gamma Mean
gamma_mean = cumsum(gamma_mcmc)/seq_along(gamma_mcmc)
plot2 = plot(seq_along(gamma_mean), gamma_mean, xlab = 'Time', ylab = 'gamma', main = paste("Mean of gamma MCMC chain, True gamma = ",gammaX))
print(plot2)

#r0
plot.ts(r0_total_mcmc,  ylab = 'r0', main = paste("R0 total - MCMC Super spreading model, true total r0 = ", true_tot_r0))
#r0 mean
r0_tot_mean = cumsum(r0_total_mcmc)/seq_along(r0_total_mcmc)
plot2 = plot(seq_along(r0_tot_mean), r0_tot_mean, xlab = 'Time', ylab = 'r0 total', main = paste("Mean of R0 total MCMC chain, True R0 total = ", true_tot_r0))
print(plot2)

#r0 hist
#hist(r0_total_mcmc, prob = TRUE, breaks = 80,main = paste("Histogram of R0_total MCMC samples, True R0 total = ", true_tot_r0))
#Hist - density
hist2 <- hist(r0_total_mcmc, breaks = 80)
hist2$counts <- hist2$counts/sum(hist2$counts)
hist3 = plot(hist2, xlab = 'R0 total', ylab = 'Density', 
             main = 'Empirical density of R0 total - MCMC samples')
print(hist3)