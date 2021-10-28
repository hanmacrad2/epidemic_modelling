#Description
#Super spreading model- Monte Carlo
library(MASS)

#Setup
setwd("~/GitHub/epidemic_modelling/Model_super_spreading")
source("functions.R")
#par(mar=c(1,1,1,1))

#*MCMC - Investigate one parameter at a time
mcmc_ss_investigate_r0_total <- function(data, n, sigma, alphaX, betaX, gammaX, flag_infer, true_R0, x0 = 1) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Initialise parameters
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_total_vec <- vector('numeric', n)
  r0_total_vec[1] <- x0
  count_accept = 0
  
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
        count_accept = count_accept + 1
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
      count_accept = count_accept + 1
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
    log_accept_prob = logl_new - logl_prev #+ prior1 - prior2
    
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      gamma_vec[i] <- gamma_dash
      count_accept = count_accept + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
    }
    } else {
      gamma_vec[i] <- gamma_vec[i - 1]
    }
    
    #Calculate r0 value
    r0_total_vec[i] = alpha_vec[i] + beta_vec[i]*gamma_vec[i]
  }
  #Final stats - Acceptance Rate
  accept_rate = 100*count_accept/n
  data_10_pc = 0.1*n
  
  #Mean of parameter infered
  if (flag_infer == 'alpha') {
    param_mean_val = mean(alpha_vec[n-data_10_pc:n])
  } else if (flag_infer == 'beta') {
    param_mean_val = mean(beta_vec[n-data_10_pc:n])
  } else if (flag_infer == 'gamma') {
    param_mean_val = mean(gamma_vec[n-data_10_pc:n])
  }
  
  #Results - Create dataframe
  df_results <- data.frame(
    alpha = alphaX,
    beta = betaX,
    gamma = gammaX,
    param_infer = flag_infer,
    param_mean = param_mean_val,
    R0 = true_R0, 
    R0_mean_MCMC = mean(r0_total_vec[n-data_10_pc:n]),
    accept_rate = accept_rate) 
  
  print(df_results)
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_total_vec, accept_rate))
}

#********
#Implement
num_days = 50

#Gamma params for infectiousness curve (lambda) distribution
shape_gamma = 6
scale_gamma = 1 

#INSERT PARAMETERS
param_infer = 'alpha'
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
mcmc_params_ad = mcmc_ss_investigate_r0_total(sim_data, n, sigma, alphaX, betaX, gammaX, param_infer, true_tot_r0)
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

accept_rate = mcmc_params_ad[5]

#*****************************
#Plot results
plot_mcmc_results <- function(param_infer, dist_type, true_tot_r0, param_mcmc, r0_total_vec){
  
  #Sim data
  plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste("Daily Infections, Super Spreading Events ", dist_type, ", r0 = ", true_tot_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

  #Paramter (alpha, beta, gamma)
  plot.ts(param_mcmc, xlab = 'Time', ylab = 'alpha',
          main = paste("MCMC chain of ", param_infer),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Paramter mean
  dat_10_pc = 0.1*n
  param_mean_val = mean(param_mcmc[n-dat_10_pc:n])
  param_mean = cumsum(param_mcmc)/seq_along(param_mcmc)
  plot(seq_along(param_mean), param_mean, xlab = 'Time', ylab = param_infer,
       main = paste("Mean of ", param_infer, ". Final Mean = ", round(param_mean_val, 2), ", True value = ", alphaX),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #R0
  plot.ts(r0_total_mcmc,  xlab = 'Time', ylab = 'r0',
          main = paste("R0 total, MCMC Super spreading model, true total R0 = ", true_tot_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #R0 mean
  r0_tot_mean = cumsum(r0_total_mcmc)/seq_along(r0_total_mcmc)
  plot2 = plot(seq_along(r0_tot_mean), r0_tot_mean,
               xlab = 'Time', ylab = 'r0 total',
               main = paste("Mean of R0 total MCMC chain, True R0 total = ", true_tot_r0),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)

  #Histogram
  hist(r0_total_mcmc, freq = FALSE, xlab = 'R0 total', ylab = 'Density', 
       main = 'Empirical density of R0 total MCMC samples',
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
}

#Apply
dist_type = 'Neg Bin dist'
plot_mcmc_results(param_infer, dist_type, true_tot_r0, alpha_mcmc, r0_total_mcmc)

#Done
#alpha 1, 