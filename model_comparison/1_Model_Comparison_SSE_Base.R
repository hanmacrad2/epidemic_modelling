#Model Comparision

setwd("~/GitHub/epidemic_modelling")
source("epidemic_functions.R")
source("helper_functions.R")

#Epidemic params
num_days = 50
shape_gamma = 6 #Gamma params for infectiousness curve (lambda) distribution
scale_gamma = 1 
seed_count = 1 

##############################
#1. MCMC
rjmcmc_sse_base <- function(data, n, sigma, x0 = 1, prior = TRUE) { #thinning_factor, burn_in
  
  'Returns mcmc samples of alpha & acceptance rate'
  print('MCMC SUPERSPREADING')
  
  #Initialise params
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  alpha_vec[1] <- x0; beta_vec[1] <- x0;
  gamma_vec[1] <- x0; r0_vec[1] <- x0;
  
  #Extract params
  sigma_a = sigma[1]; sigma_b = sigma[2]
  sigma_g = sigma[3]; sigma_bg = sigma[4];
  
  #Result vectors
  count_accept1 = 0; count_accept2 = 0;
  count_accept3 = 0; count_accept4 = 0; 
  count_accept5 = 0; count_accept6 = 0;
  #flag_true = FALSE
  
  #Create folder for mcmc results 
  #folder_mcmc = paste0(folder_results, '/mcmc')
  #ifelse(!dir.exists(file.path(folder_mcmc)), dir.create(file.path(folder_mcmc), recursive = TRUE), FALSE)
  
  #MCMC chain
  for(i in 2:n) { #1:n
    
    #******************************************************
    #ALPHA
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma_a) 
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    #log alpha
    logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev  #+ prior1 - prior
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - alpha_dash + alpha_vec[i-1]
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
    }
    
    #************************************************************************
    #BETA (Only if B > 0)
    if (beta_vec[i-1] > 0){ #& (i > 1)
      
      beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma_b) 
      if(beta_dash < 0){
        beta_dash = abs(beta_dash)
      }
      #loglikelihood
      logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1])
      logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
      log_accept_prob = logl_new - logl_prev
      #Priors
      if (prior){
        log_accept_prob = log_accept_prob - beta_dash + beta_vec[i-1]
      }
      #Metropolis Acceptance Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        count_accept2 = count_accept2 + 1
      } else {
        beta_vec[i] <- beta_vec[i-1]
      }
      
      #************************************************************************
      #GAMMA
      gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma_g) 
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
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        gamma_vec[i] <- gamma_dash
        count_accept3 = count_accept3 + 1
      } else {
        gamma_vec[i] <- gamma_vec[i-1]
      }
      
      #R0
      r0_vec[i] = alpha_vec[i] + beta_vec[i]*gamma_vec[i]
      
      #*****************************************************
      #GAMMA-BETA
      gamma_dash <- gamma_vec[i] + rnorm(1, sd = sigma_bg)#Alter sigma_bg depending on acceptance rate. 
      #Acc rate too big -> Make sigma bigger. Acc rate too small -> make sigma smaller
      if(gamma_dash < 1){ #If less then 1
        gamma_dash = 2 - gamma_dash #abs(gamma_dash)
      }
      #New beta 
      beta_new = (r0_vec[i] - alpha_vec[i])/gamma_dash #Proposing new Gamma AND Beta. Beta_dash = f(R0 & gamma_dash)
      
      if(beta_new >= 0){ #Only accept values of beta > 0
        
        logl_new = log_like_ss_lse(data, alpha_vec[i], beta_new, gamma_dash)
        logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i])
        log_accept_prob = logl_new - logl_prev  
        #Priors
        if (prior){
          log_accept_prob = log_accept_prob - beta_dash + beta_vec[i]
        }
        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          beta_vec[i] <- beta_new
          count_accept4 = count_accept4 + 1
        } 
      }
    }
    
    #RJMCMC Step 
    #Reverse of proposals. Prob of proposing 0 when not 0, = 1. Therefore one of qs is 1. While othre 
    if ((beta_vec[i-1] > 0) | (gamma_vec[i-1] > 0)){ #Proposal. 
      #Not a random walk metropolis - as not using current values to decide the next. As current values are zero - choosing non zero.
      #Indpendence Sampler 
      beta_dash = 0
      gamma_dash = 0
      #Everything cancels
      logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_dash)
      print(paste0('New log_l = ', logl_new))
      logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
      print(paste0('Prev log_l = ', logl_prev))
      log_accept_prob = logl_new - logl_prev 
      
      #Metropolis Step
      print('B 0 proposal')
      unif_var = runif(1)
      print(paste0('unif_var = ', unif_var))
      print(paste0('log_accept_prob = ', log_accept_prob))
      
      if(!(is.na(log_accept_prob)) && log(unif_var) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        gamma_vec[i] <- gamma_dash
        count_accept5 = count_accept5 + 1
      } 
      
    } else { #This acceptance prob will be the reverse of the first version
      #Independence sampler - Propose from prior. If VERY lucky value is accepted to be able to jump between models. 
      beta_dash = rexp(1) #q - the proposal distribution is equal to the prior disribution. Reason: Acc prob = like*prior*q/(like*prior*q)
      gamma_dash = rexp(1) + 1
      
      #Everything cancels
      logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_dash)
      logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
      log_accept_prob = logl_new - logl_prev 
      print('B Independ. proposal')
      unif_var = runif(1)
      print(paste0('unif_var = ', unif_var))
      print(paste0('log_accept_prob = ', log_accept_prob))
      
      #Metropolis Step
      if(!(is.na(log_accept_prob)) && log(unif_var) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        gamma_vec[i] <- gamma_dash
        count_accept6 = count_accept6 + 1
      } 
    }
    
  }
  
  #Final stats
  accept_rate1 = 100*count_accept1/n
  accept_rate2 = 100*count_accept2/n
  accept_rate3 = 100*count_accept3/n
  accept_rate4 = 100*count_accept4/n
  accept_rate5 = 100*count_accept5/n
  accept_rate6 = 100*count_accept6/n
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3, accept_rate4, accept_rate5, accept_rate6))
}


############# --- INSERT PARAMETERS! --- ######################################
n_mcmc = 100 #5500

#### - MCMC params - ######
alphaX = 0.8 
betaX = 0.1 
gammaX = 10 
true_r0 = alphaX + betaX*gammaX
true_r0
model_params = c(alphaX, betaX, gammaX, true_r0)

#MCMC - sigma
sigma_a = 0.4*alphaX
sigma_b = 1.0*betaX 
sigma_g = 0.85*gammaX
sigma_bg = 1.5*gammaX
sigma = c(sigma_a, sigma_b, sigma_g, sigma_bg)
#sigma_base = 0.25 #0.5

#RUN MCMC
start_time = Sys.time()
print('Start time:')
print(start_time)
mcmc_params = rjmcmc_sse_base(sim_data, n, sigma)
end_time = Sys.time()
print('End time:')
print(end_time)
time_elap = get_time(start_time, end_time)

true_r0 = alphaX + betaX*gammaX
true_r0
model_params = c(alphaX, betaX, gammaX, true_r0)

#MCMC - sigma
sigma_a = 0.4*alphaX
sigma_b = 1.0*betaX 
sigma_g = 0.85*gammaX
sigma_bg = 1.5*gammaX
sigma = c(sigma_a, sigma_b, sigma_g, sigma_bg)

print(seed_count)
set.seed(seed_count)

#Epidemic data - Neg Bin
sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#MCMC 
#SOMETHING WEIRD HAPPEING 

#Plotting 
dist_type = 'Neg Bin,'
plot_mcmc_grid(n, sim_data, mcmc_params, true_r0, dist_type, time_elap, seed_count)
  
  
#Seed
#seed_count = seed_count + 1
seed_count
