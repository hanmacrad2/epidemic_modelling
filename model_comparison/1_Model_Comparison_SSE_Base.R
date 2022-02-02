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
#1. MCMC - Loglikelihood
log_like_ss_lse_B0 <- function(x, alphaX, betaX, gammaX, flagB0){
  
  #Params
  num_days = length(x)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)),
                                                                                         shape = shape_gamma, scale = scale_gamma)
  logl = 0 
  
  for (t in 2:num_days) {
    
    #print(t)
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    
    if(x[t] == 0){ #y_t also equal to zero
      
      #L(x_i) for y_t, x_t = 0
      logl = logl -(alphaX*lambda_t) - 
        (betaX*lambda_t*log(gammaX +1))
      
      #print(paste0('logl 1 = ', logl))
      
    } else {
      
      #Terms in inner sum
      inner_sum_vec <- vector('numeric', x[t])
      
      #Check if Beta & gamma == 0
      if (flagB0) {
        
        for (y_t in 0:x[t]){ #Sum for all values of y_t up to x_t
          
          inner_sum_Lx = (-(alphaX*lambda_t) - lfactorial(y_t) + y_t*log(alphaX*lambda_t) +
                            lgamma((x[t] - y_t) + (betaX*lambda_t)) #lgamma(betaX*lambda_t) - 
                          - lfactorial(x[t] - y_t)) #- (betaX*lambda_t*log(gammaX +1)) 
          #+ (x[t] - y_t)*log(gammaX) -(x[t] - y_t)*log(gammaX + 1))
          
          #Check inf
          if (is.infinite(inner_sum_Lx)){
            #print(paste0('inner_sum_Lx is inf:', inner_sum_Lx))
            #print(paste0('yt value = ', y_t))
          } else {
            inner_sum_vec[y_t + 1] = inner_sum_Lx
          }
          
        }
        
        #ELSE IF NOT == 0
      } else {
        
        for (y_t in 0:x[t]){ #Sum for all values of y_t up to x_t
          #print(paste0('y_t  = ', y_t))
          #Store inner L(x_i) term in vector position
          inner_sum_Lx = (-(alphaX*lambda_t) - lfactorial(y_t) + y_t*log(alphaX*lambda_t) +
                            lgamma((x[t] - y_t) + (betaX*lambda_t)) - lgamma(betaX*lambda_t) -
                            lfactorial(x[t] - y_t) - (betaX*lambda_t*log(gammaX +1)) +
                            (x[t] - y_t)*log(gammaX) -(x[t] - y_t)*log(gammaX + 1))
          
          #Check inf
          if (is.infinite(inner_sum_Lx)){
            #print(paste0('inner_sum_Lx is inf:', inner_sum_Lx))
            #print(paste0('yt value = ', y_t))
          } else {
            inner_sum_vec[y_t + 1] = inner_sum_Lx
          }
        }
      }
      
      #Calculate max element in inner vector, for all y_t for a given t, x[t]
      #print(paste0('inner_sum_vec = ', inner_sum_vec))
      lx_max = max(inner_sum_vec)
      #print(paste0('lx_max = ', lx_max))
      
      #Calculate lse
      lse = lx_max + log(sum(exp(inner_sum_vec - lx_max) ))
      #print(paste0('lse = ', lse))
      
      #Add to overall log likelihood 
      logl = logl + lse 
      
    }
    
  }
  
  logl
  
}

##################
#RJMCMC  
rjmcmc_sse_base <- function(data, n, sigma, x0 = 1, prior = TRUE) { #thinning_factor, burn_in
  
  'Returns mcmc samples for sse model w/ rjmcmc & acceptance rates'
  
  #Initialise params
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  alpha_vec[1] <- x0; beta_vec[1] <- x0;
  gamma_vec[1] <- x0; r0_vec[1] <- x0;
  flagB0 = FALSE
  
  #Extract params
  sigma_a = sigma[1]; sigma_b = sigma[2]
  sigma_g = sigma[3]; sigma_bg = sigma[4];
  
  #Result vectors
  count_accept1 = 0; count_accept2 = 0;
  count_accept3 = 0; count_accept4 = 0; 
  count_accept5 = 0; count_accept6 = 0;
  
  #MCMC chain
  for(i in 2:n) {
    
    #******************************************************
    #ALPHA
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma_a) 
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    
    #log alpha
    logl_new = log_like_ss_lse_B0(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1], flagB0)
    logl_prev = log_like_ss_lse_B0(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1], flagB0)
    prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev  #+ prior1 - prior
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - alpha_dash + alpha_vec[i-1]
    }
    
    #Print
    #print(paste0('alpha_dash = ', alpha_dash))
    if (is.na(log_accept_prob)){
      #print(paste0('logl_new = ', logl_new))
      #print(paste0('logl_prev = ', logl_prev))
    } else {
      #print(paste0('log_accept_prob = ', log_accept_prob))
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
    }
    
    #************************************************************************
    #BETA (ONLY IF B > 0)
    if (beta_vec[i-1] > 0){ 
      
      beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma_b) 
      if(beta_dash < 0){
        beta_dash = abs(beta_dash)
      }
      #loglikelihood
      logl_new = log_like_ss_lse_B0(data, alpha_vec[i], beta_dash, gamma_vec[i-1], flagB0)
      logl_prev = log_like_ss_lse_B0(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1], flagB0)
      log_accept_prob = logl_new - logl_prev
      #Priors
      if (prior){
        log_accept_prob = log_accept_prob - beta_dash + beta_vec[i-1]
      }
      #Metropolis Acceptance Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        count_accept2 = count_accept2 + 1
        flagB0 = FALSE
      } else {
        beta_vec[i] <- beta_vec[i-1]
      }
      #Print
      #print(paste0('beta_dash = ', beta_dash))
      #print(paste0('log_accept_prob = ', log_accept_prob))
      
      #************************************************************************
      #GAMMA
      gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma_g) 
      if(gamma_dash < 1){
        gamma_dash = 2 - gamma_dash #abs(gamma_dash)
      }
      #Acceptance Probability
      logl_new = log_like_ss_lse_B0(data, alpha_vec[i], beta_vec[i], gamma_dash, flagB0)
      logl_prev = log_like_ss_lse_B0(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1], flagB0)
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
      
      #Print
      #print(paste0('gamma_dash = ', gamma_dash))
      #print(paste0('log_accept_prob = ', log_accept_prob))
      
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
        
        logl_new = log_like_ss_lse_B0(data, alpha_vec[i], beta_new, gamma_dash, flagB0)
        logl_prev = log_like_ss_lse_B0(data, alpha_vec[i], beta_vec[i], gamma_vec[i], flagB0)
        log_accept_prob = logl_new - logl_prev  
        #Priors
        if (prior){
          log_accept_prob = log_accept_prob - beta_dash + beta_vec[i]
        }
        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          beta_vec[i] <- beta_new
          count_accept4 = count_accept4 + 1
          flagB0 = FALSE
        } 
        #Print
        #print(paste0('2 gamma_dash = ', gamma_dash))
        #print(paste0('log_accept_prob = ', log_accept_prob))
      } 
    }
    
    #RJMCMC Step 
    if ((beta_vec[i-1] > 0) | (gamma_vec[i-1] > 0)){ 
      print('B 0 proposal')
      beta_dash = 0
      gamma_dash = 0
      flagB0 = TRUE
      #Acceptance probability (everything cancels)
      logl_new = log_like_ss_lse_B0(data, alpha_vec[i], beta_dash, gamma_dash, flagB0)
      logl_prev = log_like_ss_lse_B0(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1], flagB0)
      log_accept_prob = logl_new - logl_prev 
      
      #Metropolis Step
      unif_var = runif(1)
      print(paste0('log_accept_prob = ', log_accept_prob))
      print(paste0('log(unif_var) = ', log(unif_var)))
      
      if(!(is.na(log_accept_prob)) && log(unif_var) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        gamma_vec[i] <- gamma_dash
        count_accept5 = count_accept5 + 1
      } else {
        flagB0 = FALSE 
      }
    } else { 
      print('B Independ. proposal')
      
      #Independence sampler - Propose from prior. If VERY lucky value is accepted to be able to jump between models. 
      beta_dash = rexp(1) 
      gamma_dash = rexp(1) + 1
      
      #Everything cancels
      logl_new = log_like_ss_lse_B0(data, alpha_vec[i], beta_dash, gamma_dash, flagB0)
      logl_prev = log_like_ss_lse_B0(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1], flagB0)
      log_accept_prob = logl_new - logl_prev 
      
      #Metropolis Step
      unif_var = runif(1)
      print(paste0('log_accept_prob = ', log_accept_prob))
      print(paste0('log(unif_var) = ', log(unif_var)))
      
      if(!(is.na(log_accept_prob)) && log(unif_var) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        gamma_vec[i] <- gamma_dash
        count_accept6 = count_accept6 + 1
        flagB0 = FALSE
        
      } } }
  
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
n_mcmc = 5000 #00 #20 #5 #0 #5 #15 #00 #5500

#### - MCMC params - ######
alphaX = 0.8 
betaX = 0.1 
gammaX = 10 
true_r0 = alphaX + betaX*gammaX
true_r0
#model_params = c(alphaX, betaX, gammaX, true_r0)

#MCMC - sigma
sigma_a = 0.4*alphaX
sigma_b = 1.0*betaX 
sigma_g = 0.85*gammaX
sigma_bg = 1.5*gammaX
sigma = c(sigma_a, sigma_b, sigma_g, sigma_bg)
#sigma_base = 0.25 #0.5

print(seed_count)
set.seed(seed_count)

#Epidemic data - Neg Bin
#sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#RUN MCMC
start_time = Sys.time()
print('Start time:')
print(start_time)
mcmc_params = rjmcmc_sse_base(sim_data, n_mcmc, sigma)
end_time = Sys.time()
print('End time:')
print(end_time)
time_elap = get_time(start_time, end_time)

#Plotting 
dist_type = 'Neg Bin,'
plot_mcmc_grid(n_mcmc, sim_data, mcmc_params, true_r0, dist_type, time_elap, seed_count)
  
#Seed
#seed_count = seed_count + 1
seed_count
