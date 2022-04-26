#SSI MODEL
setwd("~/GitHub/epidemic_modelling") 
source("epidemic_functions.R") 

#DATA
num_days = 50
#Lambda params
shape_g = 6
scale_g = 1
#SSI data params
aX = 0.8 #0.7 #0.8 1.1 #Without ss event, ~r0.
bX = 0.1 #0.15 #0.1 #0.2
cX = 10 #8#10 #8
true_r0 = aX + bX*cX
true_r0
model_params = c(aX, bX, cX, true_r0)

#MCMC params - sigma
sigma_a = 0.4*aX
sigma_b = 1.0*bX #0.1
sigma_c = 0.85*cX
sigma_bg = 1.5*cX
sigma = c(sigma_a, sigma_b, sigma_c, sigma_bg)
time_elap = 0

#MCMC
gamma_prior = FALSE
gamma_priors = c(0,0)
RJMCMCX = FALSE
alpha_transformX = FALSE

#MODEL - SSI
#Log Likelihood 
LOG_LIKE_SSI <- function(sim_data, aX, bX, cX){
  
  #Data
  n = sim_data[[1]]
  s = sim_data[[2]]
  
  #Params
  num_days = length(n)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete g) i,e Prob less than x2 - prob less than x1; the area in between 
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  logl = 0
  
  for (t in 1:num_days) { #*1 or 2
    
    lambda_t = sum((n[1:(t-1)] + cX*s[1:(t-1)])*rev(prob_infect[1:(t-1)]))
    #inner_sum_xt = 0 #
    
    logl = logl - lambda_t*(aX + bX) + n[t]*(log(aX) + log(lambda_t)) + s[t]*(log(bX) + log(lambda_t)) 
    + 2*log(1) - lfactorial(n[t]) - lfactorial(s[t])
    
  }
  
  logl
}

#APPLY
loglike = LOG_LIKE_SSI(sim_data, aX, bX, cX)
loglike

##############################
#1. MCMC
MCMC_SSI <- function(data, n_mcmc, sigma, model_params, flag_gam_prior_on_b, gam_priors_on_b,
                     x0 = 1, prior = TRUE, a_transform = FALSE,
                     DATA_AUG = TRUE) {#thinning_factor, burn_in
  
  'Returns MCMC samples of SSE model parameters (a, b, c, r0 = a + b*c) 
  w/ acceptance rates. Includes a transform, b-c transform' 
  print('MCMC SUPERSPREADER INDIVIDUALS')
  
  'Priors
  p(a) = exp(1) = rate*exp(-rate*x) = 1*exp(-1*a) = exp(-a). log(exp(-a)) = - a
  p(b) = exp(1) or p(b) = g(shape, scale), for e.g g(3, 2)
  p(c) = exp(1) + 1 = 1 + exp(-c) = exp(c - 1)'
  
  #DATA
  n_data = data[[1]]; s_data = data[[2]]
  x_data = n_data + s_data
  #x = get_total_x(data)
  print('Passed')
  
  #Initialise params
  a_vec <- vector('numeric', n_mcmc); b_vec <- vector('numeric', n_mcmc)
  c_vec <- vector('numeric', n_mcmc); r0_vec <- vector('numeric', n_mcmc)
  log_like_vec <- vector('numeric', n_mcmc)
  s_vec <- vector('numeric', n_mcmc); n_vec <- vector('numeric', n_mcmc)
  
  a_vec[1] <- model_params[1]; b_vec[1] <- model_params[2] 
  c_vec[1] <- model_params[3]; r0_vec[1] <- model_params[4];
  log_like_vec[1] <- LOG_LIKE_SSI(data, a_vec[1], b_vec[1],  c_vec[1]) 
  #Data aug params
  s_vec[1] <- x0;  n_vec[1] <- x0
  
  a = a_vec[1]; b =  b_vec[1]; 
  c = c_vec[1]; log_like = log_like_vec[1]
  s = s_vec[1]; n <- n_vec[1]
  
  #Extract params
  sigma_a = sigma[1]; sigma_b = sigma[2]
  sigma_c = sigma[3]; sigma_bg = sigma[4];
  
  #Result vectors
  count_accept1 = 0; 
  count_accept2 = 0; count_reject2 = 0;
  count_accept3 = 0; count_reject3 = 0;
  count_accept4 = 0; count_reject4 = 0;
  count_accept5 = 0; #DATA AUG
  
  #MCMC chain #
  for(i in 2:n_mcmc) {
    
    #****************************************************** 
    #a
    a_dash <- a + rnorm(1, sd = sigma_a) 
    if(a_dash < 0){
      a_dash = abs(a_dash)
    }
    
    #log a
    logl_new = LOG_LIKE_SSI(data, a_dash, b, c)
    log_accept_prob = logl_new - log_like  #+ prior1 - prior
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - a_dash + a
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      a <- a_dash
      count_accept1 = count_accept1 + 1
      log_like = logl_new
    } 
    
    #************************************************************************
    #b 
    #Only if (b > 0){ WHY? Took it out
    b_dash <- b + rnorm(1, sd = sigma_b) 
    if(b_dash < 0){
      b_dash = abs(b_dash)
    }
    #loglikelihood
    logl_new = LOG_LIKE_SSI(data, a, b_dash, c)
    log_accept_prob = logl_new - log_like #logl_prev
    
    #Priors
    if (flag_gam_prior_on_b){
      log_accept_prob = log_accept_prob + log_gamma_dist(b_dash,gam_priors_on_b) - log_gamma_dist(b,gam_priors_on_b) 
    } else {
      log_accept_prob = log_accept_prob - b_dash + b 
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      b <- b_dash
      log_like = logl_new
      count_accept2 = count_accept2 + 1
    } else {
      count_reject2 = count_reject2 + 1
    }
    
    #************************************************************************
    #c
    c_dash <- c + rnorm(1, sd = sigma_c) 
    if(c_dash < 1){
      c_dash = 2 - c_dash #Prior on c - gt 1
    }
    #Acceptance Probability
    logl_new = LOG_LIKE_SSI(data, a, b, c_dash)
    log_accept_prob = logl_new - log_like 
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - c_dash + c
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      c <- c_dash
      log_like <- logl_new
      count_accept3 = count_accept3 + 1
    } else {
      count_reject3 = count_reject3 + 1
    }
    
    #*****************************************************
    #c-b
    c_dash <- c + rnorm(1, sd = sigma_bg) #Alter sigma_bg depending on acceptance rate.
    #Acc rate too big -> Make sigma bigger. Acc rate too small -> make sigma smaller
    if(c_dash < 1){ #If less then 1
      c_dash = 2 - c_dash #abs(c_dash)
    }
    #New b
    r0 = a + b*c 
    b_new = (r0 - a)/c_dash #Proposing new c AND b. b_dash = f(R0 & c_dash)
    
    if(b_new >= 0){ #Only accept values of b > 0
      
      logl_new = LOG_LIKE_SSI(data, a, b_new, c_dash)
      log_accept_prob = logl_new - log_like 
      
      #Priors: c or Exp
      if (flag_gam_prior_on_b){
        log_accept_prob = log_accept_prob + log_gamma_dist(b_new, gam_priors_on_b) - log_gamma_dist(b, gam_priors_on_b)
      } else {  
        log_accept_prob = log_accept_prob - b_new + b
      }
      
      #Metropolis Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        b <- b_new
        c <- c_dash
        count_accept4 = count_accept4 + 1
      } else {
        count_reject4 = count_reject4 + 1
      }
    }
    
    #************************************
    #DATA AUGMENTATION
    #************************************
    if (DATA_AUG){
      
      #FOR EACH S_T
      for(t in 1:length(x_data)){
        
        #Here make a copy of 
        data_dash = data 
        #s_data called s_data_dash 
        
        #PROPOSAL for s
        u = runif(1)
        if (u < 0.5) {
          st_dash = data[[2]][t] + 1
        } else {
          st_dash = data[[2]][t] - 1 
        }
        
        
        #ACCEPTANCE PROBABILITY
        #DATA
        #s_data[t] = st_dash; n_data[t] = x_data[t] - st_dash;
        #data_aug = list(n_data, s_data)
        data_dash[[2]][t] = st_dash
        data_dash[[1]][t] =  data[[1]][t] + data[[2]][t] - st_dash 
        
        #CHECKS
        if((data_dash[[2]][t] < 0) || (data_dash[[1]][t] < 0)){
          next  
        }
        
        logl_new = LOG_LIKE_SSI(data_dash, a, b, c)
        log_accept_prob = logl_new - log_like  #+ prior1 - prior
        u_var = log(runif(1))
        
        #PRINT LOG_LIKE + ACCEPT PROB
        if (i%%100 == 0){ #%% - Modulus
          print(paste0('loglike = ', log_like))
          print(paste0('loglike new = ', logl_new))
          print(paste0('log_accept_prob new = ', log_accept_prob))
          print(paste0(' log(runif(1)) = ',  u_var))
          print('**********')
        }

        #ACCEPTANCE STEP
        if(!(is.na(log_accept_prob)) && u_var < log_accept_prob) {
          data <- data_dash
          log_like <- logl_new
          count_accept5 = count_accept5 + 1 #MAKE ACCEPT_RATE VECTOR OF LENGTH T. increase t_th count for every i
        }
      }
    }
    
    #POPPULATE MODEL PARAMETERS W/ CURRENT VALUES
    a_vec[i] <- a; b_vec[i] <- b
    c_vec[i] <- c; r0_vec[i] <- a + b*c
    log_like_vec[i] <- log_like
  }
      
  #Final stats
  accept_rate1 = 100*count_accept1/(n_mcmc-1)
  accept_rate2 = 100*count_accept2/(count_accept2 + count_reject2)
  accept_rate3 = 100*count_accept3/(count_accept3 + count_reject3)
  accept_rate4 = 100*count_accept4/(count_accept4 + count_reject4)
  accept_rate5 = 100*count_accept5/(n_mcmc-1)
  
  #Return a, acceptance rate
  return(list(a_vec, b_vec, c_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3, accept_rate4,
              count_accept2, count_accept3, count_accept4, accept_rate5, data))
}

#****************
#DATA
seed_count = 2 #seed_count = seed_count + 1 #print(paste0('i mcmc = ', i))
set.seed(seed_count)
sim_data = simulation_super_spreaders(num_days, shape_g, scale_g, aX, bX, cX)
#PLOTS
par(mfrow=c(2,1))
nt = sim_data[[1]]
plot.ts(nt, ylab = 'Daily Infections count', main = 'Super Spreaders Model - Daily Infections count')
st = sim_data[[2]]
plot.ts(st, ylab = 'Daily Infections count', main = 'Super Spreaders Model - Daily Infections count')

#Total
sim_dataX = nt + st
plot.ts(sim_dataX, ylab = 'Daily Infections count', main = 'Super Spreaders Model - Daily Infections count')

##**********
#DATA AUG
#Plot
sim_dataX = sim_data[[1]] + sim_data[[2]]
model_typeX = 'SSI'; time_elap = 0

#APPLY MCMC
n_mcmc = 100000 #0 #5000
mcmc_params = MCMC_SSI(sim_data, n_mcmc, sigma, model_params,
                       sgamma_prior, gamma_priors, DATA_AUG = FALSE)
#Plot
plot_mcmc_grid(n_mcmc, sim_dataX, mcmc_params, true_r0, time_elap, seed_count, model_type = model_typeX,
               flag_gam_prior_on_b = gamma_prior, gam_priors_on_b = gamma_priors, rjmcmc = RJMCMCX,
               mod_par_names = c('a', 'b', 'c'))

#**********************
#DATA AUG
n_mcmc = 5000
mcmc_params_da = MCMC_SSI(sim_data, n_mcmc, sigma, model_params, gamma_prior,
                       gamma_priors, DATA_AUG = TRUE)

#Plot
plot_mcmc_grid(n_mcmc, sim_dataX, mcmc_params_da, true_r0, time_elap, seed_count, model_type = model_typeX,
               flag_gam_prior_on_b = gamma_prior, gam_priors_on_b = gamma_priors, rjmcmc = RJMCMCX,
               data_aug = TRUE,
               mod_par_names = c('a', 'b', 'c'))

#Compare data
data_aug = mcmc_params[13]
n2 = data_aug[[1]]
data_aug[[2]]

