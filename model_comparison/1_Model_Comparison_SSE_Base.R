#Model Comparision
"~/GitHub/epidemic_modelling"

#setwd("~/GitHub/epidemic_modelling")
source("epidemic_functions.R")
source("helper_functions.R")

#source("~/GitHub/epidemic_modelling/epidemic_functions.R")
#source("~/GitHub/epidemic_modelling/helper_functions.R")

#Epidemic params
num_days = 50
shape_gamma = 6 #Gamma params for infectiousness curve (lambda) distribution
scale_gamma = 1 
seed_count = 1 

##############################
#1. MCMC - Loglikelihood
log_like_B0 <- function(y, alphaX) {
  
  #Params
  num_days = length(y)
  shape_gamma = 6
y  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days - 1)),                                                                                    shape = shape_gamma, scale = scale_gamma)
  logl = 0
  
  for (t in 2:num_days) {
    
    #Data
    y_t = y[t]
    lambda_t = sum(y[1:(t - 1)] * rev(prob_infect[1:(t - 1)]))
    
    if((y[t] == 0) | is.na(y[t])) {
      logl = logl - (alphaX * lambda_t)
      
    } else {
      #Add to log likelihood
      logl = logl +  y_t * log(alphaX * lambda_t) - (alphaX * lambda_t) -
        lfactorial(y_t) #- lfactorial(x[t] - y_t))
      
    }
    
  }
  
  logl
  
}

#***************************************************************#####
#GAMMA PRIOR ON BETA
log_gamma_dist <- function(param, gamma_priors){
  
  #Params
  shapeX = gamma_priors[1]
  scaleX = gamma_priors[2]
  
  log_gamma_dist = (1/lgamma(shapeX)*shapeX*log(scaleX))*(shapeX - 1)*log(param)*(-param/shapeX)
  
  log_gamma_dist
}

##################
#RJMCMC  
rjmcmc_sse_base_II <- function(data, n, sigma, model_params, gamma_prior, gamma_priors,
                            x0 = 1, prior = TRUE, alpha_transform = FALSE) {#thinning_factor, burn_in
  
  'Returns MCMC samples of SSE model parameters (alpha, beta, gamma, r0 = a + b*g) 
  w/ rjmcmc & acceptance rate. Includes alpha transform, beta-gamma transfform 
  Priors
  p(alpha) = exp(1) = rate*exp(-rate*x) = 1*exp(-1*alpha) = exp(-alpha). log(exp(-alpha)) = - alpha
  p(beta) = exp(1) or p(beta) = gamma(shape, scale), for e.g gamma(3, 2)
  p(gamma) = exp(1) + 1 = 1 + exp(-gamma) '
  
  #Initialise params
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  like_vec <- vector('numeric', n)
  
  #Alpha vecs
  #alpha_vec_i = c(); alpha_vec_ii = c(); alpha_vec_iii = c();
  
  alpha_vec[1] <- model_params[1]; beta_vec[1] <- model_params[2] #0.5 #x0;
  gamma_vec[1] <- model_params[3]; r0_vec[1] <- model_params[4];
  like_vec[1] <- log_like_ss_lse(data, alpha_vec[1], beta_vec[1],  gamma_vec[1])   
  
  #Extract params
  sigma_a = sigma[1]; sigma_b = sigma[2]
  sigma_g = sigma[3]; sigma_bg = sigma[4];
  
  #Result vectors
  count_accept1 = 0; 
  count_accept2 = 0; count_reject2 = 0;
  count_accept3 = 0; count_reject3 = 0;
  count_accept4 = 0; count_reject4 = 0;
  count_accept5 = 0; count_accept6 = 0;
  count_reject5 = 0; count_reject6 = 0;
  
  #MCMC chain
  for(i in 2:n) {
    
    #******************************************************
    #ALPHA
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma_a) 
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    
    #log alpha
    logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    log_accept_prob = logl_new - like_vec[i-1]  #+ prior1 - prior
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - alpha_dash + alpha_vec[i-1]
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
      like_vec[i] = logl_new
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
      like_vec[i] = like_vec[i-1]
    }
    #Alpha check
    #alpha_vec_i = c(alpha_vec_i,  alpha_vec[i])
    
    #if(i == 2){print(paste0('alpha_vec[2] = ')) }
    #************************************************************************
    #BETA (ONLY IF B > 0)
    if (beta_vec[i-1] > 0){ 
      
      beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma_b) 
      if(beta_dash < 0){
        beta_dash = abs(beta_dash)
      }
      #loglikelihood
      logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1])
      #logl_prev = log_like_ss_lse_B0(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
      log_accept_prob = logl_new - like_vec[i] #logl_prev
      
      #Priors
      if (gamma_prior){
        log_accept_prob = log_accept_prob + log_gamma_dist(beta_dash, gamma_priors) - log_gamma_dist(beta_vec[i-1], gamma_priors) 
      } else {
        log_accept_prob = log_accept_prob - beta_dash + beta_vec[i-1] 
      }
      
      #Metropolis Acceptance Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        like_vec[i] = logl_new
        count_accept2 = count_accept2 + 1
      } else {
        beta_vec[i] <- beta_vec[i-1]
        count_reject2 = count_reject2 + 1
      }
      
      #************************************************************************
      #GAMMA
      gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma_g) 
      if(gamma_dash < 1){
        gamma_dash = 2 - gamma_dash #abs(gamma_dash)
      }
      #Acceptance Probability
      logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash)
      #logl_prev = log_like_ss_lse_B0(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
      log_accept_prob = logl_new - like_vec[i] #logl_prev 
      #Priors
      if (prior){
        log_accept_prob = log_accept_prob - gamma_dash + gamma_vec[i-1]
      }
      #Metropolis Acceptance Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        gamma_vec[i] <- gamma_dash
        like_vec[i] <- logl_new
        count_accept3 = count_accept3 + 1
      } else {
        gamma_vec[i] <- gamma_vec[i-1]
        count_reject3 = count_reject3 + 1
      }
      
      #R0
      r0_current = alpha_vec[i] + beta_vec[i]*gamma_vec[i] #r0_vec[i] #Stays as r0 until updated again 
      
      #*****************************************************
      #GAMMA-BETA
      gamma_dash <- gamma_vec[i] + rnorm(1, sd = sigma_bg) #Alter sigma_bg depending on acceptance rate. 
      #Acc rate too big -> Make sigma bigger. Acc rate too small -> make sigma smaller
      if(gamma_dash < 1){ #If less then 1
        gamma_dash = 2 - gamma_dash #abs(gamma_dash)
      }
      #New beta  r0_vec[i]
      beta_new = (r0_current - alpha_vec[i])/gamma_dash #Proposing new Gamma AND Beta. Beta_dash = f(R0 & gamma_dash)
        
      if(beta_new >= 0){ #Only accept values of beta > 0
        
        logl_new = log_like_ss_lse(data, alpha_vec[i], beta_new, gamma_dash)
        #logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i])
        log_accept_prob = logl_new - like_vec[i] #logl_prev  
        
        #Priors
        if (gamma_prior){
          log_accept_prob = log_accept_prob + log_gamma_dist(beta_new, gamma_priors) - log_gamma_dist(beta_vec[i-1], gamma_priors) 
        } else{
          log_accept_prob = log_accept_prob - beta_new + beta_vec[i-1] 
        }
        
        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          beta_vec[i] <- beta_new
          like_vec[i] <- logl_new
          #R0
          r0_current = alpha_vec[i] + beta_vec[i]*gamma_dash
          count_accept4 = count_accept4 + 1
        } else {
          count_reject4 = count_reject4 + 1
        }
      } 
    } else {
      count_reject2 = count_reject2 + 1
      count_reject3 = count_reject3 + 1
      count_reject4 = count_reject4 + 1
    } #end of if b[i-1] > 0
    
    #************************************************************
    #RJMCMC STEP 
    #************************************************************
    
    #************************************************************
    #* M_I *#
    if ((beta_vec[i] > 0) | (gamma_vec[i] > 0)){ #Look to it 

      #print('B 0 proposal')
      beta_dash = 0
      gamma_dash = 0
      
      #alpha
      if (alpha_transform) { #ro = alpha + beta*gamma. alpha_dash_base (R0_base) = alpha_sse + beta_sse*gamma_sse (R0_SSE)
        alpha_dash = alpha_vec[i] + beta_vec[i]*gamma_vec[i] #Increase. as alpha_dash is actually the new R_0            #r0_current #r0_vec[i] #- beta_new*gamma_dash #Line added 
      } else alpha_dash = alpha_vec[i]
      
      #Alpha check
      #alpha_vec_ii = c(alpha_vec_ii, alpha_dash)
      
      #Check
      #if(is.nan(alpha_dash)){
      #print(paste0('i = ', i)); print(paste0('alpha_dash = ', alpha_dash));
      #print(paste0('beta_i = ', beta_vec[i])); print(paste0('gamma_vec[i] = ', gamma_vec[i]))
      
      #Check alpha postive==
      if (alpha_dash > 0) { #Automatically satisfied as we've increased alpha. *Remove
        
      #Acceptance probability (everything cancels)
      logl_new = log_like_B0(data, alpha_dash)
      log_accept_prob = logl_new - like_vec[i] - alpha_dash + alpha_vec[i]
      
      #logl_prev. #Multiply by 100 for example. Increase prior ratio so adequate   
      
      #Metropolis Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) { #+log(1000)
        
        beta_vec[i] <- beta_dash
        gamma_vec[i] <- gamma_dash
        alpha_vec[i] <- alpha_dash
        like_vec[i] <- logl_new
        count_accept5 = count_accept5 + 1
        #print('0s accepted')
        
      } else count_reject5 = count_reject5 + 1
  } else count_reject5 = count_reject5 + 1
      
    } else { 
      
      #************************************************************
      #* M_II 
      
      #Independence sampler - Propose from prior. If VERY lucky value is accepted to be able to jump between models. 
      beta_dash = rexp(1) 
      gamma_dash = rexp(1) + 1 
      
      #alpha
      if (alpha_transform) { #alpha_sse = ro_base (alpha_base) - beta_sse*gamma_sse
        alpha_dash = alpha_vec[i] - beta_dash*gamma_dash #(alpha_vec[i] - (beta_vec[i]*gamma_vec[i])) - beta_dash*gamma_dash #Preserves alpha, beta, gamma. Will we need the Jacobian?
      } else alpha_dash = alpha_vec[i]
      
      #print(paste0('i = ', i)); print(paste0('alpha_dash = ', alpha_dash));
      #print(paste0('beta_i = ', beta_vec[i])); print(paste0('gamma_vec[i] = ', gamma_vec[i]))
      
      #Check alpha positive==
      if (alpha_dash > 0) {
        
        #Everything cancels
        logl_new = log_like_ss_lse(data, alpha_dash, beta_dash, gamma_dash)
        log_accept_prob = logl_new - like_vec[i] - alpha_dash + alpha_vec[i]
        
        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          beta_vec[i] <- beta_dash
          gamma_vec[i] <- gamma_dash
          alpha_vec[i] <- alpha_dash
          like_vec[i] = logl_new 
          count_accept6 = count_accept6 + 1
        
        } else count_reject6 = count_reject6 + 1
      } else count_reject6 = count_reject6 + 1
      
      #Alpha check
      #alpha_vec_iii = c(alpha_vec_iii, alpha_dash)
    }
    
  }
  
  #Bayes Factor
  beta_pc0 = (length(which(beta_vec == 0)))/length(beta_vec) #Check beta_mcmc
  bayes_factor = beta_pc0/(1-beta_pc0); bayes_factor = round(bayes_factor, 6)
  
  #Final stats
  accept_rate1 = 100*count_accept1/(n-1)
  accept_rate2 = 100*count_accept2/(count_accept2 + count_reject2)
  accept_rate3 = 100*count_accept3/(count_accept3 + count_reject3)
  accept_rate4 = 100*count_accept4/(count_accept4 + count_reject4)
  #RJMCMC Steps 
  accept_rate5 = 100*count_accept5/(count_accept5 + count_reject5) #Check count_accept + count_reject = n_mcmc 
  accept_rate6 = 100*count_accept6/(count_accept6 + count_reject6)
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3, accept_rate4,
              accept_rate5, accept_rate6, count_accept5, count_accept6,
              count_reject5, count_reject6, count_accept2, count_accept3, count_accept4, beta_pc0, bayes_factor))
}


##################
#RJMCMC - Without alpha/beta transform
rjmcmc_sse_base <- function(data, n, sigma, model_params, gamma_prior, gamma_priors,
                               x0 = 1, prior = TRUE, alpha_transform = FALSE) {#thinning_factor, burn_in
  
  'Returns MCMC samples of SSE model parameters (alpha, beta, gamma, r0 = a + b*g) 
  w/ rjmcmc & acceptance rate. Includes alpha transform, beta-gamma transfform 
  Priors
  p(alpha) = exp(1) = rate*exp(-rate*x) = 1*exp(-1*alpha) = exp(-alpha). log(exp(-alpha)) = - alpha
  p(beta) = exp(1) or p(beta) = gamma(shape, scale), for e.g gamma(3, 2)
  p(gamma) = exp(1) + 1 = 1 + exp(-gamma) = exp(gamma - 1)'
  
  #Initialise params
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  like_vec <- vector('numeric', n)
  
  #Alpha vecs
  #alpha_vec_i = c(); alpha_vec_ii = c(); alpha_vec_iii = c();
  
  alpha_vec[1] <- model_params[1]; beta_vec[1] <- model_params[2] #0.5 #x0;
  gamma_vec[1] <- model_params[3]; r0_vec[1] <- model_params[4];
  like_vec[1] <- log_like_ss_lse(data, alpha_vec[1], beta_vec[1],  gamma_vec[1])   
  
  #Extract params
  sigma_a = sigma[1]; sigma_b = sigma[2]
  sigma_g = sigma[3]; sigma_bg = sigma[4];
  
  #Result vectors
  count_accept1 = 0; 
  count_accept2 = 0; count_reject2 = 0;
  count_accept3 = 0; count_reject3 = 0;
  count_accept4 = 0; count_reject4 = 0;
  count_accept5 = 0; count_accept6 = 0;
  count_reject5 = 0; count_reject6 = 0;
  
  #MCMC chain
  for(i in 2:n) {
    
    #******************************************************
    #ALPHA
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma_a) 
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    
    #log alpha
    logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    log_accept_prob = logl_new - like_vec[i-1]  #+ prior1 - prior
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - alpha_dash + alpha_vec[i-1]
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
      like_vec[i] = logl_new
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
      like_vec[i] = like_vec[i-1]
    }
    
    #************************************************************************
    #BETA (ONLY IF B > 0) ????
    if (beta_vec[i-1] > 0){ 
      
      beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma_b) 
      if(beta_dash < 0){
        beta_dash = abs(beta_dash)
      }
      #loglikelihood
      logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1])
      #logl_prev = log_like_ss_lse_B0(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
      log_accept_prob = logl_new - like_vec[i] #logl_prev
      
      #Priors
      if (gamma_prior){
        log_accept_prob = log_accept_prob + log_gamma_dist(beta_dash, gamma_priors) - log_gamma_dist(beta_vec[i-1], gamma_priors) 
      } else {
        log_accept_prob = log_accept_prob - beta_dash + beta_vec[i-1] 
      }
      
      #Metropolis Acceptance Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        like_vec[i] = logl_new
        count_accept2 = count_accept2 + 1
      } else {
        beta_vec[i] <- beta_vec[i-1]
        count_reject2 = count_reject2 + 1
      }
      
      #************************************************************************
      #GAMMA
      gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma_g) 
      if(gamma_dash < 1){
        gamma_dash = 2 - gamma_dash #abs(gamma_dash)
      }
      #Acceptance Probability
      logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash)
      #logl_prev = log_like_ss_lse_B0(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
      log_accept_prob = logl_new - like_vec[i] #logl_prev 
      #Priors
      if (prior){
        log_accept_prob = log_accept_prob - gamma_dash + gamma_vec[i-1]
      }
      #Metropolis Acceptance Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        gamma_vec[i] <- gamma_dash
        like_vec[i] <- logl_new
        count_accept3 = count_accept3 + 1
      } else {
        gamma_vec[i] <- gamma_vec[i-1]
        count_reject3 = count_reject3 + 1
      }
      
      #R0
      r0_current = alpha_vec[i] + beta_vec[i]*gamma_vec[i] #r0_vec[i] #Stays as r0 until updated again 
      
      #*****************************************************
      #GAMMA-BETA
      gamma_dash <- gamma_vec[i] + rnorm(1, sd = sigma_bg) #Alter sigma_bg depending on acceptance rate.
      #Acc rate too big -> Make sigma bigger. Acc rate too small -> make sigma smaller
      if(gamma_dash < 1){ #If less then 1
        gamma_dash = 2 - gamma_dash #abs(gamma_dash)
      }
      #New beta  r0_vec[i]
      beta_new = (r0_current - alpha_vec[i])/gamma_dash #Proposing new Gamma AND Beta. Beta_dash = f(R0 & gamma_dash)

      if(beta_new >= 0){ #Only accept values of beta > 0

        logl_new = log_like_ss_lse(data, alpha_vec[i], beta_new, gamma_dash)
        #logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i])
        log_accept_prob = logl_new - like_vec[i] #logl_prev

        #Priors
        if (gamma_prior){
          log_accept_prob = log_accept_prob + log_gamma_dist(beta_new, gamma_priors) - log_gamma_dist(beta_vec[i-1], gamma_priors)
        } else{
          log_accept_prob = log_accept_prob - beta_new + beta_vec[i-1]
        }

        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          beta_vec[i] <- beta_new
          beta_vec[i] <- beta_new
          gamma_vec[i] <- gamma_dash
          #R0
          r0_current = alpha_vec[i] + beta_vec[i]*gamma_dash
          count_accept4 = count_accept4 + 1
        } else {
          count_reject4 = count_reject4 + 1
        }
      }
    } else {
      count_reject2 = count_reject2 + 1
      count_reject3 = count_reject3 + 1
      #count_reject4 = count_reject4 + 1
    } #end of if b[i-1] > 0
    
    #************************************************************
    #RJMCMC STEP 
    #************************************************************
    
    #************************************************************
    #* M_I *#
    if ((beta_vec[i] > 0) | (gamma_vec[i] > 0)){ #Look to it 
      
      #print('B 0 proposal')
      beta_dash = 0
      gamma_dash = 0
      
      #alpha
      if (alpha_transform) { #ro = alpha + beta*gamma. alpha_dash_base (R0_base) = alpha_sse + beta_sse*gamma_sse (R0_SSE)
        alpha_dash = alpha_vec[i] + beta_vec[i]*gamma_vec[i] #Increase. as alpha_dash is actually the new R_0            #r0_current #r0_vec[i] #- beta_new*gamma_dash #Line added 
      } else alpha_dash =  alpha_vec[i] #alpha_vec[i] + rnorm(1, sd = sigma_a) #alpha_vec[i] #alpha_vec[i] + rnorm(1, sd = sigma_a) #
      
      #Check alpha positive ==
      if (alpha_dash > 0) { #Automatically satisfied as we've increased alpha. *Remove
        
        #Acceptance probability (everything cancels)
        logl_new = log_like_B0(data, alpha_dash)
        log_accept_prob = logl_new - like_vec[i] - alpha_dash + alpha_vec[i]
        
        #logl_prev. #Multiply by 100 for example. Increase prior ratio so adequate   
        
        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) { #+log(1000)
          
          beta_vec[i] <- beta_dash
          gamma_vec[i] <- gamma_dash
          alpha_vec[i] <- alpha_dash
          like_vec[i] <- logl_new
          count_accept5 = count_accept5 + 1
          
        } else count_reject5 = count_reject5 + 1
      } else count_reject5 = count_reject5 + 1
      
    } else { 
      
      #************************************************************
      #* M_II 
      
      #Independence sampler - Propose from prior. If VERY lucky value is accepted to be able to jump between models. 
      beta_dash = rexp(1) 
      gamma_dash = rexp(1) + 1 
      
      #alpha
      if (alpha_transform) { #alpha_sse = ro_base (alpha_base) - beta_sse*gamma_sse
        alpha_dash = alpha_vec[i] - beta_dash*gamma_dash #(alpha_vec[i] - (beta_vec[i]*gamma_vec[i])) - beta_dash*gamma_dash #Preserves alpha, beta, gamma. Will we need the Jacobian?
      } else alpha_dash = alpha_vec[i] #alpha_vec[i] + rnorm(1, sd = sigma_a) #alpha_vec[i] #alpha_vec[i] + rnorm(1, sd = sigma_a) #alpha_vec[i]
      
      #Check alpha positive==
      if (alpha_dash > 0) {
        
        #Everything cancels
        logl_new = log_like_ss_lse(data, alpha_dash, beta_dash, gamma_dash)
        log_accept_prob = logl_new - like_vec[i] - alpha_dash + alpha_vec[i]
        
        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          beta_vec[i] <- beta_dash
          gamma_vec[i] <- gamma_dash
          alpha_vec[i] <- alpha_dash
          like_vec[i] = logl_new 
          count_accept6 = count_accept6 + 1
          
        } else count_reject6 = count_reject6 + 1
      } else count_reject6 = count_reject6 + 1
      
      #Alpha check
      #alpha_vec_iii = c(alpha_vec_iii, alpha_dash)
    }
  }
  
  #Bayes Factor
  beta_pc0 = (length(which(beta_vec == 0)))/length(beta_vec) #Check beta_mcmc
  bayes_factor = beta_pc0/(1-beta_pc0); bayes_factor = round(bayes_factor, 6)
  
  #Final stats
  accept_rate1 = 100*count_accept1/(n-1)
  accept_rate2 = 100*count_accept2/(count_accept2 + count_reject2)
  accept_rate3 = 100*count_accept3/(count_accept3 + count_reject3)
  accept_rate4 = 100*count_accept4/(count_accept4 + count_reject4)
  #RJMCMC Steps 
  accept_rate5 = 100*count_accept5/(count_accept5 + count_reject5) #Check count_accept + count_reject = n_mcmc 
  accept_rate6 = 100*count_accept6/(count_accept6 + count_reject6)
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3, accept_rate4,
              accept_rate5, accept_rate6, count_accept5, count_accept6,
              count_reject5, count_reject6, count_accept2, count_accept3, count_accept4, beta_pc0, bayes_factor))
}


##################
#MCMC 
mcmc_sse_gprior <- function(data, n_mcmc, sigma, model_params, gamma_prior, gamma_priors,
                            x0 = 1, prior = TRUE) {
  
  'Returns MCMC samples of SSE model parameters (alpha, beta, gamma, r0 = a + b*g) 
  w/ rjmcmc & acceptance rate
  Priors
  p(alpha) = exp(1) = rate*exp(-rate*x) = 1*exp(-1*alpha) = exp(-alpha). log(exp(-alpha)) = - alpha
  p(beta) = exp(1) or p(beta) = gamma(shape, scale), for e.g gamma(3, 2)
  p(gamma) = exp(1) + 1 = 1 + exp(-gamma) '
  
  #Initialise params
  alpha_vec <- vector('numeric', n_mcmc); beta_vec <- vector('numeric', n_mcmc)
  gamma_vec <- vector('numeric', n_mcmc); r0_vec <- vector('numeric', n_mcmc)
  like_vec <- vector('numeric', n_mcmc)
  
  alpha_vec[1] <- model_params[1]; beta_vec[1] <- model_params[2] #x0;
  gamma_vec[1] <- model_params[3]; r0_vec[1] <- model_params[4];
  like_vec[1] <- log_like_ss_lse(data, alpha_vec[1], beta_vec[1],  gamma_vec[1])   
  
  #Extract params
  sigma_a = sigma[1]; sigma_b = sigma[2]
  sigma_g = sigma[3]; sigma_bg = sigma[4];
  
  #Result vectors
  count_accept1 = 0; 
  count_accept2 = 0; count_reject2 = 0;
  count_accept3 = 0; count_reject3 = 0;
  count_accept4 = 0; count_reject4 = 0;
  
  #MCMC chain
  for(i in 2:n_mcmc) {
    
    #******************************************************
    #ALPHA
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma_a) 
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    
    #log alpha
    logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    log_accept_prob = logl_new - like_vec[i-1]  #+ prior1 - prior
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - alpha_dash + alpha_vec[i-1]
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
      like_vec[i] = logl_new
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
      like_vec[i] = like_vec[i-1]
    }
    
    #************************************************************************
    #BETA 
    beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma_b) 
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }
    #loglikelihood
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1])
    #logl_prev = log_like_ss_lse_B0(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    log_accept_prob = logl_new - like_vec[i] #logl_prev
    
    #Priors
    if (gamma_prior){
      log_accept_prob = log_accept_prob + log_gamma_dist(beta_dash, gamma_priors) - log_gamma_dist(beta_vec[i-1], gamma_priors) 
    } else {
      log_accept_prob = log_accept_prob - beta_dash + beta_vec[i-1] 
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      beta_vec[i] <- beta_dash
      like_vec[i] = logl_new
      count_accept2 = count_accept2 + 1
    } else {
      beta_vec[i] <- beta_vec[i-1]
      count_reject2 = count_reject2 + 1
    }
    
    #************************************************************************
    #GAMMA
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma_g) 
    if(gamma_dash < 1){
      gamma_dash = 2 - gamma_dash #abs(gamma_dash)
    }
    #Acceptance Probability
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash)
    #logl_prev = log_like_ss_lse_B0(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
    log_accept_prob = logl_new - like_vec[i] #logl_prev 
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - gamma_dash + gamma_vec[i-1]
    }
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      gamma_vec[i] <- gamma_dash
      like_vec[i] <- logl_new
      count_accept3 = count_accept3 + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
      count_reject3 = count_reject3 + 1
    }
    
    #R0
    r0_current = alpha_vec[i] + beta_vec[i]*gamma_vec[i] #r0_vec[i] #Stays as r0 until updated again 
    
    #*****************************************************
    #GAMMA-BETA
    gamma_dash <- gamma_vec[i] + rnorm(1, sd = sigma_bg) #Alter sigma_bg depending on acceptance rate. 
    #Acc rate too big -> Make sigma bigger. Acc rate too small -> make sigma smaller
    if(gamma_dash < 1){ #If less then 1
      gamma_dash = 2 - gamma_dash #abs(gamma_dash)
    }
    #New beta  r0_vec[i]
    beta_new = (r0_current - alpha_vec[i])/gamma_dash #Proposing new Gamma AND Beta. Beta_dash = f(R0 & gamma_dash)
    
    if(beta_new >= 0){ #Only accept values of beta > 0
      
      logl_new = log_like_ss_lse(data, alpha_vec[i], beta_new, gamma_dash)
      #logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i])
      log_accept_prob = logl_new - like_vec[i] #logl_prev  
      
      #Priors
      if (gamma_prior){
        log_accept_prob = log_accept_prob + log_gamma_dist(beta_new, gamma_priors) - log_gamma_dist(beta_vec[i-1], gamma_priors) 
      } else{
        log_accept_prob = log_accept_prob - beta_new + beta_vec[i-1] 
      }
      
      #Metropolis Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        beta_vec[i] <- beta_new
        gamma_vec[i] <- gamma_dash
        like_vec[i] <- logl_new
        count_accept4 = count_accept4 + 1
      } else {
        count_reject4 = count_reject4 + 1
      }
    }
  }
  
  #Final stats
  accept_rate1 = 100*count_accept1/(n_mcmc -1)
  accept_rate2 = 100*count_accept2/(count_accept2 + count_reject2)
  accept_rate3 = 100*count_accept3/(count_accept3 + count_reject3)
  accept_rate4 = 100*count_accept4/(count_accept4 + count_reject4)
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3, accept_rate4))
}


#FIX.
# mcmc_sse_gprior  <- function(data, n, sigma, model_params, gamma_prior, gamma_priors,
#                             x0 = 1, prior = TRUE) {
#   
#   'Returns MCMC samples of SSE model parameters (alpha, beta, gamma, r0 = a + b*g) & acceptance rates'
#   print('MCMC SUPERSPREADING')
#   
#   #Initialise params
#   alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
#   gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
#   
#   alpha_vec[1] <- x0; beta_vec[1] <- x0;
#   gamma_vec[1] <- x0; r0_vec[1] <- x0;
#   
#   #Extract params
#   sigma_a = sigma[1]; sigma_b = sigma[2]
#   sigma_g = sigma[3]; sigma_bg = sigma[4];
#   
#   #Result vectors
#   count_accept1 = 0; count_accept2 = 0;
#   count_accept3 = 0; count_accept4 = 0;
#   flag_true = FALSE
#   
#   #Create folder for mcmc results 
#   #folder_mcmc = paste0(folder_results, '/mcmc')
#   #ifelse(!dir.exists(file.path(folder_mcmc)), dir.create(file.path(folder_mcmc), recursive = TRUE), FALSE)
#   
#   #MCMC chain
#   for(i in 2:n) {
#     
#     #******************************************************
#     #ALPHA
#     alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma_a) 
#     if(alpha_dash < 0){
#       alpha_dash = abs(alpha_dash)
#     }
#     #log alpha
#     logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
#     logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
#     prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
#     prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
#     log_accept_prob = logl_new - logl_prev  #+ prior1 - prior
#     #Priors
#     if (prior){
#       log_accept_prob = log_accept_prob - alpha_dash + alpha_vec[i-1]
#     }
#     
#     #Metropolis Acceptance Step
#     if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
#       alpha_vec[i] <- alpha_dash
#       count_accept1 = count_accept1 + 1
#     } else {
#       alpha_vec[i] <- alpha_vec[i-1]
#     }
#     
#     #************************************************************************
#     #BETA 
#     
#     beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma_b) 
#     #(Only if B > 0)
#     if(beta_dash < 0){
#       beta_dash = abs(beta_dash)
#     }
#     
#     #loglikelihood
#     logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1])
#     logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
#     log_accept_prob = logl_new - logl_prev
#     
#     #Priors
#     if (gamma_prior){
#       log_accept_prob = log_accept_prob + log_gamma_dist(beta_dash, gamma_priors) 
#       - log_gamma_dist(beta_vec[i-1], gamma_priors) 
#     } else {
#       log_accept_prob = log_accept_prob - beta_dash + beta_vec[i-1] 
#     }
# 
#     #Metropolis Acceptance Step
#     if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
#       beta_vec[i] <- beta_dash
#       count_accept2 = count_accept2 + 1
#     } else {
#       beta_vec[i] <- beta_vec[i-1]
#     }
#     
#     #************************************************************************
#     #GAMMA
#     gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma_g) 
#     if(gamma_dash < 1){
#       gamma_dash = 2 - gamma_dash #abs(gamma_dash)
#     }
#     #Acceptance Probability
#     logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash) 
#     logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
#     log_accept_prob = logl_new - logl_prev 
#     #Priors
#     if (prior){
#       log_accept_prob = log_accept_prob - gamma_dash + gamma_vec[i-1]
#     }
#     #Metropolis Acceptance Step
#     if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
#       gamma_vec[i] <- gamma_dash
#       count_accept3 = count_accept3 + 1
#     } else {
#       gamma_vec[i] <- gamma_vec[i-1]
#     }
#     
#     #R0
#     r0_vec[i] = alpha_vec[i] + beta_vec[i]*gamma_vec[i]
#     
#     #*****************************************************
#     #GAMMA-BETA
#     gamma_dash <- gamma_vec[i] + rnorm(1, sd = sigma_bg)#Alter sigma_bg depending on acceptance rate. 
#     #Acc rate too big -> Make sigma bigger. Acc rate too small -> make sigma smaller
#     if(gamma_dash < 1){ #If less then 1
#       gamma_dash = 2 - gamma_dash #abs(gamma_dash)
#     }
#     #New beta 
#     beta_new = (r0_vec[i] - alpha_vec[i])/gamma_dash #Proposing new Gamma AND Beta. Beta_dash = f(R0 & gamma_dash)
#     
#     if(beta_new >= 0){ #Only accept values of beta > 0
#       
#       logl_new = log_like_ss_lse(data, alpha_vec[i], beta_new, gamma_dash)
#       logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i])
#       log_accept_prob = logl_new - logl_prev  
#       #Priors
#       if (prior){
#         log_accept_prob = log_accept_prob - beta_dash + beta_vec[i]
#       }
#       #Metropolis Step
#       if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
#         beta_vec[i] <- beta_new
#         count_accept4 = count_accept4 + 1
#       } 
#     }
#     #End of if 
#     #A final 5th move -> B = 0 if not 0, and vice versa
#   }
#   
#   #Final stats
#   accept_rate1 = 100*count_accept1/n
#   accept_rate2 = 100*count_accept2/n
#   accept_rate3 = 100*count_accept3/n
#   accept_rate4 = 100*count_accept4/n
#   
#   #Return alpha, acceptance rate
#   return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
#               accept_rate1, accept_rate2, accept_rate3, accept_rate4))
# }


#Epidemic data
#SSE DATA
#sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
#model_typeX = 'SSE'
#plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#BASE DATA
#sim_data_base = simulate_branching(num_days, true_r0, shape_gamma, scale_gamma)
#model_typeX = 'Base'
#plot.ts(sim_data_base, ylab = 'Daily Infections count', main = 'Daily Infections count')

#RUN MCMC
# start_time = Sys.time()
# print('Start time:')
# print(start_time)
# #mcmc_params = rjmcmc_sse_base(sim_data, n_mcmc, sigma, model_params)
# end_time = Sys.time()
# print('End time:')
# print(end_time)
# time_elap = get_time(start_time, end_time)
# 
# #Plotting
# #plot_mcmc_grid(n_mcmc, sim_data, mcmc_params, true_r0, time_elap, seed_count, model_type = model_typeX)

#************
#GAMMA PRIORS
#RUN MCMC
# gamma_prior = TRUE
# gamma_priorsX = c(2, 2.5)
# start_time = Sys.time()
# print('Start time:')
# print(start_time)
# mcmc_params = rjmcmc_sse_base(sim_data, n_mcmc, sigma, model_params, gamma_prior, gamma_priors)
# end_time = Sys.time()
# print('End time:')
# print(end_time)
# time_elap = get_time(start_time, end_time)
# 
# # #Plotting
# plot_mcmc_grid(n_mcmc, sim_data, mcmc_params, true_r0, time_elap, seed_count, g_prior = TRUE, g_priorsX = gamma_priorsX, model_type = model_typeX)
# 
# #Bayes Factor
# bf = mcmc_params[[19]]
# bf
# 
# beta_mcmc = mcmc_params[[2]]

#Which
#b <- c(1,5,8,4,5, 6, 5, 5)
#a <- c(1, 3, 3, 4, 5, 6, 7, 111)
#ind = which(b == 5)
#a[ind]

#Test

