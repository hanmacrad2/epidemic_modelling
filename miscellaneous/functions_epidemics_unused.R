#Functions epidemic - unused 

############################
# LOG LIKELIHOODS


#Log Likelihood - log-exp-sum trick 
log_like_ss_lse_B0 <- function(x, alphaX, betaX, gammaX){
  
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
      
      #CUT FOR LOOP; REPLACE ALL YT WITH XT (ONE ITERATION) LET YT = XT. 
      for (y_t in 0:x[t]){ #Sum for all values of y_t up to x_t
        #print(paste0('y_t  = ', y_t))
        #Store inner L(x_i) term in vector position
        # inner_sum_Lx = (-(alphaX*lambda_t) - lfactorial(y_t) + y_t*log(alphaX*lambda_t) +
        #                   lgamma((x[t] - y_t) + (betaX*lambda_t)) - lgamma(betaX*lambda_t) - 
        #                   lfactorial(x[t] - y_t) - (betaX*lambda_t*log(gammaX +1)) + 
        #                   (x[t] - y_t)*log(gammaX) -(x[t] - y_t)*log(gammaX + 1))
        
        inner_sum_Lx = (-(alphaX*lambda_t) - lfactorial(y_t) + y_t*log(alphaX*lambda_t) +
                          lgamma((x[t] - y_t) + (betaX*lambda_t)) #lgamma(betaX*lambda_t) - 
                        - lfactorial(x[t] - y_t)) #- (betaX*lambda_t*log(gammaX +1)) 
        #+ (x[t] - y_t)*log(gammaX) -(x[t] - y_t)*log(gammaX + 1))
        
        #print(paste0('inner_sum_Lx  = ', inner_sum_Lx))
        
        #Check inf
        if (is.infinite(inner_sum_Lx)){
          #print(paste0('inner_sum_Lx is inf:', inner_sum_Lx))
          #print(paste0('yt value = ', y_t))
        } else {
          inner_sum_vec[y_t + 1] = inner_sum_Lx
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

######################
#MCMC

##################
#RJMCMC  
rjmcmc_sse_base <- function(data, n, sigma, model_params, x0 = 1, prior = TRUE) { #thinning_factor, burn_in
  
  'Returns mcmc samples for sse model w/ rjmcmc & acceptance rates'
  
  #Initialise params
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  like_vec <- vector('numeric', n)
  
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
  count_accept4 = 0;
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
    
    #Old
    #logl_prev = log_like_ss_lse_B0(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    #prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    #prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
    
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
      #Beta prior = exp(1) = rate*exp(-rate*x) = 1*exp(-1*beta) = exp(-beta)
      #log(exp(-beta)) = -beta
      if (prior){
        log_accept_prob = log_accept_prob - beta_dash + beta_vec[i-1] 
      }
      #Metropolis Acceptance Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        like_vec[i] = logl_new
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
        if (prior){
          log_accept_prob = log_accept_prob - beta_dash + beta_vec[i]
        }
        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          beta_vec[i] <- beta_new
          like_vec[i] <- logl_new
          count_accept4 = count_accept4 + 1
        } 
      } 
    } #end of if b[i-1] > 0
    
    #************************************************************
    #RJMCMC Step 
    if ((beta_vec[i] > 0) | (gamma_vec[i] > 0)){ #Look to it 
      
      #print('B 0 proposal')
      beta_dash = 0
      gamma_dash = 0
      alpha_dash = r0_current #r0_vec[i] #- beta_new*gamma_dash #Line added 
      #Acceptance probability (everything cancels)
      logl_new = log_like_B0(data, alpha_dash)
      log_accept_prob = logl_new - like_vec[i] #logl_prev. #Multiply by 100 for example. Increase prior ratio so adquate  
      
      #Metropolis Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        gamma_vec[i] <- gamma_dash
        alpha_vec[i] <- alpha_dash
        like_vec[i] <- logl_new
        count_accept5 = count_accept5 + 1
        #print('0s accepted')
        
      } else {
        count_reject5 = count_reject5 + 1
      }
    } else { 
      #print('B Independ. proposal')
      
      #Independence sampler - Propose from prior. If VERY lucky value is accepted to be able to jump between models. 
      beta_dash = rexp(1) 
      gamma_dash = rexp(1) + 1 
      alpha_dash = r0_current - beta_dash*gamma_dash #Preserves alpha, beta, gamma. Will we need the Jacobian? 
      
      #Check alpha postive
      if (alpha_dash > 0) {
        
        #Everything cancels
        logl_new = log_like_ss_lse(data, alpha_dash, beta_dash, gamma_dash)
        log_accept_prob = logl_new - like_vec[i]  #logl_prev #Jacobian 
        
        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          beta_vec[i] <- beta_dash
          gamma_vec[i] <- gamma_dash
          alpha_vec[i] <- alpha_dash
          like_vec[i] = logl_new 
          count_accept6 = count_accept6 + 1
          print('B independent accepted')
          
          #Print
          print(paste0('beta_dash = ', beta_dash))
          print(paste0('gamma_dash = ', gamma_dash))
          print(paste0('alpha_dash = ', alpha_dash)) 
          print(paste0('logl_new = ', logl_new)) 
          
        } else {
          count_reject6 = count_reject6 + 1
        }
      }
      
    }
    
    r0_vec[i] = r0_current
  }
  
  #Final stats
  accept_rate1 = 100*count_accept1/(n-1)
  accept_rate2 = 100*count_accept2/(n-1) #Modify - need a count_reject for beta and gamma 
  accept_rate3 = 100*count_accept3/(n-1)
  accept_rate4 = 100*count_accept4/(n-1)
  #RJMCMC Steps 
  accept_rate5 = 100*count_accept5/(count_accept5 + count_reject5) # + 1)
  accept_rate6 = 100*count_accept6/(count_accept6 + count_reject6) # + 1)
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3, accept_rate4,
              accept_rate5, accept_rate6, count_accept5, count_accept6,
              count_reject5, count_reject6))
}


##*********************************************************************
#RJMCMC function 
rjmcmc_sse_base_g_prior <- function(data, n, sigma, model_params, gamma_prior, gamma_priors,
                                    x0 = 1, prior = TRUE) { #thinning_factor, burn_in
  
  'Returns mcmc samples for sse model w/ rjmcmc & acceptance rates. Gamma(shape, scale) prior on beta.'
  
  #Initialise params
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  like_vec <- vector('numeric', n)
  
  alpha_vec[1] <- model_params[1]; beta_vec[1] <- model_params[2] #x0;
  gamma_vec[1] <- model_params[3]; r0_vec[1] <- model_params[4];
  like_vec[1] <- log_like_ss_lse(data, alpha_vec[1], beta_vec[1],  gamma_vec[1]) 
  
  #Extract params
  sigma_a = sigma[1]; sigma_b = sigma[2]
  sigma_g = sigma[3]; sigma_bg = sigma[4];
  
  #Priors
  if (gamma_prior){
    prior = FALSE
  }
  #Result vectors
  count_accept1 = 0; count_accept2 = 0;
  count_accept3 = 0; count_accept4 = 0; 
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
    
    #Old
    #logl_prev = log_like_ss_lse_B0(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    #prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    #prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
    
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
      }
      
      else if (prior){
        log_accept_prob = log_accept_prob - beta_dash + beta_vec[i-1] 
      }
      #Metropolis Acceptance Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        like_vec[i] = logl_new
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
      }
      
      #R0
      r0_current = alpha_vec[i] + beta_vec[i]*gamma_vec[i] #r0_vec[i] #Stays as r0 until updated again 
      #PRINT
      #print('****************************')
      #print(paste0('1 RO_i = ', r0_current))
      
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
        if (prior){
          log_accept_prob = log_accept_prob - beta_dash + beta_vec[i]
        }
        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          beta_vec[i] <- beta_new
          like_vec[i] <- logl_new
          count_accept4 = count_accept4 + 1
        } 
      } 
    } #end of if b[i-1] > 0
    
    #************************************************************
    #RJMCMC Step 
    if ((beta_vec[i] > 0) | (gamma_vec[i] > 0)){ #Look to it 
      
      #print('B 0 proposal')
      beta_dash = 0
      gamma_dash = 0
      alpha_dash = r0_current #r0_vec[i] #- beta_new*gamma_dash #Line added 
      #Acceptance probability (everything cancels)
      logl_new = log_like_B0(data, alpha_dash)
      #print(paste0('like_vec[i] = ', like_vec[i]))
      #print(paste0('logl_new = ', logl_new))
      log_accept_prob = logl_new - like_vec[i] #logl_prev 
      
      #Metropolis Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        gamma_vec[i] <- gamma_dash
        alpha_vec[i] <- alpha_dash
        like_vec[i] <- logl_new
        count_accept5 = count_accept5 + 1
        #print('0s accepted')
        
      } else {
        count_reject5 = count_reject5 + 1
      }
    } else { 
      #print('B Independ. proposal')
      
      #Independence sampler - Propose from prior. If VERY lucky value is accepted to be able to jump between models. 
      beta_dash = rexp(1) 
      gamma_dash = rexp(1) + 1 
      alpha_dash = r0_current - beta_dash*gamma_dash #Preserves alpha, beta, gamma. Will we need the Jacobian? 
      
      #Print
      # print(paste0('beta_dash = ', beta_dash))
      # print(paste0('gamma_dash = ', gamma_dash))
      # print(paste0('alpha_dash = ', alpha_dash))
      
      #Check alpha postive
      if (alpha_dash > 0) {
        
        #Everything cancels
        logl_new = log_like_ss_lse(data, alpha_dash, beta_dash, gamma_dash)
        log_accept_prob = logl_new - like_vec[i]  #logl_prev #Jacobian 
        
        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          beta_vec[i] <- beta_dash
          gamma_vec[i] <- gamma_dash
          alpha_vec[i] <- alpha_dash
          like_vec[i] = logl_new 
          count_accept6 = count_accept6 + 1
          #print('B independent accepted')
          
        } else {
          count_reject6 = count_reject6 + 1
        }
      }
    }
    
    r0_vec[i] = r0_current
    #Print 
    # print(paste0('beta_vec[i] = ',  beta_vec[i]))
    # print(paste0('gamma_vec[i] = ',  beta_vec[i]))
  }
  
  #Final stats
  accept_rate1 = 100*count_accept1/n 
  accept_rate2 = 100*count_accept2/n
  accept_rate3 = 100*count_accept3/n
  accept_rate4 = 100*count_accept4/n
  #RJMCMC Steps 
  accept_rate5 = 100*count_accept5/(count_accept5 + count_reject5 + 1)
  accept_rate6 = 100*count_accept6/(count_accept6 + count_reject6 + 1)
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3, accept_rate4, accept_rate5, accept_rate6))
}

#############
#PLOTS

################################################################################
# GRID PLOT
################################################################################

plot_mcmc_grid_3x4 <- function(sim_data, mcmc_params, true_r0, dist_type,
                               total_time, seed_count, prior = TRUE, flag5 = FALSE){
  
  #Plot Set up
  plot.new()
  par(mfrow=c(3,4))
  
  #Extract params
  alpha_mcmc = mcmc_params[1]; alpha_mcmc = unlist(alpha_mcmc)
  like_a = mcmc_params[8]; like_a = unlist(like_a)
  prior_a = mcmc_params[11]; prior_a = unlist(prior_a)
  
  beta_mcmc = mcmc_params[2]; beta_mcmc = unlist(beta_mcmc)
  like_b = mcmc_params[9]; like_b = unlist(like_b)
  prior_b = mcmc_params[12]; prior_b = unlist(prior_b)
  
  gamma_mcmc = mcmc_params[3]; gamma_mcmc = unlist(gamma_mcmc)
  like_g = mcmc_params[10]; like_g = unlist(like_g)
  prior_g = mcmc_params[13]; prior_g = unlist(prior_g)
  
  r0_mcmc = mcmc_params[4]; r0_mcmc = unlist(r0_mcmc)
  
  #Cumulative means + param sample limits
  #r0
  r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
  r0_lim = max(true_r0, max(r0_mcmc))
  r0_lim2 = max(true_r0, r0_mean)
  
  #alpha
  alpha_mean = cumsum(alpha_mcmc)/seq_along(alpha_mcmc)
  like_a_mean = cumsum(like_a)/seq_along(like_a)
  prior_a_mean = cumsum(prior_a)/seq_along(prior_a)
  a_lim =  max(alphaX, max(alpha_mcmc))
  a_lim2 =  max(alphaX, alpha_mean)
  
  #beta
  beta_mean = cumsum(beta_mcmc)/seq_along(beta_mcmc)
  like_b_mean = cumsum(like_b)/seq_along(like_b)
  prior_b_mean = cumsum(prior_b)/seq_along(prior_b)
  b_lim = max(betaX, max(beta_mcmc))
  b_lim2 = max(betaX, beta_mean)
  
  #gamma
  gamma_mean = cumsum(gamma_mcmc)/seq_along(gamma_mcmc)
  like_g_mean = cumsum(like_g)/seq_along(like_g)
  prior_g_mean = cumsum(prior_g)/seq_along(prior_g)
  g_lim =  max(gammaX, max(gamma_mcmc))
  g_lim2 =  max(gammaX, gamma_mean) 
  
  #***********
  #* Plots *
  
  #i.Infections
  plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste(seed_count, "Day Infts SS Evnts", dist_type, "r0 = ", true_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #ii. MCMC Trace Plots
  plot.ts(alpha_mcmc, ylab = 'alpha', ylim=c(0, a_lim),
          main = paste("SS Evnts, true a = ", alphaX, '(green), ll(red), prior(blue)'),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = alphaX, col = 'green', lwd = 2) #True = green
  lines(seq_along(like_a), like_a, col = 'red')
  lines(seq_along(prior_a), prior_a, col = 'blue')
  
  plot.ts(beta_mcmc, ylab = 'beta', ylim=c(0, b_lim),
          main = paste("MCMC SS Events, true beta = ", betaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = betaX, col = 'green', lwd = 2) #True = green
  lines(seq_along(like_b), like_b, col = 'red')
  lines(seq_along(prior_b), prior_b, col = 'blue')
  
  plot.ts(gamma_mcmc,  ylab = 'gamma', ylim=c(0,g_lim),
          main = paste("MCMC SS Events, true gamma = ", gammaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = gammaX, col = 'green', lwd = 2) #True = green
  lines(seq_along(like_g), like_g, col = 'red')
  lines(seq_along(prior_g), prior_g, col = 'blue')
  
  #plot.ts(r0_mcmc,  ylab = 'r0', main = paste("MCMC SS Events, true r0 = ", r0_true))
  
  #Title
  #text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
  #     line2user(line=2, side=3), paste('MCMC SS, True R0:', true_r0, 'Prior = ', prior), xpd=NA, cex=2, font=2)
  
  #iii. Cumulative mean plots
  #r0 Mean
  plot2 = plot(seq_along(r0_mean), r0_mean,
               ylim=c(0, r0_lim),
               xlab = 'Time', ylab = 'R0', main = paste("R0 MCMC Mean, True R0 = ", true_r0),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = true_r0, col = 'orange', lwd = 2)
  
  #alpha mean
  plot2 = plot(seq_along(alpha_mean), alpha_mean,
               ylim=c(0, a_lim),
               xlab = 'Time', ylab = 'alpha', main = paste("Alpha MCMC mean, True alpha = ",alphaX),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = alphaX, col = 'green', lwd = 2)
  lines(seq_along(like_a_mean), like_a_mean, col = 'red')
  lines(seq_along(prior_a_mean), prior_a_mean, col = 'blue')
  
  #beta mean
  plot2 = plot(seq_along(beta_mean), beta_mean,
               ylim=c(0, b_lim),
               xlab = 'Time', ylab = 'beta', main = paste("Beta MCMC mean, True beta = ",betaX),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = betaX, col = 'green', lwd = 2)
  lines(seq_along(like_b_mean), like_b_mean, col = 'red')
  lines(seq_along(prior_b_mean), prior_b_mean, col = 'blue')
  
  #gamma Mean
  plot2 = plot(seq_along(gamma_mean), gamma_mean,
               xlab = 'Time', ylab = 'gamma', main = paste("Gamma MCMC mean, True gamma = ",gammaX),
               ylim=c(0, g_lim),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = gammaX, col = 'green', lwd = 2)
  lines(seq_along(like_g_mean), like_g_mean, col = 'red')
  lines(seq_along(prior_g_mean), prior_g_mean, col = 'blue')
  
  #iv. Param Histograms (Plots 9,11,12)
  hist(r0_mcmc, freq = FALSE, breaks = 100,
       xlab = 'R0 total', #ylab = 'Density', 
       main = paste('R0 total MCMC samples. Prior = ', prior),
       xlim=c(0, r0_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = true_r0, col = 'orange', lwd = 2)
  
  #v. Beta vs gamma
  plot(beta_mcmc, gamma_mcmc,
       xlab = 'beta', ylab = 'gamma', main = 'Beta vs Gamma',
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Hist Beta 
  hist(beta_mcmc, freq = FALSE, breaks = 100,
       xlab = 'beta', #ylab = 'Density', 
       main = paste("Beta, True beta = ", betaX), 
       xlim=c(0, b_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = betaX, col = 'blue', lwd = 2)
  
  #Hist Gamma 
  hist(gamma_mcmc, freq = FALSE, breaks = 100,
       xlab = 'gamma', #ylab = 'Density', 
       main = paste("Gamma, True gamma = ", gammaX),
       xlim=c(0, g_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = gammaX, col = 'green', lwd = 2)
  
  #Final Mean Stats
  data_10_pc = 0.1*n
  a_mcmc_mean = round(mean(alpha_mcmc[n-data_10_pc:n]), 2)
  b_mcmc_mean = round(mean(beta_mcmc[n-data_10_pc:n]), 2)
  g_mcmc_mean = round(mean(gamma_mcmc[n-data_10_pc:n]), 2)
  r0_mcmc_mean = round(mean(r0_mcmc[n-data_10_pc:n]), 2)
  
  #Results
  if (flag5){
    print(paste0('flag5 = '), flag5)
    df_results <- data.frame(
      alpha = alphaX, a_mc = a_mcmc_mean,
      beta = betaX, b_mc = b_mcmc_mean,
      gamma = gammaX, g_mc = g_mcmc_mean,
      R0 = true_r0, R0_mc = r0_mcmc_mean,
      accept_rate_a = round(mcmc_params[[5]],2),
      a_rte_b = round(mcmc_params[[6]], 2),
      a_rte_g = round(mcmc_params[[7]],2),
      a_rte_rj = round(mcmc_params[[8]],2),
      tot_time = total_time) 
    
  } else {
    df_results <- data.frame(
      alpha = alphaX, a_mc = a_mcmc_mean,
      beta = betaX, b_mc = b_mcmc_mean,
      gamma = gammaX, g_mc = g_mcmc_mean,
      R0 = true_r0, R0_mc = r0_mcmc_mean,
      accept_rate_a = round(mcmc_params[[5]],2),
      a_rte_b = round(mcmc_params[[6]], 2),
      a_rte_g = round(mcmc_params[[7]],2),
      tot_time = total_time) 
  }
  
  print(df_results)
  
}

###########
#FUNCTION II

#****************
#MCMC Plots 4x4
plot_mcmc_x4_II <- function(sim_data, mcmc_params, true_r0, dist_type,
                            total_time, seed_count, prior, max_sum_val,
                            model_crit = TRUE, joint = TRUE, flag5 = TRUE){
  
  #Plot Set up
  plot.new()
  
  if (joint){
    par(mfrow=c(4,4))
  } else {
    par(mfrow=c(3,4))
  }
  
  
  #Extract params
  alpha_mcmc = mcmc_params[1]
  alpha_mcmc = unlist(alpha_mcmc)
  
  beta_mcmc = mcmc_params[2]
  beta_mcmc = unlist(beta_mcmc)
  
  gamma_mcmc = mcmc_params[3]
  gamma_mcmc = unlist(gamma_mcmc)
  
  r0_mcmc = mcmc_params[4]
  r0_mcmc = unlist(r0_mcmc)
  
  #Cumulative means + param sample limits
  #r0
  r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
  r0_lim = max(true_r0, max(r0_mcmc))
  r0_lim2 = max(true_r0, r0_mean)
  
  #alpha
  alpha_mean = cumsum(alpha_mcmc)/seq_along(alpha_mcmc)
  a_lim =  max(alphaX, max(alpha_mcmc))
  a_lim2 =  max(alphaX, alpha_mean)
  
  #beta
  beta_mean = cumsum(beta_mcmc)/seq_along(beta_mcmc)
  b_lim = max(betaX, max(beta_mcmc))
  b_lim2 = max(betaX, beta_mean)
  
  #gamma
  gamma_mean = cumsum(gamma_mcmc)/seq_along(gamma_mcmc)
  g_lim =  max(gammaX, max(gamma_mcmc))
  g_lim2 =  max(gammaX, gamma_mean) 
  
  
  #***********
  #* Plots *
  
  #i.Infections
  plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste(seed_count, "Day Infts SS Evnts", dist_type, "r0 = ", true_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #ii. MCMC Trace Plots
  plot.ts(alpha_mcmc, ylab = 'alpha', ylim=c(0, a_lim),
          main = paste("MCMC SS Events, true alpha = ", alphaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = alphaX, col = 'red', lwd = 2) #True = green
  #lines(seq_along(like_a), like_a, col = 'red')
  #lines(seq_along(prior_a), prior_a, col = 'blue')
  
  plot.ts(beta_mcmc, ylab = 'beta', ylim=c(0, b_lim),
          main = paste("MCMC SS Events, true beta = ", betaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = betaX, col = 'blue', lwd = 2) #True = green
  
  plot.ts(gamma_mcmc,  ylab = 'gamma', ylim=c(0,g_lim),
          main = paste("MCMC SS Events, true gamma = ", gammaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = gammaX, col = 'green', lwd = 2) #True = green
  
  #plot.ts(r0_mcmc,  ylab = 'r0', main = paste("MCMC SS Events, true r0 = ", r0_true))
  
  #Title
  #text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
  #     line2user(line=2, side=3), paste('MCMC SS, True R0:', true_r0, 'Prior = ', prior), xpd=NA, cex=2, font=2)
  
  #iii. Cumulative mean plots
  #r0 Mean
  plot2 = plot(seq_along(r0_mean), r0_mean,
               ylim=c(0, r0_lim),
               xlab = 'Time', ylab = 'R0', main = paste("R0 MCMC Mean, True R0 = ", true_r0),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = true_r0, col = 'orange', lwd = 2)
  
  #alpha mean
  plot2 = plot(seq_along(alpha_mean), alpha_mean,
               ylim=c(0, a_lim),
               xlab = 'Time', ylab = 'alpha', main = paste("Alpha MCMC mean, True alpha = ",alphaX),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = alphaX, col = 'red', lwd = 2)
  #lines(seq_along(like_a_mean), like_a_mean, col = 'red')
  #lines(seq_along(prior_a_mean), prior_a_mean, col = 'blue')
  
  #beta mean
  plot2 = plot(seq_along(beta_mean), beta_mean,
               ylim=c(0, b_lim),
               xlab = 'Time', ylab = 'beta', main = paste("Beta MCMC mean, True beta = ",betaX),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = betaX, col = 'blue', lwd = 2)
  
  #gamma Mean
  plot2 = plot(seq_along(gamma_mean), gamma_mean,
               xlab = 'Time', ylab = 'gamma', main = paste("Gamma MCMC mean, True gamma = ",gammaX),
               ylim=c(0, g_lim),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = gammaX, col = 'green', lwd = 2)
  
  #iv. Param Histograms (Plots 9,11,12)
  hist(r0_mcmc, freq = FALSE, breaks = 100,
       xlab = 'R0 total', #ylab = 'Density', 
       main = paste('R0 total MCMC samples. Prior = ', prior),
       xlim=c(0, r0_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = true_r0, col = 'orange', lwd = 2)
  
  # #v. Beta vs gamma
  # plot(beta_mcmc, gamma_mcmc,
  #      xlab = 'beta', ylab = 'gamma', main = 'Beta vs Gamma',
  #      cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Hist alpha 
  hist(alpha_mcmc, freq = FALSE, breaks = 100,
       xlab = 'alpha', #ylab = 'Density', 
       main = paste("alpha, True alpha = ", alphaX), 
       xlim=c(0, a_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = alphaX, col = 'red', lwd = 2)
  #Prior
  # x <- seq(from = 0, to = 20, by = 0.05)
  # exp1 = dexp(x, 1)
  # lines(seq_along(exp1), exp1, type = 'l')
  
  
  #Hist Beta 
  hist(beta_mcmc, freq = FALSE, breaks = 100,
       xlab = 'beta', #ylab = 'Density', 
       main = paste("Beta, True beta = ", betaX), 
       xlim=c(0, b_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = betaX, col = 'blue', lwd = 2)
  #Prior
  # x <- seq(from = 0, to = 20, by = 0.05)
  # exp1 = dexp(x, 1)
  # lines(seq_along(exp1), exp1, type = 'l')
  
  
  #Hist Gamma 
  hist(gamma_mcmc, freq = FALSE, breaks = 100,
       xlab = 'gamma', #ylab = 'Density', 
       main = paste("Gamma, True gamma = ", gammaX),
       xlim=c(0, g_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = gammaX, col = 'green', lwd = 2)
  #Prior
  # x <- seq(from = 0, to = 20, by = 0.05)
  # exp1 = 1 + dexp(x, 1)
  # lines(seq_along(exp1), exp1, type = 'l')
  
  #Final Mean Stats
  data_10_pc = 0.5*n #50%
  a_mcmc_mean = round(mean(alpha_mcmc[n-data_10_pc:n]), 2)
  b_mcmc_mean = round(mean(beta_mcmc[n-data_10_pc:n]), 2)
  g_mcmc_mean = round(mean(gamma_mcmc[n-data_10_pc:n]), 2)
  r0_mcmc_mean = round(mean(r0_mcmc[n-data_10_pc:n]), 2)
  
  #Joint distrbutions
  if (model_crit){
    
    model_criticism(mcmc_params, sim_data, max_sum_val)
  }
  
  #Joint distrbutions
  if (joint){
    
    # #v. r0 vs beta
    # plot(beta_mcmc, r0_mcmc,
    #      xlab = 'beta', ylab = 'R0', main = 'Beta vs R0',
    #      cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    
    #v. alpha vs beta
    plot(alpha_mcmc, beta_mcmc,
         xlab = 'alpha', ylab = 'beta', main = 'alpha vs Beta',
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    
    #v. alpha vs gamma
    plot(alpha_mcmc, gamma_mcmc,
         xlab = 'alpha', ylab = 'gamma', main = 'alpha vs gamma',
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    
    #v. beta vs gamma
    plot(beta_mcmc, gamma_mcmc,
         xlab = 'beta', ylab = 'gamma', main = 'Beta vs gamma',
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  }
  
  #Results
  if (flag5){
    df_results <- data.frame(
      alpha = alphaX, a_mc = a_mcmc_mean,
      beta = betaX, b_mc = b_mcmc_mean,
      gamma = gammaX, g_mc = g_mcmc_mean,
      R0 = true_r0, R0_mc = r0_mcmc_mean,
      accept_rate_a = round(mcmc_params[[5]],2),
      a_rte_b = round(mcmc_params[[6]], 2),
      a_rte_g = round(mcmc_params[[7]],2),
      a_rte_rj = round(mcmc_params[[8]],2),
      tot_time = total_time) 
    
  } else {
    df_results <- data.frame(
      alpha = alphaX, a_mc = a_mcmc_mean,
      beta = betaX, b_mc = b_mcmc_mean,
      gamma = gammaX, g_mc = g_mcmc_mean,
      R0 = true_r0, R0_mc = r0_mcmc_mean,
      accept_rate_a = round(mcmc_params[[5]],2),
      a_rte_b = round(mcmc_params[[6]], 2),
      a_rte_g = round(mcmc_params[[7]],2),
      tot_time = total_time) 
  }
  
}

################################################################################
# Plot MCMC Total Results
################################################################################

plot_mcmc_results <- function(sim_data, mcmc_params, true_r0, 
                              dist_type, total_time, seed_count, prior, beta_gamma_flag){
  
  #Plot Set up
  plot.new()
  par(mfrow=c(3,4))
  
  #Extract params
  alpha_mcmc = mcmc_params[1]
  alpha_mcmc = unlist(alpha_mcmc)
  
  beta_mcmc = mcmc_params[2]
  beta_mcmc = unlist(beta_mcmc)
  
  gamma_mcmc = mcmc_params[3]
  gamma_mcmc = unlist(gamma_mcmc)
  
  r0_mcmc = mcmc_params[4]
  r0_mcmc = unlist(r0_mcmc)
  
  #Cumulative means + param sample limits
  #r0
  r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
  r0_lim = max(true_r0, max(r0_mcmc))
  r0_lim2 = max(true_r0, r0_mean)
  
  #alpha
  alpha_mean = cumsum(alpha_mcmc)/seq_along(alpha_mcmc)
  a_lim =  max(alphaX, max(alpha_mcmc))
  a_lim2 =  max(alphaX, alpha_mean)
  
  #beta
  beta_mean = cumsum(beta_mcmc)/seq_along(beta_mcmc)
  b_lim = max(betaX, max(beta_mcmc))
  b_lim2 = max(betaX, beta_mean)
  
  #gamma
  gamma_mean = cumsum(gamma_mcmc)/seq_along(gamma_mcmc)
  g_lim =  max(gammaX, max(gamma_mcmc))
  g_lim2 =  max(gammaX, gamma_mean) 
  
  #Beta-gamma
  if (beta_gamma_flag){
    b_g_a_rate = mcmc_params[14]
    b_g_a_rate = round(unlist(b_g_a_rate), 2)
  } else {
    b_g_a_rate = 0
  }
  
  #***********
  #* Plots *
  
  #i.Infections
  plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste(seed_count, "Day Infts SS Evnts", dist_type, "r0 = ", true_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #ii. MCMC Trace Plots
  plot.ts(alpha_mcmc, ylab = 'alpha', ylim=c(0, a_lim),
          main = paste("MCMC SS Events, true alpha = ", alphaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = alphaX, col = 'green', lwd = 2)
  
  plot.ts(beta_mcmc, ylab = 'beta', ylim=c(0, b_lim),
          main = paste("MCMC SS Events, true beta = ", betaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = betaX, col = 'green', lwd = 2)
  
  plot.ts(gamma_mcmc,  ylab = 'gamma', ylim=c(0,g_lim),
          main = paste("MCMC SS Events, true gamma = ", gammaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = gammaX, col = 'green', lwd = 2)
  #plot.ts(r0_mcmc,  ylab = 'r0', main = paste("MCMC SS Events, true r0 = ", r0_true))
  
  #iii. Cumulative mean plots
  #r0 Mean
  plot2 = plot(seq_along(r0_mean), r0_mean,
               ylim=c(0, r0_lim),
               xlab = 'Time', ylab = 'R0', main = paste("R0 MCMC Mean, True R0 = ", true_r0),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = true_r0, col = 'orange', lwd = 2)
  
  #alpha mean
  plot2 = plot(seq_along(alpha_mean), alpha_mean,
               ylim=c(0, a_lim),
               xlab = 'Time', ylab = 'alpha', main = paste("Alpha MCMC mean, True alpha = ",alphaX),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = alphaX, col = 'red', lwd = 2)
  
  #beta mean
  plot2 = plot(seq_along(beta_mean), beta_mean,
               ylim=c(0, b_lim),
               xlab = 'Time', ylab = 'beta', main = paste("Beta MCMC mean, True beta = ",betaX),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = betaX, col = 'blue', lwd = 2)
  
  #gamma Mean
  plot2 = plot(seq_along(gamma_mean), gamma_mean,
               xlab = 'Time', ylab = 'gamma', main = paste("Gamma MCMC mean, True gamma = ",gammaX),
               ylim=c(0, g_lim),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = gammaX, col = 'green', lwd = 2)
  
  #iv. Param Histograms (Plots 9,11,12)
  hist(r0_mcmc, freq = FALSE, breaks = 100,
       xlab = 'R0 total', #ylab = 'Density', 
       main = paste('R0 total MCMC samples. Prior = ', prior),
       xlim=c(0, r0_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = true_r0, col = 'orange', lwd = 2)
  
  #v. Beta vs gamma
  plot(beta_mcmc, gamma_mcmc,
       xlab = 'beta', ylab = 'gamma', main = 'Beta vs Gamma',
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Hist Beta 
  hist(beta_mcmc, freq = FALSE, breaks = 100,
       xlab = 'beta', #ylab = 'Density', 
       main = paste("Beta, True beta = ", betaX), 
       xlim=c(0, b_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = betaX, col = 'blue', lwd = 2)
  
  #Hist Gamma 
  hist(gamma_mcmc, freq = FALSE, breaks = 100,
       xlab = 'gamma', #ylab = 'Density', 
       main = paste("Gamma, True gamma = ", gammaX),
       xlim=c(0, g_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = gammaX, col = 'green', lwd = 2)
  
  #Final Mean Stats
  data_10_pc = 0.1*n
  a_mcmc_mean = round(mean(alpha_mcmc[n-data_10_pc:n]), 2)
  b_mcmc_mean = round(mean(beta_mcmc[n-data_10_pc:n]), 2)
  g_mcmc_mean = round(mean(gamma_mcmc[n-data_10_pc:n]), 2)
  r0_mcmc_mean = round(mean(r0_mcmc[n-data_10_pc:n]), 2)
  
  #Results
  df_results <- data.frame(
    alpha = alphaX,
    a_mc = a_mcmc_mean,
    beta = betaX,
    b_mc = b_mcmc_mean,
    gamma = gammaX,
    g_mc = g_mcmc_mean,
    R0 = true_r0, 
    R0_mc = r0_mcmc_mean,
    accept_rate_a = round(mcmc_params[[5]],2),
    a_rte_b = round(mcmc_params[[6]], 2),
    a_rte_g = round(mcmc_params[[7]],2),
    b_g_a_rate = b_g_a_rate, 
    tot_time = total_time) 
  
  print(df_results)
  
}

################
#Printing
g_priorsX = c(2,2.5)
beta_prior = paste0('gamma(',  g_priorsX[1], ', ', g_priorsX[2], ')')
beta_prior

#
bf_exp1 = c()
bf_exp1 = c(bf_exp1, 2)
bf_exp1
