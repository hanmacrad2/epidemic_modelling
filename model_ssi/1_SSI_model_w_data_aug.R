#****************************************************************
#1. SSI MODEL MCMC + DATA AUGMENTATION
#****************************************************************

#SETUP
# setwd("~/GitHub/epidemic_modelling") 
# source("epidemic_functions.R") 
# source("plot_functions.R") 
# source("helper_functions.R") 
# 
# #DATA SIMULATION PARAMS
# num_days = 50
# shape_g = 6; scale_g = 1 #Infectious pressure (lambda) - gamma params
# 
# #SSI specific (*TO DO: DESIGN OF EXPERIMENTS FOR PARAM COMBINATIONS)
# aX = 0.8; bX = 0.1; cX = 10 
# true_r0 = aX + bX*cX
# true_r0
# model_params = list(m1 = aX, m2 = bX, m3 = cX, true_r0 = true_r0)
# 
# #MCMC PARAMS  
# n_mcmc = 1000
# #SIGMA
# sigma_a = 0.4*aX; sigma_b = 1.0*bX #0.1 #SHOULD SIGMA BE DEFINED MORE RIGOROUS
# sigma_c = 0.85*cX; sigma_bc = 1.5*cX
# sigma = list(sigma_a = sigma_a, sigma_b = sigma_b, #Acc rate too big -> Make sigma bigger. 
#              sigma_c = sigma_c, sigma_bc = sigma_bc) #Acc rate too small -> make sigma smaller
# time_elap = 0
# mcmc_inputs = list(n_mcmc = n_mcmc, sigma = sigma, 
#                    model_params = model_params, x0 = 1, seed_count = 3)

#****************************************************************
#1. MODEL SSI - LOG LIKELIHOOD
#****************************************************************
LOG_LIKE_SSI <- function(sim_data, aX, bX, cX){
  
  #Data
  n = sim_data[[1]]; s = sim_data[[2]]
  
  #Params
  num_days = length(n)
  shape_gamma = 6; scale_gamma = 1
  logl = 0
  
  #INFECTIOUSNESS  - Difference of 2 GAMMA distributions. Discretized 
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) -
    pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  for (t in 1:num_days) { #*1 or 2
    
    #INFECTIOUS PRESSURE - SUM OF ALL INDIVIDUALS INFECTIOUSNESS 
    lambda_t = sum((n[1:(t-1)] + cX*s[1:(t-1)])*rev(prob_infect[1:(t-1)]))
    
    #LOG-LIKELIHOOD 
    logl = logl - lambda_t*(aX + bX) + n[t]*(log(aX) + log(lambda_t)) + s[t]*(log(bX) + log(lambda_t))  + 2*log(1) - lfactorial(n[t]) - lfactorial(s[t])
  }
  
  logl
}


#************************************************************************
#1. SSI MCMC                              (W/ DATA AUGMENTATION OPTION)
#************************************************************************
MCMC_SSI <- function(data,
                     mcmc_inputs = list(n_mcmc = n_mcmc, sigma = sigma,
                                        initial_pos = list(aX = 2, bX = 0.05, cX = 15)), #THINNING FACTOR, burn_in  
                     priors_list = list(a_prior = c(1, 0), b_prior = c(10, 1/100), b_prior_exp = c(1,0),
                                        c_prior = c(10, 1), c_prior_exp = c(0.1,0)),
                     FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                       PRIOR = TRUE,
                                       B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                       FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE)) { 
  
  'Returns MCMC samples of SSI model parameters (a, b, c, r0 = a + b*c) 
  w/ acceptance rates.
  INCLUDES; DATA AUGMENTATION, B-C transform' 
  #print('MCMC SUPERSPREADER INDIVIDUALS')
  #print(data[1])
  
  'Priors
  p(a) = exp(rate) = rate*exp(-rate*x). log(r*exp(-r*x)) = log(r) - rx
      -> E.g exp(1) = 1*exp(-1*a) = exp(-a). log(exp(-a)) = - a
  p(b) = exp(1) or p(b) = g(shape, scale), for e.g g(3, 2)
  p(c) = exp(1) + 1 = 1 + exp(-c) = exp(c - 1)'
  
  
  #**********************************************
  #INITIALISE PARAMS (THINK ABOUT STARTING POINTS)
  #**********************************************
  time = length(data[[1]]); 
  n_mcmc = mcmc_inputs$n_mcmc; sigma = mcmc_inputs$sigma
  a_vec <- vector('numeric', n_mcmc); b_vec <- vector('numeric', n_mcmc)
  c_vec <- vector('numeric', n_mcmc); r0_vec <- vector('numeric', n_mcmc)
  log_like_vec <- vector('numeric', n_mcmc)
  
  #INITIALISE: MCMC[1] of MCMC VECTORS &  RUNNING PARAMS
  a_vec[1] <- mcmc_inputs$initial_pos$aX;  a = a_vec[1]
  b_vec[1] <- mcmc_inputs$initial_pos$bX; b = b_vec[1]
  c_vec[1] <-mcmc_inputs$initial_pos$cX; c = c_vec[1]
  r0_vec[1] <- 1;
  log_like_vec[1] <- LOG_LIKE_SSI(data, a_vec[1], b_vec[1], c_vec[1])
  
  #INITIALISE: ACCEPTANCE COUNTS 
  list_accept_counts = list(count_accept1 = 0, count_accept2 = 0, count_accept3 = 0,
                            count_accept4 = 0, count_accept5 = 0)
  
  mat_count_da = matrix(0, n_mcmc, time) #i x t
  non_ss = matrix(0, n_mcmc, time) #USE THINNING FACTOR
  ss = matrix(0, n_mcmc, time) #USE THINNING FACTOR
  
  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2:n_mcmc) {
    
    #****************************************************** 
    #a
    a_dash <- a + rnorm(1, sd = sigma$sigma_a) 
    if(a_dash < 0){
      a_dash = abs(a_dash)
    }
    
    #log a
    logl_new = LOG_LIKE_SSI(data, a_dash, b, c)
    log_accept_prob = logl_new - log_like  #+ prior1 - prior
    #Priors
    if (FLAGS_LIST$PRIOR){
      log_accept_prob = log_accept_prob - a_dash + a
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      a <- a_dash
      list_accept_counts$count_accept1 = list_accept_counts$count_accept1 + 1
      log_like = logl_new
    } 
    
    #************************************************************************ Only if (b > 0){ ?
    #b  
    b_dash <- b + rnorm(1, sd = sigma$sigma_b) 
    if(b_dash < 0){
      b_dash = abs(b_dash)
    }
    #loglikelihood
    logl_new = LOG_LIKE_SSI(data, a, b_dash, c)
    log_accept_prob = logl_new - log_like
    
    #Priors
    if (FLAGS_LIST$B_PRIOR_GAMMA){
      log_accept_prob = log_accept_prob +
        dgamma(b_dash, shape = priors_list$b_prior[1], scale = priors_list$b_prior[2], log = TRUE) -
        dgamma(b, shape = priors_list$b_prior[1], scale = priors_list$b_prior[2], log = TRUE)
    } else {
      log_accept_prob = log_accept_prob - b_dash + b 
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      b <- b_dash
      log_like = logl_new
      list_accept_counts$count_accept2 = list_accept_counts$count_accept2 + 1
    } 
    
    #************************************************************************
    #c
    c_dash <- c + rnorm(1, sd = sigma$sigma_c) 
    if(c_dash < 1){
      c_dash = 2 - c_dash #Prior on c: > 1
    }
    #Acceptance Probability
    logl_new = LOG_LIKE_SSI(data, a, b, c_dash)
    log_accept_prob = logl_new - log_like 
    
    #Priors
    if(FLAGS_LIST$C_PRIOR_GAMMA){
      log_accept_prob = log_accept_prob + dgamma(c_dash, shape = priors_list$c_prior[1], scale = priors_list$c_prior[1], log = TRUE) -
        dgamma(c, shape = priors_list$c_prior[1], scale = priors_list$c_prior[2], log = TRUE)
    } else {
      log_accept_prob = log_accept_prob - priors_list$c_prior_exp[1]*c_dash + priors_list$c_prior_exp[1]*c
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      c <- c_dash
      log_like <- logl_new
      list_accept_counts$count_accept3 = list_accept_counts$count_accept3 + 1
    }
    
    #*****************************************************
    #B-C TRANSFORM
    if(FLAGS_LIST$BC_TRANSFORM){
      
      c_dash <- c + rnorm(1, sd = sigma$sigma_bc)
      #Prior > 1
      if(c_dash < 1){
        c_dash = 2 - c_dash
      }
      #New b
      b_transform = ((a + b*c) - a)/c_dash #b = (r0 - a)c
      
      if(b_transform >= 0){ #Only accept values of b > 0
        
        logl_new = LOG_LIKE_SSI(data, a, b_transform, c_dash)
        log_accept_prob = logl_new - log_like
        
        #PRIORS
        #b prior
        if (FLAGS_LIST$B_PRIOR_GAMMA) {
          tot_b_prior = dgamma(b_transform, shape = priors_list$b_prior[1], scale = priors_list$b_prior[2], log = TRUE) -
            dgamma(b, shape = priors_list$b_prior[1], scale = priors_list$b_prior[2], log = TRUE)
        } else { 
          tot_b_prior = - b_transform + b #exp(1) piror
        }
        
        #c prior
        if (FLAGS_LIST$C_PRIOR_GAMMA) {
          tot_c_prior = dgamma(c_dash, shape = priors_list$c_prior[1], scale = priors_list$c_prior[2], log = TRUE) -
            dgamma(c, shape = priors_list$c_prior[1], scale = priors_list$c_prior[2], log = TRUE)
        } else { 
          tot_c_prior = priors_list$c_prior_exp[1]*c_dash + priors_list$c_prior_exp[1]*c 
        }
        
        #LOG ACCEPT PROB
        log_accept_prob = log_accept_prob + tot_b_prior + tot_c_prior 
        
        #Metropolis Step
        if (!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          b <- b_transform
          c <- c_dash
          log_like <- logl_new
          list_accept_counts$count_accept4 = list_accept_counts$count_accept4 + 1
        }
      }
    }
    
    #************************************
    #DATA AUGMENTATION 
    #************************************
    if (FLAGS_LIST$DATA_AUG){
      
      if(i == 2){
        print('DATA AUG CHECK')
      }
      #FOR EACH S_T
      for(t in 1:time){
        
        #Copy of data (or update as necessary)
        data_dash = data
        
        #STOCHASTIC PROPOSAL for s
        if (runif(1) < 0.5) {
          st_dash = data[[2]][t] + 1
        } else {
          st_dash = data[[2]][t] - 1 
        }
        
        #ACCEPTANCE PROBABILITY
        data_dash[[2]][t] = st_dash #s_t = st_dash 
        data_dash[[1]][t] =  data[[1]][t] + data[[2]][t] - st_dash #n_t = x_t - s_t
        
        #CRITERIA FOR S_T & N_T  
        if((data_dash[[2]][t] < 0) || (data_dash[[1]][t] < 0)){
          #print(paste0(i, t, 'WARNING'))
          #Store
          non_ss[i, t] = data[[1]][t]
          ss[i, t] = data[[2]][t]
          next  
        } 
        #print(paste0(i, t, 'WARNING - SHOULDNT MATCH'))
        
        logl_new = LOG_LIKE_SSI(data_dash, a, b, c)
        log_accept_prob = logl_new - log_like  
        
        #PRINT LOG_LIKE + ACCEPT PROB
        # if (t == 1) {
        #   if ((i == 2) || (i%%50 == 0)){ # i%%100 == 0 #Modulus 
        #     print(paste0('i: ', i,  ', t: ', t))
        #     nt_dash =  data[[1]][t] + data[[2]][t] - st_dash
        #     print(paste0('st: ', data[[2]][t],  ', st_dash: ', st_dash))
        #     print(paste0('nt: ', data[[1]][t],  ', nt_dash: ', nt_dash))
        #     print(paste0('loglike = ', log_like))
        #     print(paste0('loglike new = ', logl_new))
        #     print(paste0('log_accept_prob new = ', log_accept_prob))
        #     #print(paste0(' log(runif(1)) = ', u_var))
        #     print('**********')
        #   }
        # }
          
        #METROPOLIS ACCEPTANCE STEP
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          
          #PRINT LOG_LIKE + ACCEPT PROB
          # if (t == 1) {
          #   #if ((i == 1) || (i%%50 == 0)){ # i%%100 == 0 #Modulus 
          #   print(paste0('i: ', i,  ', t: ', t))
          #   nt_dash =  data[[1]][t] + data[[2]][t] - st_dash
          #   print(paste0('st: ', data[[2]][t],  ', st_dash: ', st_dash))
          #   print(paste0('nt: ', data[[1]][t],  ', nt_dash: ', nt_dash))
          #   print(paste0('loglike = ', log_like))
          #   print(paste0('loglike new = ', logl_new))
          #   print(paste0('log_accept_prob new = ', log_accept_prob))
          #   #print(paste0(' log(runif(1)) = ', u_var))
          #   print('**********')
          #   #}
          # }
          
          #ACCEPT
          data <- data_dash
          log_like <- logl_new
          mat_count_da[i, t] = mat_count_da[i, t] + 1
          list_accept_counts$count_accept5 = list_accept_counts$count_accept5 + 1
        }
        
        #Store
        non_ss[i, t] = data[[1]][t] #TAKE MEAN ACROSS MCMC DIMENSION (PLOT 0 > 50)
        ss[i, t] = data[[2]][t]
      }
    }
    
    #Loglikelihood Check (Passing - no error)
    if (log_like!=LOG_LIKE_SSI(data, a, b, c)) print(paste0('ERROR! logl diff = ', log_like - LOG_LIKE_SSI(data, a, b, c)))
    
    #POPPULATE MODEL PARAMETERS W/ CURRENT VALUES
    a_vec[i] <- a; b_vec[i] <- b
    c_vec[i] <- c; r0_vec[i] <- a + b*c
    log_like_vec[i] <- log_like
    }
  
  #NON-SS & SS OUTPUT
  
  #Final stats
  accept_rate1 = 100*list_accept_counts$count_accept1/(n_mcmc-1)
  accept_rate2 = 100*list_accept_counts$count_accept2/(n_mcmc-1) #(list_accept_counts$count_accept2 + list_reject_counts$count_accept2)
  accept_rate3 = 100*list_accept_counts$count_accept3/(n_mcmc-1) 
  accept_rate4 = 100*list_accept_counts$count_accept4/(n_mcmc-1)
  accept_rate5 = 100*list_accept_counts$count_accept5/((n_mcmc-1)*time) #i x t
  
  #Acceptance rates 
  list_accept_rates = list(accept_rate1 = accept_rate1,
                           accept_rate2 = accept_rate2, accept_rate3 = accept_rate3,
                           accept_rate4 = accept_rate4, accept_rate5 = accept_rate5)
  print(list_accept_rates)
  #Return a, acceptance rate
  return(list(a_vec = a_vec, b_vec = b_vec, c_vec = c_vec, r0_vec = r0_vec,
              list_accept_rates = list_accept_rates, 
              data = data, mat_count_da = mat_count_da, #13, 14
              non_ss = non_ss, ss = ss, #15, 16
              non_ss_mean = colMeans(non_ss),
              ss_mean = colMeans(ss))) #16 
}

