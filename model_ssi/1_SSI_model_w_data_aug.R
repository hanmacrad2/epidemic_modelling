#****************************************************************
#1. SSI MODEL MCMC + DATA AUGMENTATION
#****************************************************************

#SETUP
setwd("~/GitHub/epidemic_modelling") 
source("epidemic_functions.R") 
source("helper_functions.R") 

#DATA SIMULATION PARAMS
num_days = 50
shape_g = 6; scale_g = 1 #Infectious pressure (lambda) - gamma params

#SSI specific (*TO DO: DESIGN OF EXPERIMENTS FOR PARAM COMBINATIONS)
aX = 0.8; bX = 0.1; cX = 10 
true_r0 = aX + bX*cX
true_r0
model_params = list(m1 = aX, m2 = bX, m3 = cX, true_r0 = true_r0)

#MCMC PARAMS  
#SIGMA
sigma_a = 0.4*aX; sigma_b = 1.0*bX #0.1 #SHOULD SIGMA BE DEFINED MORE RIGOROUS
sigma_c = 0.85*cX; sigma_bc = 1.5*cX
sigma = list(sigma_a = sigma_a, sigma_b = sigma_b, #Acc rate too big -> Make sigma bigger. 
             sigma_c = sigma_c, sigma_bc = sigma_bc) #Acc rate too small -> make sigma smaller
time_elap = 0
mcmc_inputs = list(n_mcmc = n_mcmc, sigma = sigma, 
                   model_params = model_params, x0 = 1)

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

#APPLY
#loglike = LOG_LIKE_SSI(sim_data, aX, bX, cX)
#loglike

#************************************************************************
#1. SSI MCMC                              (W/ DATA AUGMENTATION OPTION)
#************************************************************************
MCMC_SSI <- function(data,
                     mcmc_inputs = list(n_mcmc = n_mcmc, sigma = sigma, 
                                        model_params = model_params, x0 = 1), #THINNING FACTOR, burn_in  
                     priors_list = list(a_prior = c(1, 0), b_prior = c(10, 1/100), b_prior_exp = c(1,0),
                                        c_prior = c(10, 1), c_prior_exp = c(0.1,0)),
                     FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                       PRIOR = TRUE,
                                       B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                       FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE)) { 
  
  'Returns MCMC samples of SSI model parameters (a, b, c, r0 = a + b*c) 
  w/ acceptance rates.
  INCLUDES; DATA AUGMENTATION, B-C transform' 
  print('MCMC SUPERSPREADER INDIVIDUALS')
  
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
  
  #INITIALISE: MCMC[1] of MCMC VECTORS
  a_vec[1] <- model_params$m1; b_vec[1] <-model_params$m2
  c_vec[1] <- model_params$m3; r0_vec[1] <- model_params$true_r0;
  log_like_vec[1] <- LOG_LIKE_SSI(data, a_vec[1], b_vec[1], c_vec[1])
  
  #INITIALISE: RUNNING PARAMS
  a = model_params$m1; b =  model_params$m2; 
  c = model_params$m3; log_like = log_like_vec[1]
  
  #INITIALISE: ACCEPTANCE COUNTS 
  list_accept_counts = list(count_accept1 = 0, count_accept2 = 0, count_accept3 = 0,
                            count_accept4 = 0, count_accept5 = 0)
  list_reject_counts = list(count_reject2 = 0, count_reject3 = 0,
                            count_reject4 = 0)
  
  mat_count_da = matrix(0, n_mcmc, time) #i x t
  n_non_super_spreaders = matrix(0, n_mcmc, time) #USE THINNING FACTOR
  s_super_spreaders = matrix(0, n_mcmc, time) #USE THINNING FACTOR
  
  #******************************
  #MCMC CHAIN
  #******************************
  for(i in 2: n_mcmc) {
    
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
    
    #************************************************************************
    #b 
    #Only if (b > 0){ ?
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
      
    } else {
      list_reject_counts$n_reject2 = list_reject_counts$n_reject2 + 1
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
      list_accept_counts$count_accept3 =  list_accept_counts$count_accept3 + 1
    } else {
      list_reject_counts$n_reject3 = list_reject_counts$n_reject3 + 1
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
            dgamma(b, shape = priors_list$c_prior[1], scale = priors_list$c_prior[2], log = TRUE)
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
        } else {
          list_reject_counts$n_reject4 = list_reject_counts$n_reject4 + 1
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
        data_dash = data#*****

        #MAKE STARTING POINTS EXTREME; BEST WAY TO DO??
        if(FLAGS_LIST$FLAG_NS_DATA_AUG){ 
          data[[1]][t] = data[[1]][t] + data[[2]][t]
          data[[2]][t] = 0
        } else if (FLAGS_LIST$FLAG_SS_DATA_AUG){
          data[[2]][t] = data[[1]][t] + data[[2]][t]
          data[[1]][t] = 0
        }
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
          
          #Store
          n_non_super_spreaders[i, t] = data[[1]][t]
          s_super_spreaders[i, t] = data[[2]][t]
          next  
        } 
        
        logl_new = LOG_LIKE_SSI(data_dash, a, b, c)
        log_accept_prob = logl_new - log_like  
        
        #METROPOLIS ACCEPTANCE STEP
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          
          #ACCEPT
          data <- data_dash
          log_like <- logl_new
          mat_count_da[i, t] = mat_count_da[i, t] + 1
          list_accept_counts$count_accept5 = list_accept_counts$count_accept5 + 1
        }
        
        #Store
        n_non_super_spreaders[i, t] = data[[1]][t]
        s_super_spreaders[i, t] = data[[2]][t]
      }
    }
    
    #Loglikelihood Check (Passing - no error)
    if (log_like!=LOG_LIKE_SSI(data, a, b, c)) print(paste0('ERROR! logl diff = ', log_like - LOG_LIKE_SSI(data, a, b, c)))
    
    #POPPULATE MODEL PARAMETERS W/ CURRENT VALUES
    a_vec[i] <- a; b_vec[i] <- b
    c_vec[i] <- c; r0_vec[i] <- a + b*c
    log_like_vec[i] <- log_like
  }
  
  #Final stats
  accept_rate1 = 100*list_accept_counts$count_accept1/(n_mcmc-1)
  accept_rate2 = 100*list_accept_counts$count_accept2/(list_accept_counts$count_accept2 + list_reject_counts$count_accept2)
  accept_rate3 = 100*list_accept_counts$count_accept3/(list_accept_counts$count_accept3 + list_reject_counts$count_accept3)
  accept_rate4 = 100*list_accept_counts$count_accept4/(list_accept_counts$count_accept4 + list_reject_counts$count_accept4)
  accept_rate5 = 100*list_accept_counts$count_accept5/((n_mcmc-1)*time) #i x t
  
  #Acceptance rates 
  list_accept_rates = list(accept_rate1 = accept_rate1,
                           accept_rate2 = accept_rate2, accept_rate3 = accept_rate3,
                           accept_rate4 = accept_rate4, accept_rate5 = accept_rate5)
  #Return a, acceptance rate
  return(list(a_vec = a_vec, b_vec = b_vec, c_vec = c_vec, r0_vec = r0_vec,
              list_accept_rates = list_accept_rates, 
              data = data, mat_count_da = mat_count_da, #13, 14
              n_non_super_spreaders = n_non_super_spreaders, #15
              s_super_spreaders = s_super_spreaders)) #16 
}


#****************************************************************
#DATASET - GENERATED USING SIMULATION FUNCTIONS
#****************************************************************
seed_count = 3 #seed_count = seed_count + 1 #print(paste0('i mcmc = ', i))
set.seed(seed_count)
sim_data = simulation_super_spreaders(num_days, shape_g, scale_g, aX, bX, cX)

#PLOTS
par(mfrow=c(2,1))
non_ss = sim_data[[1]]
plot.ts(non_ss, ylab = 'Daily Infections count', main = 'Non Super-Spreaders' )
ss = sim_data[[2]]
plot.ts(ss, ylab = 'Daily Infections count', main = 'Super-Spreaders')

#Total
sim_dataX = non_ss + ss
plot.ts(sim_dataX, ylab = 'Daily Infections count', main = 'Total - Super Spreaders Model, Daily Infections count')

#****************************************************************
# APPLY MCMC SSI MODEL   
#***************************************************************
n_mcmc = 1000 #100000 #100000 
mcmc_inputs = list(n_mcmc = n_mcmc, sigma = sigma, 
                   model_params = model_params, x0 = 1)

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))

mcmc_params_da1 = MCMC_SSI(sim_data, mcmc_inputs = mcmc_inputs,
                           FLAGS_LIST = list(DATA_AUG = FALSE, BC_TRANSFORM = TRUE,
                                             PRIOR = TRUE,
                                             B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                             FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)

#PLOT RESULTS
model_typeX = 'SSI'; 
plot_mcmc_grid(n_mcmc, sim_dataX, mcmc_params_da1, true_r0, time_elap, seed_count,
               model_params,
               model_type = model_typeX,
               FLAG_G_PRIOR_B = TRUE, gam_priors_on_b = c(10, 1/100),
               rjmcmc = RJMCMCX, data_aug = TRUE,
               mod_par_names = c('a', 'b', 'c'))

#****************************************************************
# APPLY MCMC SSI MODEL + DATA AUGMENTATION  
#***************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_params_da2 = MCMC_SSI(sim_data, mcmc_inputs = mcmc_inputs)

end_time = Sys.time()
time_elap = get_time(start_time, end_time)

#PLOT RESULTS
model_typeX = 'SSI'; 
plot_mcmc_grid(n_mcmc, sim_dataX, mcmc_params_da2, true_r0, time_elap, seed_count,
               model_params,
               model_type = model_typeX,
               FLAG_G_PRIOR_B = TRUE, gam_priors_on_b = c(10, 1/100),
               rjmcmc = RJMCMCX, data_aug = TRUE,
               mod_par_names = c('a', 'b', 'c'))

#****************************************************************
# APPLY MCMC SSI MODEL + NON-SS EXTREME CASE
#***************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))

mcmc_params_da3 = MCMC_SSI(sim_data, mcmc_inputs = mcmc_inputs,
                           FLAGS_LIST = list(DATA_AUG = FALSE, BC_TRANSFORM = TRUE,
                                             PRIOR = TRUE,
                                             B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                             FLAG_NS_DATA_AUG = TRUE, FLAG_SS_DATA_AUG = FALSE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)

#PLOT RESULTS
model_typeX = 'SSI'; 
plot_mcmc_grid(n_mcmc, sim_dataX, mcmc_params_da3, true_r0, time_elap, seed_count,
               model_params,
               model_type = model_typeX,
               FLAG_G_PRIOR_B = TRUE, gam_priors_on_b = c(10, 1/100),
               rjmcmc = RJMCMCX, data_aug = TRUE,
               mod_par_names = c('a', 'b', 'c'))

#****************************************************************
# APPLY MCMC SSI MODEL + SS EXTREME CASE
#***************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))

mcmc_params_da4 = MCMC_SSI(sim_data, mcmc_inputs = mcmc_inputs,
                           FLAGS_LIST = list(DATA_AUG = FALSE, BC_TRANSFORM = TRUE,
                                             PRIOR = TRUE,
                                             B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                             FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = TRUE))

end_time = Sys.time()
time_elap = get_time(start_time, end_time)

#PLOT RESULTS
model_typeX = 'SSI'; 
plot_mcmc_grid(n_mcmc, sim_dataX, mcmc_params_da4, true_r0, time_elap, seed_count,
               model_params,
               model_type = model_typeX,
               FLAG_G_PRIOR_B = TRUE, gam_priors_on_b = c(10, 1/100),
               rjmcmc = RJMCMCX, data_aug = TRUE,
               mod_par_names = c('a', 'b', 'c'))

#**************************************#**************************************#**************************************
#2. PRIORS; GAMMA PRIORS ON B & C

#MCMC RUNS
#TO DO
#INCLUDE GAMMA PRIOR ON BETA (HUMP SHAPED) - DONE :D
#RUN WITH DIFFERENT STARTING VALUES :) **TO DO
# - +Gelmans
#INCLUDE PRIORS ON PLOT (HIST PLOT) **TO DO


#SSI; DIFFERENT STARTING POINTS

