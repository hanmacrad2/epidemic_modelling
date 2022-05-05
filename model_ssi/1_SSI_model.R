#****************************************************************
#1. SSI MODEL MCMC
#****************************************************************

#SETUP
setwd("~/GitHub/epidemic_modelling") 
source("epidemic_functions.R") 

#SIM DATA
num_days = 50
shape_g = 6; scale_g = 1 #Lambda params

#SSI specific (*DESIGN OF EXPERIMENTS FOR PARAM COMBINATIONS)
aX = 0.8; bX = 0.1; cX = 10 
true_r0 = aX + bX*cX
true_r0
model_params = c(aX, bX, cX, true_r0)

#MCMC params 
sigma_a = 0.4*aX; sigma_b = 1.0*bX #0.1
sigma_c = 0.85*cX; sigma_bc = 1.5*cX
sigma = c(sigma_a, sigma_b, sigma_c, sigma_bc)
gamma_prior = FALSE; gamma_priors = c(0,0)
RJMCMCX = FALSE; alpha_transformX = FALSE
time_elap = 0

#****************************************************************
#1. MODEL SSI - LOG LIKELIHOOD
#****************************************************************
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
#loglike = LOG_LIKE_SSI(sim_data, aX, bX, cX)

#****************************************************************
#1. SSI MCMC
#****************************************************************
MCMC_SSI <- function(data, n_mcmc, sigma, model_params,
                     flag_gam_prior_on_b, gam_priors_on_b, x0 = 1, 
                     prior = TRUE, DATA_AUG = TRUE, BC_TRANSFORM = FALSE) { #THINNING FACTOR, burn_in
  
  'Returns MCMC samples of SSI model parameters (a, b, c, r0 = a + b*c) 
  w/ acceptance rates.
  INCLUDES; DATA AUGMENTATION, B-C transform' 
  print('MCMC SUPERSPREADER INDIVIDUALS')
  
  'Priors
  p(a) = exp(rate) = rate*exp(-rate*x) --> exp(1) = 1*exp(-1*a) = exp(-a). log(exp(-a)) = - a
  p(b) = exp(1) or p(b) = g(shape, scale), for e.g g(3, 2)
  p(c) = exp(1) + 1 = 1 + exp(-c) = exp(c - 1)'
  
  #EXTRACT DATA + PARAMS
  time = length(data[[1]]);
  #SIGMA. Alter depending on acceptance rate
  sigma_a = sigma[1]; sigma_b = sigma[2]  #Acc rate too big -> Make sigma bigger. 
  sigma_c = sigma[3]; sigma_bc = sigma[4]; #Acc rate too small -> make sigma smaller
  
  #**********************************************
  #INITIALISE PARAMS
  #**********************************************
  a_vec <- vector('numeric', n_mcmc); b_vec <- vector('numeric', n_mcmc)
  c_vec <- vector('numeric', n_mcmc); r0_vec <- vector('numeric', n_mcmc)
  log_like_vec <- vector('numeric', n_mcmc)
  
  #Initialise 1st elements of mcmc vectors 
  a_vec[1] <- model_params[1]; b_vec[1] <- model_params[2] 
  c_vec[1] <- model_params[3]; r0_vec[1] <- model_params[4];
  log_like_vec[1] <- LOG_LIKE_SSI(data, a_vec[1], b_vec[1],  c_vec[1]) 
  
  #Initialise running parameters 
  a = a_vec[1]; b =  b_vec[1]; 
  c = c_vec[1]; log_like = log_like_vec[1]
  
  #Initialise Result vectors
  count_accept1 = 0; 
  count_accept2 = 0; count_reject2 = 0;
  count_accept3 = 0; count_reject3 = 0;
  count_accept4 = 0; count_reject4 = 0;
  count_accept5 = 0 
  mat_count_da = matrix(0, n_mcmc, time) #i x t
  b_count = 0; index_b_count = 2;
  n_non_super_spreaders = matrix(0, n_mcmc, time) #USE THINNING FACTOR
  s_super_spreaders = matrix(0, n_mcmc, time) #USE THINNING FACTOR
  
  #******************************
  #MCMC CHAIN
  #******************************
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
    log_accept_prob = logl_new - log_like
    
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
      b_count = b_count + 1
      
      #PRINT Beta values 
      # if (b_count < 10){ # i%%100 == 0 #Modulus 
      #   print(paste0('count i: ', i))
      #   print(paste0('b: ', b))
      #   print(paste0('loglike = ', log_like))
      #   print(paste0('log_accept_prob new = ', log_accept_prob))
      #   print('**********')
      #   index_b_count = i + 1
      # }
      
    } else {
      count_reject2 = count_reject2 + 1
    }
    
    #PRINT Beta + 1 value
    if (i == index_b_count){ # i%%100 == 0 #Modulus 
      print(paste0('count i: ', i))
      print(paste0('b: ', b))
      print(paste0('loglike = ', log_like))
      print(paste0('log_accept_prob new = ', log_accept_prob))
      print('**********')
    }
    
    #************************************************************************
    #c
    c_dash <- c + rnorm(1, sd = sigma_c) 
    if(c_dash < 1){
      c_dash = 2 - c_dash #Prior on c: > 1
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
    #b-c
    if(BC_TRANSFORM){
      c_dash <- c + rnorm(1, sd = sigma_bc)
      if(c_dash < 1){
        c_dash = 2 - c_dash
      }
      #New b
      b_new = ((a + b*c) - a)/c_dash
      
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
      
    }
    
    #************************************
    #DATA AUGMENTATION
    #************************************
    if (DATA_AUG){
      
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
          
          #Store
          n_non_super_spreaders[i, t] = data[[1]][t]
          s_super_spreaders[i, t] = data[[2]][t]
          next  
        }
        
        logl_new = LOG_LIKE_SSI(data_dash, a, b, c)
        log_accept_prob = logl_new - log_like  
        u_var = log(runif(1))

        #METROPOLIS ACCEPTANCE STEP
        if(!(is.na(log_accept_prob)) && u_var < log_accept_prob) {
          
          #PRINT LOG_LIKE + ACCEPT PROB
          if (t == 1) {
            if ((i == 1) || (i%%500 == 0)){ # i%%100 == 0 #Modulus 
              print(paste0('i: ', i,  ', t: ', t))
              nt_dash =  data[[1]][t] + data[[2]][t] - st_dash
              print(paste0('st: ', data[[2]][t],  ', st_dash: ', st_dash))
              print(paste0('nt: ', data[[1]][t],  ', nt_dash: ', nt_dash))
              print(paste0('loglike = ', log_like))
              print(paste0('loglike new = ', logl_new))
              print(paste0('log_accept_prob new = ', log_accept_prob))
              print(paste0(' log(runif(1)) = ', u_var))
              print('**********')
            }
          }
          
          #ACCEPT
          data <- data_dash
          log_like <- logl_new
          mat_count_da[i, t] = mat_count_da[i, t] + 1
          count_accept5 = count_accept5 + 1 #MAKE ACCEPT_RATE VECTOR OF LENGTH T. increase t_th count for every i=
        }
        
        #Store
        n_non_super_spreaders[i, t] = data[[1]][t]
        s_super_spreaders[i, t] = data[[2]][t]
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
  accept_rate5 = 100*count_accept5/((n_mcmc-1)*time) #i x t
  
  #Return a, acceptance rate
  return(list(a_vec, b_vec, c_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3, accept_rate4,
              count_accept2, count_accept3, count_accept4, accept_rate5, 
              data, mat_count_da, n_non_super_spreaders, #13, 14, 15, 16
              s_super_spreaders))
}

#****************************************************************
#DATA
#****************************************************************
seed_count = 2 
set.seed(seed_count)
sim_data = simulation_super_spreaders(num_days, shape_g, scale_g, aX, bX, cX)

#PLOTS
par(mfrow=c(2,1))
nt = sim_data[[1]]
plot.ts(nt, ylab = 'Daily Infections count', main = 'Non Super-Spreaders' )
st = sim_data[[2]]
plot.ts(st, ylab = 'Daily Infections count', main = 'Super-Spreaders')

#Total
sim_dataX = nt + st
plot.ts(sim_dataX, ylab = 'Daily Infections count', main = 'Total - Super Spreaders Model, Daily Infections count')

#****************************************************************
# APPLY MCMC SSI MODEL
#****************************************************************
#n_mcmc = 100 #100000 #5000
mcmc_paramsI = MCMC_SSI(sim_data, n_mcmc, sigma, model_params,
                       gamma_prior, gamma_priors, DATA_AUG = FALSE)
#PLOT RESULTS
model_typeX = 'SSI'; time_elap = 0
plot_mcmc_grid(n_mcmc, sim_dataX, mcmc_paramsI, true_r0, time_elap, seed_count, model_type = model_typeX,
               flag_gam_prior_on_b = gamma_prior, gam_priors_on_b = gamma_priors, rjmcmc = RJMCMCX,
               mod_par_names = c('a', 'b', 'c'))


#****************************************************************
# APPLY MCMC SSI MODEL + DATA AUG
#***************************************************************
n_mcmc = 100000 
mcmc_params_da3 = MCMC_SSI(sim_data, n_mcmc, sigma, model_params, gamma_prior,
                       gamma_priors, DATA_AUG = TRUE)

#PLOT RESULTS
model_typeX = 'SSI'; time_elap = 0
plot_mcmc_grid(n_mcmc, sim_dataX, mcmc_params_da3, true_r0, time_elap, seed_count, model_type = model_typeX,
               flag_gam_prior_on_b = gamma_prior, gam_priors_on_b = gamma_priors, rjmcmc = RJMCMCX,
               data_aug = TRUE,
               mod_par_names = c('a', 'b', 'c'))

#DATA AUG OUPUT
data_augmented = mcmc_params_da3[[13]]
n_final = data_augmented[[1]]
s_final = data_augmented[[2]]

#DATA AUG ACCEPTED COUNTS FOR i x t
mat_da = mcmc_params_da3[[14]]
colSums(mat_da3)

#N & S FOR ALL i x t 
non_ss = mcmc_params_da[[15]]

ss = mcmc_params_da[[16]]

#**************************************
#TESTING: MCMC SIZE
n_mcmc = 10 
mcmc_params_daZ = MCMC_SSI(sim_data, n_mcmc, sigma, model_params, gamma_prior,
                           gamma_priors, DATA_AUG = TRUE)

#PLOT RESULTS
model_typeX = 'SSI'; time_elap = 0
plot_mcmc_grid(n_mcmc, sim_dataX, mcmc_params_daZ, true_r0, time_elap, seed_count, model_type = model_typeX,
               flag_gam_prior_on_b = gamma_prior, gam_priors_on_b = gamma_priors, rjmcmc = RJMCMCX,
               data_aug = TRUE,
               mod_par_names = c('a', 'b', 'c'))

#***************************************************************
#DATA AUG WITH B-C TRANSFORM
mcmc_params_da3 = MCMC_SSI(sim_data, n_mcmc, sigma, model_params, gamma_prior,
                           gamma_priors, DATA_AUG = TRUE, BC_TRANSFORM = TRUE)

#PLOT RESULTS
model_typeX = 'SSI'; time_elap = 0
plot_mcmc_grid(n_mcmc, sim_dataX, mcmc_params_da3, true_r0, time_elap, seed_count, model_type = model_typeX,
               flag_gam_prior_on_b = gamma_prior, gam_priors_on_b = gamma_priors, rjmcmc = RJMCMCX,
               data_aug = TRUE,
               mod_par_names = c('a', 'b', 'c'))

#****************************************************************#****************************************************************

#****************************************************************
#DATASET II
#****************************************************************
seed_count = 3 #seed_count = seed_count + 1 #print(paste0('i mcmc = ', i))
set.seed(seed_count)
sim_data3 = simulation_super_spreaders(num_days, shape_g, scale_g, aX, bX, cX)

#PLOTS
par(mfrow=c(1,1))
nt3 = sim_data3[[1]]
plot.ts(nt3, ylab = 'Daily Infections count', main = 'Non Super-Spreaders' )
st3 = sim_data3[[2]]
plot.ts(st3, ylab = 'Daily Infections count', main = 'Super-Spreaders')

#Total
sim_dataX3 = nt3 + st3
plot.ts(sim_dataX3, ylab = 'Daily Infections count', main = 'Total - Super Spreaders Model, Daily Infections count')

#****************************************************************
# APPLY MCMC SSI MODEL + DATA AUG
#***************************************************************
n_mcmc = 100000 
mcmc_params_da33 = MCMC_SSI(sim_data3, n_mcmc, sigma, model_params, gamma_prior,
                           gamma_priors, DATA_AUG = TRUE)

#PLOT RESULTS
model_typeX = 'SSI'; time_elap = 0
plot_mcmc_grid(n_mcmc, sim_dataX3, mcmc_params_da33, true_r0, time_elap, seed_count, model_type = model_typeX,
               flag_gam_prior_on_b = gamma_prior, gam_priors_on_b = gamma_priors, rjmcmc = RJMCMCX,
               data_aug = TRUE,
               mod_par_names = c('a', 'b', 'c'))

#DATA AUG OUPUT
data_augmented3 = mcmc_params_da33[[13]]
n_final3 = data_augmented3[[1]]
n_final3
s_final3 = data_augmented3[[2]]
s_final3

#DATA AUG ACCEPTED COUNTS FOR i x t
mat_da3 = mcmc_params_da33[[14]]
colSums(mat_da3)

#N & S FOR ALL i x t 
non_ss3 = mcmc_params_da33[[15]]
non_ss3
#colSums(non_ss)

ss3 = mcmc_params_da33[[16]]
ss3
