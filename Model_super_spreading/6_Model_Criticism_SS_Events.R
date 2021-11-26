#Model Criticism
#Super spreading Events model- MCMC Inference + Model Criticism
#Setup
library(MASS)
library(pracma)
library(tidyverse)
library(tibble)
setwd("~/GitHub/epidemic_modelling/Model_super_spreading")
source("functions.R")
#source("plotting_functions.R")

#Epidemic params
num_days = 50
#Gamma params for infectiousness curve (lambda) distribution
shape_gamma = 6
scale_gamma = 1 
#seed_count = 1
#par(mar=c(1,1,1,1))

############# --- INSERT PARAMETERS! --- ######################################
alphaX = 0.8 #0.7 #0.8 #0.7 #0.8 #0.7 #0.7 #1.1 #0.8 #1.1 #0.8 #1.1 # 0.8 #2 #0.9 #2 #2 #Without ss event, ~r0.
betaX = 0.1 #0.05 #0.025 #0.2 #0.1 #0.2 #0.05 #0.1 #0.05 #0.2 #0.05 #0.2 #0.05 #0.2 #0.2 #0.05 #0.2 #0.05 #0.05
gammaX = 10 #8
true_r0 = alphaX + betaX*gammaX
true_r0
model_params = c(alphaX, betaX, gammaX, true_r0)

#Epidemic data - Neg Bin
#sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
#plot.ts(sim_data, ylab = 'Daily Infections count', main = paste('Daily Infections count, true R0 = ', true_r0))

#MCMC - get p values 
sigma_a = 0.4*alphaX
sigma_b = 1.0*betaX 
sigma_g = 0.85*gammaX
sigma_bg = 1.5*gammaX
sigma = c(sigma_a, sigma_b, sigma_g, sigma_bg)

################################################################################
# MCMC - FOUR PARAMETER UPDATES
################################################################################
mcmc_ss_mod_crit <- function(data, n, sigma, thinning_factor, x0 = 1) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
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
  count_accept3 = 0; count_accept4 = 0; count_thin = 1
  prior = TRUE; flag_true = FALSE
  
  #MCMC chain
  for(i in 2:n) {
    
    #******************************************************
    #alpha
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma_a) 
    #cat("alpha dash: ", alpha_dash, "\n")
    
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    
    #log alpha
    logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev  #+ prior1 - prior2
    
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - alpha_dash + alpha_vec[i-1]
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
    #if(log(runif(1)) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
    }
    
    #************************************************************************
    #beta
    beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma_b) 
    #cat("Beta dash: ", beta_dash, "\n")
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }
    
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
    #gamma
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
    #Gamma-Beta
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
    
    #Get Summary Statistics
    if(mod(i, thinning_factor) == 0){
      #print('GETTING SUMMARY STATS')
      
      #If df of summary stats doesn't exist - create it
      if (!exists("df_summary_stats")) {
        #print('CREATE DF')
        flag_create = TRUE
        df_summary_stats = get_summary_stats(data, alpha_vec[i], beta_vec[i], gamma_vec[i], flag_create, flag_true)
        flag_create = FALSE
        #print('df df_summary_stats')
        #print(df_summary_stats)
        
      } else {
        df_summary_stats[nrow(df_summary_stats) + 1, ] = get_summary_stats(data, alpha_vec[i],
                                                                           beta_vec[i], gamma_vec[i], flag_create, flag_true)
      }
      
      count_thin = count_thin + 1
      
      }
  }
  
  #True summary stats - set as final row for comparison 
  flag_true = TRUE
  df_summary_stats[nrow(df_summary_stats) + 1, ] = get_summary_stats(data, alpha_vec[i], beta_vec[i], gamma_vec[i], flag_create, flag_true)
  print(df_summary_stats[nrow(df_summary_stats), ])
  
  #Final stats
  #alpha
  accept_rate1 = 100*count_accept1/n
  accept_rate2 = 100*count_accept2/n
  accept_rate3 = 100*count_accept3/n
  accept_rate4 = 100*count_accept4/n
  #cat("Acceptance rate1 = ",accept_rate1, '\n')
  
  #Get p values - comparing  summary stat columns to true value 
  list_p_vals = apply(df_summary_stats, 2, FUN = function(vec) get_p_values(vec))
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3,
              accept_rate4,
              list_p_vals))
}

#Get summary stats
get_summary_stats <- function(sim_data, alpha_vec_i, beta_vec_i, gamma_vec_i, create_df_flag, flag_true){
  
  'Get summary statisitcs of the simulated data'
  
  #Simulate data
  sim_data_params = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alpha_vec_i, beta_vec_i, gamma_vec_i)
  
  'Original data as final comparion'
  if (flag_true){
    sim_data_params = sim_data
  }
  
    if (create_df_flag){
      #Df
      summary_stats_results = data.frame(
        sumX = sum(sim_data_params),
        medianX = median(sim_data_params),
        maxX = max(sim_data_params),
        stdX = std(sim_data_params),
        val_75 = quantile(sim_data_params)[4][1][1],
        val_87_5 = mean(quantile(sim_data_params)[4][1][1], quantile(sim_data_params)[5][1][1]),
        max_dif = max(abs(diff(sim_data_params))),
        med_dif = median(abs(diff(sim_data_params))),
        mean_upper_dif = mean(c(quantile(abs(diff(sim_data_params)))[4][1][1], quantile(abs(diff(sim_data_params)))[5][1][1])),
        sum_1st_half  = sum(which(sim_data_params < quantile(sim_data_params)[3][1][1])),
        sum_2nd_half =  sum(which(sim_data_params > quantile(sim_data_params)[3][1][1]))
        
      )
      
       } else {
      #List
      summary_stats_results = list(sum(sim_data_params), median(sim_data_params), max(sim_data_params),
                               std(sim_data_params), quantile(sim_data_params)[4][1][1], 
                               mean(quantile(sim_data_params)[4][1][1], quantile(sim_data_params)[5][1][1]),
                               max(abs(diff(sim_data_params))), median(abs(diff(sim_data_params))),
                               mean(c(quantile(abs(diff(sim_data_params)))[4][1][1], quantile(abs(diff(sim_data_params)))[5][1][1])),
                               sum_1st_half  = sum(which(sim_data_params < quantile(sim_data_params)[3][1][1])),
                               sum_2nd_half =  sum(which(sim_data_params > quantile(sim_data_params)[3][1][1]))
                               )
    }
    
  summary_stats_results
  
}

#Get p values - comparing  summary stat columns to true value 
get_p_values <- function(column) {
  
  #Final val
  #cat('column length:', length(column), '\n')
  #print(column[1:10])
  last_el = column[length(column)] #True value 
  #cat('last element = ', last_el, '\n')
  #P value
  lt = length(which(column <= last_el)) #Needs to be less than or equal to 
  gt = length(which(column >= last_el)) 
  min_val = min(lt, gt)
  pvalue = min_val/length(column)
  pvalue = pvalue*2
  
  #Return p value 
  #cat('p value = ', pvalue)
  pvalue
  
}

#RUN FOR MULTIPLE REPS TO GET P VALUES
get_p_values_total <- function(n, n_reps, model_params, sigma, thinning_factor, iter, flag_ss){
  
  'Run model criticism for n_reps iterations to get a sample of p values for a number of
  different summary statistics'
  
  #Get model params
  alphaX = model_params[1]; betaX = model_params[2]
  gammaX = model_params[3]; r0 = model_params[4];
  cat('r0 = ', r0, '\n'); 
  
  #Repeat for n reps
  for(rep in 1:n_reps) {
    
    cat('\n rep =', rep, '\n')
    
    #Simulate data
    if (flag_ss){
      sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
    } else {
      sim_data = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
      #sim_data = simulation_super_spreaders(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
    }
   
    #MCMC
    mcmc_params = mcmc_ss_mod_crit(sim_data, n, sigma, thinning_factor)
    list_p_vals = mcmc_params[9]
    list_p_vals = unlist(list_p_vals)
    #cat('list_p_vals:', list_p_vals, '\n length =', length(list_p_vals))
    
    if (rep == 1) { 
      cat('p value rep = ', rep)
      #Create df; sum etc
      df_p_values = data.frame(sumX = 
        list_p_vals[1],
        medianX = list_p_vals[2],
        maxX = list_p_vals[3],
        stdX = list_p_vals[4],
        val_75 = list_p_vals[5],
        val_87_5 = list_p_vals[6],
        max_dif = list_p_vals[7],
        med_dif = list_p_vals[8],
        mean_upper_dif = list_p_vals[9],
        sum_1st_half  = list_p_vals[10],
        sum_2nd_half =  list_p_vals[11]
        
      )
      print('df_p_values')
      print(df_p_values)
      
    } else {
      
      df_p_values[nrow(df_p_values) + 1, ] = list_p_vals
    }
    
    # if(mod(rep, 100) == 0){
    #   write.csv(df_p_values, file = paste('df_p_vals_', iter, '.csv'), iter, row.names = FALSE)
    #   plot.ts(sim_data, ylab = 'Daily Infections count', main = paste('Rep {}', rep, ', Daily Infections count, true R0 = ', true_r0))
    # }
    
  }
  
  #Df
  #write.csv(df_p_values, file = paste('df_p_vals_', iter, '.csv'), iter, row.names = FALSE)
  df_p_values

}

#Plot p values
plot_p_vals <- function(df_p_vals){
  
  'Plot histograms of the p values'
  par(mfrow=c(3,4)) #c(3,4)
  
  for (i in c(1:11)){
    
    hist(df_p_vals[,i], breaks = 100, #freq = FALSE, 
         #xlim = c(xmin, xmax),
         xlab = 'p value', ylab = 'Num Samples', col = 'green',
         main = paste('', toupper(colnames(df_p_vals)[i]),', R0:', true_r0),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    #abline(v = true_sum_inf, col = 'red', lwd = 2)
  }
}

############# --- RUN P VALUES --- ######################################
iter = 4
n = 10000
n_reps = 100
thinning_factor = 50 #(1/1000)*n;
flag_ss = FALSE

#Start
start_time = Sys.time()
print('Start time:')
print(start_time)
df_p_values = get_p_values_total(n, n_reps, model_params, sigma, thinning_factor, iter, flag_ss)
cat('Time elapsed:', round(Sys.time() - start_time, 2))

#Plot
plot_p_vals(df_p_values)

#Save
write.csv(df_p_values, file = paste('df_p_vals_', iter, '.csv'), iter, row.names = FALSE)




