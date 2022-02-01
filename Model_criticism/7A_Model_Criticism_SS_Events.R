#####################################################
#AIM: MODEL CRITICISM of SUPER-SPREADING EVENTS MODEL

#SETUP
library(zoo)
setwd("~/GitHub/epidemic_modelling")

#Epidemic params
num_days = 50
shape_gamma = 6 #Gamma params for infectiousness curve (lambda) distribution
scale_gamma = 1 

############# --- INSERT PARAMETERS! --- ######################################
alphaX = 0.8 #0.7 #0.8 #0.7 
betaX = 0.1 #0.05 #0.025 #0.2 #0.1 
gammaX = 10 #8
true_r0 = alphaX + betaX*gammaX
true_r0
model_params = c(alphaX, betaX, gammaX, true_r0)

#MCMC - sigma
sigma_a = 0.4*alphaX
sigma_b = 1.0*betaX 
sigma_g = 0.85*gammaX
sigma_bg = 1.5*gammaX
sigma = c(sigma_a, sigma_b, sigma_g, sigma_bg)


##############################
#1. MCMC
mcmc_ss_x4 <- function(data, n, sigma, thinning_factor, folder_results, rep, burn_in, x0 = 1, prior = TRUE) {
  
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
  flag_true = FALSE
  
  #Create folder for mcmc results 
  folder_mcmc = paste0(folder_results, '/mcmc')
  ifelse(!dir.exists(file.path(folder_mcmc)), dir.create(file.path(folder_mcmc), recursive = TRUE), FALSE)
  
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
    #End of if 
    #A final 5th move -> B = 0 if not 0, and vice versa
  }
  
  #Final stats
  accept_rate1 = 100*count_accept1/n
  accept_rate2 = 100*count_accept2/n
  accept_rate3 = 100*count_accept3/n
  accept_rate4 = 100*count_accept4/n
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3, accept_rate4))
}

##############################s
#1i. REPEAT MCMC
run_mcmc_reps_ss <- function(n, n_reps, model_params, sigma, flag_dt, base_folder, burn_in){
  
  'Run mcmc for n_reps iterations and save'
  
  #Get model params
  alphaX = model_params[1]; betaX = model_params[2]
  gammaX = model_params[3]; r0 = model_params[4];
  
  #Data_type
  flag1 = flag_dt[1]; flag2 = flag_dt[2]; flag3 = flag_dt[3] 
  cat('r0 = ', r0, '\n');
  
  if (flag1){
    start_rep = 462
  } else {
    start_rep = 1
  }
  
  #Repeat for n reps
  
  for(rep in start_rep:n_reps) {
    
    cat('\n rep =', rep, '\n')
    folder_rep = paste0(base_folder, '/rep_', rep)
    ifelse(!dir.exists(file.path(folder_rep)), dir.create(file.path(folder_rep), recursive = TRUE), FALSE)
    
    #Simulate data
    if (flag1){
      sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
      print('simulate ss events model')
      saveRDS(sim_data, file = paste0(folder_rep, '/sim_data.rds'))
    } else if (flag2){
      sim_data = simulation_super_spreaders(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
      print('simulate ss individuals model')
      saveRDS(sim_data, file = paste0(folder_rep, '/sim_data.rds'))
    } else if (flag3) {
      sim_data = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
      saveRDS(sim_data, file = paste0(folder_rep, '/sim_data.rds'))
      print('simulate base model')
    }
    
    #MCMC
    mcmc_params = mcmc_ss_x4(sim_data, n, sigma, thinning_factor, folder_rep, rep, burn_in)
      
    #SAVE MCMC PARAMS 
    saveRDS(mcmc_params, file = paste0(folder_rep, '/mcmc_params_rep_', rep, '.rds' ))
    
  }
  
}

#*******************************************************#******************************************
#*
########################################
#2.  MODEL CRITICISM - GET SUMMARY STATS

##############################
#1. P VALUES FUCNTION
get_p_values <- function(col_sum_stat, col_true_val) {
  'Get p values - comparing  summary stat columns to true value'
  
  #Final val
  num_iters = length(col_sum_stat)# - 1
  #P value
  prop_lt = length(which(col_sum_stat < col_true_val))/num_iters + 0.5*(length(which(col_sum_stat == col_true_val)))/num_iters
  prop_gt = length(which(col_sum_stat > col_true_val))/num_iters + 0.5*(length(which(col_sum_stat == col_true_val)))/num_iters
  pvalue = min(prop_lt, prop_gt)
  pvalue = pvalue*2
  
  #Return p value 
  pvalue
  
}

#Save all p vals
get_p_values_list <- function(col_sum_stat, col_true_val){
  'Save all associated p vals' 
  
  #Final val
  num_iters = length(col_sum_stat)# - 1
  #P value
  prop_lt = length(which(col_sum_stat < col_true_val))/num_iters + 0.5*(length(which(col_sum_stat == col_true_val)))/num_iters
  prop_gt = length(which(col_sum_stat > col_true_val))/num_iters + 0.5*(length(which(col_sum_stat == col_true_val)))/num_iters
  pvalue = min(prop_lt, prop_gt)
  eq_half = 0.5*(length(which(col_sum_stat == col_true_val)))
  
  #Return p value 
  return(list(prop_lt, prop_gt, pvalue, eq_half))
}

#################
#2. SUMMARY STATS
get_summary_stats <- function(data, flag_create){
  
  'Calculate summary statisitcs of the simulated data'
  #Summary stats params
  start_half2 = (length(data)/2)+1
  stop_half2 = length(data)
  
  if (flag_create){
    
    #Df of summary stats (20)
    summary_stats_results = data.frame( 
      
      sum_infects = sum(data),
      sum_1st_half  = sum(data[1:(length(data)/2)]),
      sum_2nd_half =  sum(data[start_half2:stop_half2]),
      
      median_infect = median(data),
      max_infect = max(data),
      sd_infect = sd(data),
      
      infect_q_75 = quantile(data)[4][1][1],
      infect_q_87_5 = quantile(data, probs = seq(0, 1, 0.125))[8][1][1],
      
      #Differences
      max_dif = max((diff(data))), #Change from absolute difference
      med_dif = median(diff(data)),
      max_dif2nd = max(diff(diff(data))),
      med_dif2nd = median(diff(diff(data))),
      
      #Norm Differences
      max_dif_normI = max(diff(data)/(data[1:length(data)-1]+1)),
      max_dif_normII = max(diff(data)/(rollapply(data, 2, mean, by = 1)+1)),
      max_dif_normIII = max(diff(data)/(data[1:length(data)-1]/(data[2:length(data)]+1)+1)),
      
      max_dif2nd_I = max(diff(diff(data))/(data[1:length(data)-1]+1)),
      max_dif_2ndII = max(diff(diff(data))/(rollapply(data, 2, mean, by = 1)+1)),
      
      med_dif_normI = median(diff(data)/(data[1:length(data)-1]+1)),
      med_dif_normII = median(diff(data)/(rollapply(data, 2, mean, by = 1) +1)),
      med_dif_normIII = median(diff(data)/(data[1:length(data)-1]/(data[2:length(data)]+1)+1))
      
      #med_dif2nd_I = median(diff(diff(data))/(data[1:length(data)-1]+1)),
      #med_dif_2ndII = median(diff(diff(data))/(rollapply(data, 2, mean, by = 1)+1))
      
    )
    
    
    
  } else {
    #List
    summary_stats_results = list(sum(data), sum(data[1:(length(data)/2)]), sum(data[start_half2:stop_half2]),
                                 median(data), max(data), sd(data),
                                 quantile(data)[4][1][1], quantile(data, probs = seq(0, 1, 0.125))[8][1][1],
                                 max(diff(data)), median(diff(data)), max(diff(diff(data))), median(diff(diff(data))),
                                 max(diff(data)/(data[1:length(data)-1]+1)), max(diff(data)/(rollapply(data, 2, mean, by = 1)+1)),
                                 max(diff(data)/(data[1:length(data)-1]/(data[2:length(data)]+1)+1)),
                                 max(diff(diff(data))/(data[1:length(data)-1]+1)),  max(diff(diff(data))/(rollapply(data, 2, mean, by = 1)+1)),
                                 median(diff(data)/(data[1:length(data)-1]+1)), median(diff(data)/(rollapply(data, 2, mean, by = 1)+1)),
                                 median(diff(data)/(data[1:length(data)-1]/(data[2:length(data)]+1)+1))
                                 
                                 
    )
  }
  
  summary_stats_results
  
}

get_summary_stats_I <- function(data, flag_create){
  
  'Calculate summary statisitcs of the simulated data'
  #Summary stats params
  start_half2 = (length(data)/2)+1
  stop_half2 = length(data)
  
  #Sum stats divide by point in time
  max_dif = max((diff(data))   ) #Change from absolute difference
  med_dif = median(diff(data))
  dif_dif = diff(diff(data))
  
  if (flag_create){
    
    #Df
    summary_stats_results = data.frame(
      sum_infects = sum(data),
      median_infect = median(data),
      max_infect = max(data),
      sd_infect = sd(data),
      
      infect_q_75 = quantile(data)[4][1][1],
      infect_q_87_5 = quantile(data, probs = seq(0, 1, 0.125))[8][1][1], #mean(quantile(data)[4][1][1], quantile(data)[5][1][1]),
      max_dif = max((diff(data))), #Change from absolute difference
      med_dif = median(diff(data)),
      dif_dif = diff(diff(data)),
      dif_q_87_5 = quantile(diff(data), probs = seq(0, 1, 0.125))[8][1][1],
      #mean_upper_dif = mean(c(quantile(abs(diff(data)))[4][1][1], quantile(abs(diff(data)))[5][1][1])),
      sum_1st_half  = sum(data[1:(length(data)/2)]), #sum(which(data < quantile(data)[3][1][1]))
      sum_2nd_half =  sum(data[start_half2:stop_half2]) #sum(which(data > quantile(data)[3][1][1]))
    )
    
  } else {
    #List
    summary_stats_results = list(sum(data), median(data), max(data),
                                 sd(data), quantile(data)[4][1][1], 
                                 mean(quantile(data)[4][1][1], quantile(data)[5][1][1]),
                                 max(abs(diff(data))), median(abs(diff(data))),
                                 mean(c(quantile(abs(diff(data)))[4][1][1], quantile(abs(diff(data)))[5][1][1])),
                                 sum(data[1:(length(data)/2)]),
                                 sum(data[start_half2:stop_half2]) #sum(which(data > quantile(data)[3][1][1]))
    )
  }
  
  summary_stats_results
 
}

##############################
#2B. TOTAL SUMMARY STATS 
get_sum_stats_total <- function(base_folder_current, thinning_factor, n_reps, n_mcmc){
  
  'Get summary stats and p vals for all mcmc reps'
  
  for(rep in 1:n_reps) {
    
    #Get results
    folder_rep = paste0(base_folder_current, "/rep_", rep, '/')
    cat('rep = ', rep)
    true_rep_sim = readRDS(paste0(folder_rep, '/sim_data.rds'))
    mcmc_params <- readRDS(paste0(folder_rep, '/mcmc_params_rep_', rep, '.rds' ))
    #Get true summary statistics 
    df_true_ss = get_summary_stats(true_rep_sim, TRUE)
    #Save 
    saveRDS(df_true_ss, file = paste0(folder_rep, 'df_true_sum_stats_rep_', rep, '.rds' ))
    
    #Get parameters
    alpha_mcmc = mcmc_params[1]; alpha_mcmc = unlist(alpha_mcmc)
    beta_mcmc = mcmc_params[2]; beta_mcmc = unlist(beta_mcmc)
    gamma_mcmc = mcmc_params[3]; gamma_mcmc = unlist(gamma_mcmc)
    r0_mcmc = mcmc_params[4]; r0_mcmc = unlist(r0_mcmc)
    
    #Simulate data using thinned params
    for(i in seq(burn_in, n_mcmc, by = thinning_factor)){
      #print(paste0("mcmc summary stat rep ", i))
      #Simulate data
      sim_data_model_crit = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alpha_mcmc[i], beta_mcmc[i], gamma_mcmc[i])
      #Save data
      saveRDS(sim_data_model_crit, file = paste0(folder_rep, 'mcmc/sim_data_iter_', i, '.rds' ))
      
      #Get summary stats. 
      if (i == burn_in) { #first rep
        #cat('CREATE  df_summary_stats')
        flag_create = TRUE
        df_summary_stats = get_summary_stats(sim_data_model_crit, flag_create)
        flag_create = FALSE 
        #Get indices of iterations
        list_ss_iters = c(i)
      } else {
        df_summary_stats[nrow(df_summary_stats) + 1, ] = get_summary_stats(sim_data_model_crit, flag_create)
        list_ss_iters = c(list_ss_iters, i)
      }
      
    }
    
    #Save summary stats
    saveRDS(df_summary_stats, file = paste0(folder_rep, '/df_summary_stats_', rep, ".rds"))
    #print(paste0('df_summary_stats', df_summary_stats))
    #Save ss iterations
    saveRDS(list_ss_iters, file = paste0(folder_rep, '/list_ss_iters_i', rep, '.rds'))  

    }
}

##############################
#3. GET P VALUES FOR ALL  REPS
get_p_values_total <- function(base_folder_current, n_reps){
  
  for(rep in 1:n_reps) {
    cat('rep = ', rep)
    #Get results
    folder_rep = paste0(base_folder_current, "/rep_", rep, '/')
    cat('folder_rep', folder_rep)
    true_rep_sim = readRDS(paste0(folder_rep, '/sim_data.rds'))
    #Get true summary statistics 
    df_true_ss = get_summary_stats(true_rep_sim, TRUE)
    
    #Data
    df_summary_stats_rep <- readRDS(paste0(folder_rep, '/df_summary_stats_', rep, '.rds' ))
    
    #Get p values
    list_p_vals = sapply(1:ncol(df_summary_stats_rep), function(x) get_p_values(df_summary_stats_rep[,x], df_true_ss[,x]))
    saveRDS(list_p_vals, file = paste0(folder_rep, '/list_p_vals_rep', rep, ".rds"))
    
    list_all_p_vals = sapply(1:ncol(df_summary_stats_rep), function(x) get_p_values_list(df_summary_stats_rep[,x], df_true_ss[,x]))
    saveRDS(list_all_p_vals, file = paste0(folder_rep, '/list_all_p_vals_rep_', rep, ".rds"))
    
    #Save all 
    if (!exists("df_p_values")) {
      df_p_values = data.frame(sum_infects = list_p_vals[1],
                               sum_1st_half = list_p_vals[2],
                               sum_2nd_half = list_p_vals[3],
                               median_infect = list_p_vals[4],
                               max_infect = list_p_vals[5],
                               sd_infect = list_p_vals[6],
                               infect_q_75 = list_p_vals[7],
                               infect_q_87_5 = list_p_vals[8],
                               max_dif = list_p_vals[9],
                               med_dif  = list_p_vals[10],
                               max_dif2nd =  list_p_vals[11],
                               med_dif2nd =  list_p_vals[12],
                               max_dif_normI =  list_p_vals[13],
                               max_dif_normII =  list_p_vals[14],
                               max_dif_normIII = list_p_vals[15],
                               max_dif2nd_I =  list_p_vals[16],
                               max_dif_2ndII =  list_p_vals[17],
                               med_dif_normI =  list_p_vals[18],
                               med_dif_normII =  list_p_vals[19],
                               med_dif_normIII =  list_p_vals[20]
                                
      )
      print(paste0('df_p_values', df_p_values))
      
    } else {
      df_p_values[nrow(df_p_values) + 1, ] = list_p_vals
    }
    
  }
  
  #Ensure its a df
  df_p_values = as.data.frame(df_p_values)
  
  #SaveRDS
  saveRDS(df_p_values, file = paste0(base_folder_current, '/total_p_values_iter_', iter, '.rds' ))
  
  #Return p values
  df_p_values
  
}

get_p_values_total_I <- function(base_folder_current, n_reps){

    for(rep in 1:n_reps) {
      cat('rep = ', rep)
      #Get results
      folder_rep = paste0(base_folder_current, "/rep_", rep, '/')
      cat('folder_rep', folder_rep)
      true_rep_sim = readRDS(paste0(folder_rep, '/sim_data.rds'))
      #Get true summary statistics 
      df_true_ss = get_summary_stats(true_rep_sim, TRUE)
      
      #Data
      df_summary_stats_rep <- readRDS(paste0(folder_rep, '/df_summary_stats_', rep, '.rds' ))
      
      #Get p values
      list_p_vals = sapply(1:ncol(df_summary_stats_rep), function(x) get_p_values(df_summary_stats_rep[,x], df_true_ss[,x]))
      saveRDS(list_p_vals, file = paste0(folder_rep, '/list_p_vals_rep', rep, ".rds"))
      
      list_all_p_vals = sapply(1:ncol(df_summary_stats_rep), function(x) get_p_values_list(df_summary_stats_rep[,x], df_true_ss[,x]))
      saveRDS(list_all_p_vals, file = paste0(folder_rep, '/list_all_p_vals_rep_', rep, ".rds"))
      
      #Save all 
      if (!exists("df_p_values")) {
        df_p_values = data.frame(sum_inf_counts = list_p_vals[1],
                                 median_inf_count = list_p_vals[2],
                                 max_inf_count = list_p_vals[3],
                                 sd_inf_counts = list_p_vals[4],
                                 val_75_infs_counts = list_p_vals[5],
                                 val_87_5_infs_counts = list_p_vals[6],
                                 max_dif = list_p_vals[7],
                                 med_dif = list_p_vals[8],
                                 mean_upper_dif = list_p_vals[9],
                                 sum_1st_half  = list_p_vals[10],
                                 sum_2nd_half =  list_p_vals[11]
                                 
        )
        print(paste0('df_p_values', df_p_values))
        
      } else {
        df_p_values[nrow(df_p_values) + 1, ] = list_p_vals
      }
        
    }
  
  #Ensure its a df
  df_p_values = as.data.frame(df_p_values)
  
  #SaveRDS
  saveRDS(df_p_values, file = paste0(base_folder_current, '/total_p_values_iter_', iter, '.rds' ))
  
  #Return p values
  df_p_values
  
}

#############################
#4.PLOT P VALUES
plot_p_vals0 <- function(df_p_vals){
  
  'Plot histograms of the p values'
  par(mfrow=c(4,5)) #c(3,4)
  
  #Prop lt than 0.05
  num_iters = length(df_p_vals[,1])
  
  for (i in c(1:20)){ #1:11
    
    #Prop lt 0.05
    val_05 = 0.05
    percent_lt_05 = (length(which(df_p_vals[,i] < val_05))/num_iters)*100
    
    hist(df_p_vals[,i], breaks = 100, #freq = FALSE, 
         #xlim = c(xmin, xmax),
         xlab = paste0('p value, < 0.05: ', percent_lt_05, '%'),
         ylab = 'Num Samples', col = 'green',
         main = paste('', toupper(colnames(df_p_vals)[i]),', R0:', true_r0, ', n reps:', num_iters),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(v = 0.05, col = 'red', lwd = 2)
  }
}
