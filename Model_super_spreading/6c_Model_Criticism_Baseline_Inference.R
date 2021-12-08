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

mcmc_r0_model_crit_x1_rep <- function(data, n, sigma, burn_in, thinning_factor, rep, folder_results, x0 = 1) {
  
  'Returns mcmc samples of R0'
  
  #Set up
  r0_vec <- vector('numeric', n); r0_vec[1] <- x0
  U <- runif(n)
  count_accept = 0;   flag_true = FALSE
  
  #Create folder for mcmc results 
  folder_mcmc = paste0(folder_results, '/mcmc')
  ifelse(!dir.exists(file.path(folder_mcmc)), dir.create(file.path(folder_mcmc), recursive = TRUE), FALSE)
  
  #MCMC chain
  for(i in 2:n) {
    r0_dash <- r0_vec[i-1] + rnorm(1, sd = sigma) #, mean = 0, sd = sigma_opt)
    if(r0_dash < 0){
      r0_dash = abs(r0_dash)
    }
    
    #Alpha
    log_alpha = log_like(data, r0_dash) - log_like(data, r0_vec[i-1]) - r0_dash + r0_vec[i-1] #exponential prior
    #log_alpha = log_like(data, Y) - log_like(data, r0_vec[i-1]) + dgamma(Y, shape = 1, scale = 1, log = TRUE) - dgamma(r0_vec[i-1], shape = 1, scale = 1, log = TRUE) 
    
    #Metropolis Step
    if (is.na(log_alpha)){
      print('na value')
      sprintf("r0_dash: %i", r0_dash)
    }
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      r0_vec[i] <- r0_dash
      count_accept = count_accept + 1
    } else {
      r0_vec[i] <- r0_vec[i-1]
    }
    
    #Summary Stats 
    if((i > burn_in) & (mod(i, thinning_factor) == 0)){
      #print('GETTING SUMMARY STATS')
      
      #If df of summary stats doesn't exist - create it
      if (!exists("df_summary_stats")) {
        #print('CREATE DF')
        flag_create = TRUE #Create df
        df_summary_stats = get_summary_stats_baseX(data, i, r0_vec[i], flag_create, flag_true, folder_mcmc)
        flag_create = FALSE #Now set to false - so new values just added as a list 
        
        #Get indices of iterations
        list_mcmc_iters = c(i)
        
      } else {
        df_summary_stats[nrow(df_summary_stats) + 1, ] = get_summary_stats_baseX(data, i, r0_vec[i],
                                                                                     flag_create, flag_true, folder_mcmc)
        list_mcmc_iters = c(list_mcmc_iters, i)
      }
      
    }
  }
  
  #Get Summary Statistics
  #True summary stats - set as final row for comparison 
  flag_true = TRUE
  df_summary_stats[nrow(df_summary_stats) + 1, ] = get_summary_stats_baseX(data, i, r0_vec[i],
                                                                           flag_create, flag_true, folder_mcmc)
  print(df_summary_stats[nrow(df_summary_stats), ])
  
  #Save_RDS Summary stats 
  saveRDS(df_summary_stats, file = paste0(folder_results, '/df_summary_stats_', rep, ".rds"))
  saveRDS(list_mcmc_iters, file = paste0(folder_mcmc, '/i_mcmc_list_', rep, '.rds'))
  
  #SINGLE P VALUE
  #Get p values - comparing  summary stat columns to true value 
  list_p_vals = apply(df_summary_stats, 2, FUN = function(vec) get_p_values(vec))
  #saveRDS p values for rep
  name_list_p_vals  <- paste0("list_p_vals_", rep)
  saveRDS(list_p_vals, file = paste0(folder_results, '/', name_list_p_vals, ".rds"))
  
  #ALL RELATED P VALUES
  list_p_vals_list = apply(df_summary_stats, 2, FUN = function(vec) get_p_values_list(vec))
  saveRDS(list_p_vals_list, file = paste0(folder_results, '/list_all_p_vals_rep_', rep, ".rds"))
  
  #Final stats
  accept_rate = 100*(count_accept/n)
  #cat("Acceptance rate = ", accept_rate) 
  #r0_vec = r0_vec[burn_in:n]
  
  return(list(r0_vec, accept_rate, list_p_vals, df_summary_stats))
  
}

#Get summary stats
get_summary_stats_baseX <- function(sim_data, i, r0_vec_i, create_df_flag, flag_true, folder_mcmc){
  
  'Get summary statisitcs of the simulated data'
  
  #Data
  sim_data_params = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
  saveRDS(sim_data_params, file = paste0(folder_mcmc, '/sim_data_iter_', i, '.rds' ))
  
  'Original data as final comparison'
  if (flag_true){
    sim_data_params = sim_data
  }
  
  if (create_df_flag){
    #Df
    summary_stats_results = data.frame(
      sum_inf_counts = sum(sim_data_params),
      median_inf_count = median(sim_data_params),
      max_inf_count = max(sim_data_params),
      std_inf_counts = std(sim_data_params),
      val_75_infs_counts = quantile(sim_data_params)[4][1][1],
      val_87_5_infs_counts = mean(quantile(sim_data_params)[4][1][1], quantile(sim_data_params)[5][1][1]),
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
  'Get p values - comparing  summary stat columns to true value'
  
  #Final val
  last_el = column[length(column)] #True value 
  cat('last_el;', last_el)
  print("")
  num_iters = length(column) - 1
  #P value
  prop_lt = length(which(column < last_el))/num_iters + 0.5*(length(which(column == last_el)) - 1)/num_iters
  cat('prop_lt;', prop_lt)
  print("")
  prop_gt = length(which(column > last_el))/num_iters + 0.5*(length(which(column == last_el)) - 1)/num_iters
  cat('prop_gt;', prop_gt)
  print("")
  pvalue = min(prop_lt, prop_gt)
  cat('pvalue;', pvalue)
  print("")
  cat('pvalue*2;', 2*pvalue)
  print("")
  print("****************")
  
  #Return p value 
  pvalue = pvalue*2
  pvalue
  
}

#Get p values - comparing  summary stat columns to true value 
get_p_values_orig <- function(column) {
  'Get p values - comparing  summary stat columns to true value'
  
  #Final val
  last_el = column[length(column)] #True value 
  #P value
  lt = length(which(column <= last_el)) #Needs to be less than or equal to 
  gt = length(which(column >= last_el)) #Needs to be greater than or equal to 
  min_val = min(lt, gt)
  pvalue = min_val/(length(column) - 1) #Last row is the real data
  pvalue = pvalue #*2
  
  #Return p value 
  pvalue
  
}

#Save all associated p vals
get_p_values_list <- function(column){
  'Save all associated p vals' 
  
  #Final val
  last_el = column[length(column)] #True value 
  #P value
  lt = length(which(column < last_el))/(length(column) - 1) #Needs to be less than or equal to 
  gt = length(which(column > last_el))/(length(column) - 1) #Needs to be greater than or equal to
  eq = length(which(column == last_el))/(length(column) - 1)
  
  #Return p value 
  return(list(lt, gt, eq))
}

#RUN FOR MULTIPLE REPS TO GET P VALUES
get_p_values_total_base <- function(n, n_reps, model_params, sigma, thinning_factor, flag_dt, folder_results, rep, burn_in){
  
  'Run model criticism for n_reps iterations to get a sample of p values for a number of
  different summary statistics. saveRDS result in type/iter/rep subfolder'
  
  #Get model params
  alphaX = model_params[1]; betaX = model_params[2]
  gammaX = model_params[3]; r0 = model_params[4];
  #Data_type
  flag1 = flag_dt[1]; flag2 = flag_dt[2]; flag3 = flag_dt[3] 
  cat('r0 = ', r0, '\n'); 
  
  #Repeat for n reps
  for(rep in 1:n_reps) {
    
    cat('\n rep =', rep, '\n')
    #Create folder to saveRDS results of the rep
    folder_results_rep = paste0(folder_results, '/rep_', rep)
    cat('\n folder_results_rep =', folder_results_rep, '\n')
    ifelse(!dir.exists(file.path(folder_results_rep)), dir.create(file.path(folder_results_rep), recursive = TRUE), FALSE)
    
    #Simulate data
    if (flag1){
      sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
      saveRDS(sim_data, file = paste0(folder_results_rep, '/sim_data.rds'))
    } else if (flag2){
      sim_data = simulation_super_spreaders(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
      cat('simulate ss individs')
      saveRDS(sim_data, file = paste0(folder_results_rep, '/sim_data.rds'))
    } else if (flag3) {
      sim_data = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
      saveRDS(sim_data, file = paste0(folder_results_rep, '/sim_data.rds'))
      #cat('simulate_branching')
    }
    
    #MCMC  function(data, n, sigma, burn_in, thinning_factor, rep, folder_results, x0 = 1) {
    mcmc_params = mcmc_r0_model_crit_x1_rep(sim_data, n, sigma, burn_in, thinning_factor, rep, folder_results_rep)
    #Save mcmc params 
    saveRDS(mcmc_params, file = paste0(folder_results_rep, '/mcmc_params_rep_', rep, '.rds' ))
    
    list_p_vals = mcmc_params[3]
    list_p_vals = unlist(list_p_vals)
    df_summary_stats = mcmc_params[4]
    df_summary_stats = unlist(df_summary_stats)
    
    if (rep == 1) { 
      
      #Create df; sum etc
      df_p_values = data.frame(sum_inf_counts = 
                               list_p_vals[1],
                               median_inf_count = list_p_vals[2],
                               max_inf_count = list_p_vals[3],
                               std_inf_counts = list_p_vals[4],
                               val_75_infs_counts = list_p_vals[5],
                               val_87_5_infs_counts = list_p_vals[6],
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
    
  }
  
  #Ensure its a df
  df_p_values = as.data.frame(df_p_values)
  
  #SaveRDS
  saveRDS(df_p_values, file = paste0(folder_results, '/total_p_values_iter_', iter, '.rds' ))
  
  #Return
  return(list(df_p_values, df_summary_stats))
  
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
model_type = 'ss_events' #base_sim_base_inf' #'ss_ind_sse_inf'
flags_data_type = c(TRUE, FALSE, FALSE) #1) ss_events, 2) s_spreaders, 3) basline
iter = 1
folder_results = paste0('~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism/', '', model_type, '/iter_', iter)
print(folder_results)

#Repitions 
n = 10500
n_reps = 100
burn_in = 500
thinning_factor = 50 #(1/1000)*n;

#Start
start_time = Sys.time()
print('Start time:')
print(start_time)
results = get_p_values_total_base(n, n_reps, model_params, sigma, thinning_factor, flags_data_type, folder_results, rep, burn_in)
end_time = Sys.time()
time_elap = round(end_time - start_time, 2)
print('Time elapsed:')
print(time_elap)
#cat('Time elapsed:', round(Sys.time() - start_time, 2))

############ INSPECT OUTPUT #######################
#Extract
df_p_valuesC = results[[1]]
#df_p_values = unlist(df_p_values)
plot_p_vals(df_p_valuesC)

#Get p values
#df_pvals3 <- readRDS(paste0(folder_results, '/total_p_values_iter_', iter, '.rds'))
#plot_p_vals(df_pvals3)
