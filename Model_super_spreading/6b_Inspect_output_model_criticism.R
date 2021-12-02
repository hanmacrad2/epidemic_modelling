#INSPECT OUTPUT OF MODEL CRITICISM

#Set up params
source("~/GitHub/epidemic_modelling/helper_functions.R")
results_home = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism/"

#FUNCTIONS - GET DATA
get_rep_results <- function(results_home, model_type, iter, rep, true_r0,
                            list_old, list_new, upper_quant, trim_flag, list_i){
  
  #Results inspect
  results_inspect = paste0(results_home, model_type, "/iter_", iter, "/rep_", rep, '/')
  print(results_inspect)
  
  #Data
  sim_data_rep <- readRDS(paste0(results_inspect, 'ss_data.rds'))
  cat('Sum sim data = ', sum(sim_data_rep))
  df_sum_stats <- readRDS(paste0(results_inspect, 'df_summary_stats_', rep, '.rds'))
  list_p_vals <- readRDS(paste0(results_inspect, 'list_p_vals_', rep, '.rds'))
  mcmc_params <- readRDS(paste0(results_inspect, '/mcmc_params_rep_', rep, '.rds' ))
  
  #Plot p vals & summary stats
  plot_rep_results(true_r0, model_type, sim_data_rep, df_sum_stats_new, list_p_vals, upper_quant, trim_flag) 
  
  #Plot MCMC results 
  plot_mcmc_x4_priors(sim_data_rep, mcmc_params, true_r0, 'Neg Bin,', 3.0, rep, TRUE, TRUE)
  get_mcmc_runs(results_home, sim_data_rep, mcmc_params, list_i)
  
}

#Plot results 
plot_rep_results <- function(true_r0, model_type, sim_data_rep, df_sum_stats, list_p_vals, upper_quant, trim_flag){
  
 'Plot sim data, summary stats and true summary stat for a given mcmc rep' 
  #Setup
  par(mfrow = c(3,4))
  len_data = length(list_p_vals)
  colorsX <- rainbow(len_data+1)
  colors_line <- rainbow(c(15:15+len_data+1))
  
  #Sim_data
  #i.Infections
  plot.ts(sim_data_rep, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste(rep, "Day Infts SS Evnts, SS model, ", "R0 = ", true_r0), #model_type
          col = colorsX[1],
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Columns
  for (i in c(1:len_data)){
    
    X = df_sum_stats[1:nrow(df_sum_stats)-1,i]
    
    if (trim_flag){
      X = upper_quantile(X, upper_quant)
    }
   

    hist(X, breaks = 100, #freq = FALSE, 
         #xlim = c(xmin, xmax),
         xlab = paste('', toupper(colnames(df_sum_stats)[i]), '< 99th quantile'),
         ylab = 'Num Samples',
         col = colorsX[i+1],
         main = paste('', toupper(colnames(df_sum_stats)[i]),', p value:', round(list_p_vals[i],3)),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(v = df_sum_stats[nrow(df_sum_stats),i], col = colors_line[i], lwd = 2.5) #This should be true value, not p value
  }
  
}

###################################
#INSPECT MCMC
get_mcmc_runs <- function(results_home, sim_data, mcmc_params, list_idx){
  
  #Results inspect
  sim_data_path = paste0(results_home, model_type, "/iter_", iter, "/rep_", rep, '/mcmc/')
  print(sim_data_path)
  
  #Plot
  colorsX <- rainbow(length(list_idx)+1)
  par(mfrow=c(3,4))
  
  #i. Sim infections data (True for the rep)
  plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste0(rep, ", ", model_type, ", True a:", alphaX, ', b:', betaX, ', g:', gammaX, ', R0:', true_r0), #model_type
          col = colorsX[1],
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Loop - Get mcmc rep
  for (i in seq_along(list_idx)){
    get_mcmc_i(list_idx[i], mcmc_params, sim_data_path, colorsX[i])
  }
  
}

#Plot mcmc iteration data
get_mcmc_i <- function(i, mcmc_params, sim_data_path, colourX){
  
  #Sim Data
  sim_data_mcmc_rep <- readRDS(paste0(sim_data_path, 'sim_data_iter_', i, '.rds'))

  #Extract params
  alpha_mcmc = mcmc_params[1]; alpha_mcmc = unlist(alpha_mcmc)
  alpha_i = round(alpha_mcmc[i], 2)
    
  beta_mcmc = mcmc_params[2];beta_mcmc = unlist(beta_mcmc)
  beta_i = round(beta_mcmc[i],2)
  
  gamma_mcmc = mcmc_params[3]; gamma_mcmc = unlist(gamma_mcmc)
  gamma_i = round(gamma_mcmc[i], 2)
  
  r0_mcmc = mcmc_params[4]; r0_mcmc = unlist(r0_mcmc)
  r0_i = round(r0_mcmc[i], 2)
  #colourX = rainbow(1)
  
  #Plot
  plot.ts(sim_data_mcmc_rep, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste0("It_", i, ", a:", alpha_i, ', b:', beta_i, ', g:', gamma_i, ', R0:', r0_i), #model_type
          col = colourX, lwd = 2,
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
}

#Apply function to inspect specific reps
rep = 26 #14, 73
upper_quant = 1.0
trim_flag = FALSE
list_i = seq(from = 1000, to = 10000, by = 1000)
list_i
get_rep_results(results_home, model_type, iter, rep, true_r0, list_old, list_new, upper_quant, trim_flag, list_i)
