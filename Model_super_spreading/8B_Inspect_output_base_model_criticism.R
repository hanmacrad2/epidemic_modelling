#INSPECT OUTPUT OF MODEL CRITICISM

#Set up params
setwd("~/GitHub/epidemic_modelling/Model_super_spreading")
source("functions.R")
source("~/GitHub/epidemic_modelling/helper_functions.R") 
source("7A_Model_Criticism_SS_Events.R")
source("7B_Run_model_criticism_all.R")
source("7C_Inspect_output_model_criticism.R")
#results_home = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism_II/"
results_home = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism_II_iter_II/"

#*##########################################################
#1. GET & DISPLAY TOTAL REP RESULTs - BASE MODEL
display_rep_results_base <- function(results_home, model_type, iter, rep, n_mcmc, true_r0,
                                     upper_quant, trim_flag, list_i, time_elap){
  
  #Results inspect
  results_inspect = paste0(results_home, model_type, "/iter_", iter, "/rep_", rep, '/')
  print(results_inspect)
  
  #Data
  sim_data_rep <- readRDS(paste0(results_inspect, 'sim_data.rds')) 
  mcmc_params <- readRDS(paste0(results_inspect, '/mcmc_params_rep_', rep, '.rds' ))
  df_sum_stats <- readRDS(paste0(results_inspect, 'df_summary_stats_', rep, '.rds'))
  df_true_sum_stats <- readRDS(paste0(results_inspect, 'df_true_sum_stats_rep_', rep, '.rds' ))
  list_p_vals <- readRDS(paste0(results_inspect, 'list_p_vals_rep', rep, '.rds'))
  print(df_true_sum_stats)
  
  #Plot p vals & summary stats
  plot_rep_sum_stats(true_r0, model_type, sim_data_rep, df_sum_stats, df_true_sum_stats,
                     list_p_vals, upper_quant, trim_flag) 
  
  #Plot MCMC results 
  plot_mcmc_results_r0(n_mcmc, sim_data_rep, mcmc_params, true_r0, time_elap, rep, model_type)
  get_sim_data_mcmc_runs_base(results_home, sim_data_rep, mcmc_params, list_i,model_type, rep)
  
} 

#BASE MODEL - Simulations from MCMC 
get_sim_data_mcmc_runs_base <- function(results_home, sim_data, mcmc_params, list_idx, model_type, rep){
  
  #Plot
  coloursX <- rainbow(length(list_idx)+1)
  par(mfrow=c(3,4))
  
  #i. Sim infections data (True for the rep)
  plot.ts(sim_data, xlab = 'Time', ylab = 'True Daily Infections',
          main = paste0(rep, ", True Infts, ", model_type, ", R0 = ", true_r0), #model_type
          col = coloursX[1],
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Results inspect
  sim_data_path = paste0(results_home, model_type, "/iter_", iter, "/rep_", rep, '/mcmc/')
  print(sim_data_path)
  

  #Loop - Get mcmc rep
  for (i in (seq_along(list_idx))){ #seq_along
    #Sim Data
    iter_mcmc = list_idx[i]
    print(paste0('i =', i))
    
    # if(i == length(list_idx)){
    #   iter_mcmc = list_idx[i] - 50
    # }
  
    sim_data_mcmc_rep <- readRDS(paste0(sim_data_path, 'sim_data_iter_', iter_mcmc, '.rds'))
    print(paste0('sum of sim_data_mcmc_rep = ', sum(sim_data_mcmc_rep)))
    r0_mcmc = mcmc_params[1]; r0_mcmc = unlist(r0_mcmc)
    r0_i = round(r0_mcmc[iter_mcmc], 2)
    
    #Plot
    #plot.ts(sim_data_mcmc_rep)
    plot.ts(sim_data_mcmc_rep, xlab = 'Time', ylab = 'Daily Infections count',
            main = paste0("It_", iter_mcmc, ', R0 inferred:', r0_i), #model_type
            col = coloursX[i], lwd = 2,
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  }
  
}


##############################################
#SSE MODEL - INSPECT SPECIFIC REPS
#MODELs x3 APPLY - INSPECT SPECIFIC REPS
model_type = 'base_inf_sse_sim' #'base_inf_ssi_sim'  #base_inf_base_sim' #''
iter = 1
n_mcmc = 5500
n_reps = 100
#Df p vals
df_p_vals_baseI = get_df_p_vals(results_home, model_type, iter)
plot_p_vals(df_p_vals_baseI)
upper_quant = 0.99 #1.0
trim_flag = FALSE #TRUE #
list_i = seq(from = 500, to = 5500, by = 500)

#Rep
list_i = seq(from = 500, to = 5000, by = 500)
rep = 16 
time_elap = 0.5
display_rep_results_base(results_home, model_type, iter, rep, n_mcmc, true_r0,
                    upper_quant, trim_flag, list_i, time_elap)


#
