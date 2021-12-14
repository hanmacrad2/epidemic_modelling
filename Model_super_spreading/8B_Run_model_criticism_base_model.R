#Setup
source("functions.R")
source("7A_Model_Criticism_SS_Events.R")
#results_home = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism_II/"
results_home = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism_II_iter_II/"

############################################################
# RUN

run_mcmc_base_reps <- function(n, n_reps, model_params, sigma, flag_dt, base_folder, burn_in){
  
  'Run mcmc for n_reps iterations and save'
  
  #Get model params
  alphaX = model_params[1]; betaX = model_params[2]
  gammaX = model_params[3]; r0 = model_params[4];
  
  #Data_type
  flag1 = flag_dt[1]; flag2 = flag_dt[2]; flag3 = flag_dt[3] 
  cat('r0 = ', r0, '\n'); 
  
  #Repeat for n repswhich(col_sum_stat < col_true_val) 
  
  for(rep in 1:n_reps) {
    
    cat('\n rep =', rep, '\n')
    folder_rep = paste0(base_folder, '/rep_', rep)
    ifelse(!dir.exists(file.path(folder_rep)), dir.create(file.path(folder_rep), recursive = TRUE), FALSE)
    
    #Simulate data
    if (flag1){
      sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
      saveRDS(sim_data, file = paste0(folder_rep, '/sim_data.rds'))
    } else if (flag2){
      sim_data = simulation_super_spreaders(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
      cat('simulate ss individs')
      saveRDS(sim_data, file = paste0(folder_rep, '/sim_data.rds'))
    } else if (flag3) {
      sim_data = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
      saveRDS(sim_data, file = paste0(folder_rep, '/sim_data.rds'))
      cat('simulate_branching')
    }
    
    #MCMC
    mcmc_params = mcmc_r0 ###!!
    #mcmc_ss_x4(sim_data, n, sigma, thinning_factor, folder_rep, rep, burn_in)
    
    #SAVE MCMC PARAMS 
    saveRDS(mcmc_params, file = paste0(folder_rep, '/mcmc_params_rep_', rep, '.rds' ))
    
  }
  
}


###############
#APPLY MCMC
model_type = 'sse_inf_sse_sim' #base_sim_sse_inf' #'ssi_sim_sse_inf'
flags_data_type = c(TRUE, FALSE, FALSE) #1)ss_events, 2) ss_individuals, 3) basline  
iter = 1
base_folder_current = paste0(results_home, model_type, '/iter_', iter) 
print(base_folder_current)

#Repitions 
n_mcmc = 5500
n_reps = 100
burn_in = 500
thinning_factor = 50 #(1/1000)*n;

#START MCMC
# start_time = Sys.time()
# print('Start time:')
# print(start_time)
# #run_mcmc_reps(n, n_reps, model_params, sigma, flags_data_type, base_folder_current, burn_in)
# end_time = Sys.time()
# total_time_elap = round(end_time - start_time, 2)
# print('Time elapsed:')
# print(time_elap)

###############
#APPLY SUMMARY STATS + p vals
n_reps = 100
start_time = Sys.time()
print('Start time:')
print(start_time)
get_sum_stats_total(base_folder_current, n_reps) 
df_p_valuesI = get_p_values_total(base_folder_current, n_reps) 
end_time = Sys.time()
time_elapI = round(difftime(end_time, start_time, units='hours'), 2) #round(end_time - start_time, 2)
print(paste0('Time elapsed:', time_elapI))

#PLOT
plot_p_vals(df_p_valuesI)
#plot_p_vals(df_p_vals_si)
