#Setup
source("functions.R")
source("7A_Model_Criticism_SS_Events.R")
#results_home = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism_II/"
#results_home = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism_II_iter_II/"
results_home = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism_II_iter_II/base/"
source("7A_Model_Criticism_SS_Events.R")
source("7B_Run_model_criticism_all.R")
source("7C_Inspect_output_model_criticism.R")

#params
sigma_base = 0.5

############################################################
# 1. MCMC
run_mcmc_base_reps <- function(n, n_reps, r0, sigma, flag_dt, base_folder, burn_in){
  
  'Run mcmc for n_reps iterations and save'
  
  #Data_type
  flag1 = flag_dt[1]; flag2 = flag_dt[2]; flag3 = flag_dt[3] 
  cat('r0 = ', r0, '\n'); 
  
  for(rep in 1:n_reps) {
    
    #Rep folder
    cat('\n rep =', rep, '\n')
    folder_rep = paste0(base_folder, '/rep_', rep)
    ifelse(!dir.exists(file.path(folder_rep)), dir.create(file.path(folder_rep), recursive = TRUE), FALSE)
    #MCMC folder - mcmc results 
    folder_mcmc = paste0(folder_rep, '/mcmc')
    ifelse(!dir.exists(file.path(folder_mcmc)), dir.create(file.path(folder_mcmc), recursive = TRUE), FALSE)
    
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
    mcmc_params = mcmc_r0(sim_data, n, sigma, burn_in)  ###!!
    #mcmc_ss_x4(sim_data, n, sigma, thinning_factor, folder_rep, rep, burn_in)
    
    #SAVE MCMC PARAMS 
    saveRDS(mcmc_params, file = paste0(folder_rep, '/mcmc_params_rep_', rep, '.rds' ))
    
  }
  
}

############################################################
# 2. SUMMARY STATS
get_sum_stats_base_total <- function(base_folder_current, n_reps){
  
  'Get summary stats and p vals for all mcmc reps'
  
  for(rep in 1:n_reps) {
    
    #Get results
    folder_rep = paste0(base_folder_current, "/rep_", rep, '/')
    cat('rep = ', rep)
    true_rep_sim = readRDS(paste0(folder_rep, 'sim_data.rds'))
    mcmc_params <- readRDS(paste0(folder_rep, 'mcmc_params_rep_', rep, '.rds' ))
    #Get true summary statistics 
    df_true_ss = get_summary_stats(true_rep_sim, TRUE)
    #Save 
    saveRDS(df_true_ss, file = paste0(folder_rep, 'df_true_sum_stats_rep_', rep, '.rds' ))
    
    #Get parameters
    r0_mcmc = mcmc_params[1]; r0_mcmc = unlist(r0_mcmc)
    r0_na_count = 0
    
    #MCMC folder - mcmc results #REMOVE FROM FUNCTION
    folder_mcmc = paste0(folder_rep, '/mcmc')
    ifelse(!dir.exists(file.path(folder_mcmc)), dir.create(file.path(folder_mcmc), recursive = TRUE), FALSE)
    
    #Simulate data using thinned params
    for(i in seq(burn_in, n, by = thinning_factor)){
      #cat('r0 i = ', r0_mcmc[i])
      #Check
      if (is.na(r0_mcmc[i])){
        print('r0 na value')
        r0_mcmc[i] = 0
        r0_na_count = r0_na_count + 1
      }
      
      #Simulate data
      sim_data_model_crit = simulate_branching(num_days, r0_mcmc[i], shape_gamma, scale_gamma)
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
    print(paste0('r0_na_count = ', r0_na_count))
    saveRDS(df_summary_stats, file = paste0(folder_rep, '/df_summary_stats_', rep, ".rds"))
    #print(paste0('df_summary_stats', df_summary_stats))
    #Save ss iterations
    saveRDS(list_ss_iters, file = paste0(folder_rep, '/list_ss_iters_i', rep, '.rds'))  
    
  }
}


############################################################
# RUN I

###############
#APPLY MCMC
model_type = 'base_inf_sse_sim' #base_inf_base_sim' #'base_inf_ssi_sim'
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
# run_mcmc_base_reps(n, n_reps, true_r0, sigma_base, flags_data_type, base_folder_current, burn_in)
# end_time = Sys.time()
# total_time_elap = round(end_time - start_time, 2)
# print('total_time_elap:')
# print(total_time_elap)

###############
#APPLY SUMMARY STATS + p vals
n_reps = 100
start_time = Sys.time()
print('Start time:')
print(start_time)
get_sum_stats_base_total(base_folder_current, n_reps)
df_p_valuesBI = get_p_values_total(base_folder_current, n_reps)
end_time = Sys.time()
time_elapI = round(difftime(end_time, start_time, units='hours'), 2) #round(end_time - start_time, 2)
print(paste0('Time elapsed:', time_elapI))

#PLOT
plot_p_vals(df_p_valuesBI)

############################################################
#RUN II - sse_inf_ssi_sim

#APPLY MCMC
model_type = 'base_inf_ssi_sim' #'base_inf_sse_sim' #base_inf_base_sim' #''
flags_data_type = c(FALSE, TRUE, FALSE) #1)ss_events, 2) ss_individuals, 3) basline
iter = 1
base_folder_current = paste0(results_home, model_type, '/iter_', iter)
print(base_folder_current)

#Repitions
n = 5500
n_reps = 100
burn_in = 500
thinning_factor = 50 #(1/1000)*n;
start_time = Sys.time()
print('Start time:')
print(start_time)
run_mcmc_base_reps(n, n_reps, true_r0, sigma_base, flags_data_type, base_folder_current, burn_in)
end_time = Sys.time()
time_elapA = round(end_time - start_time, 2)
print(paste0('Time elapsed mcmc :', time_elapA))

###############
#APPLY SUMMARY STATS + p vals
n_reps = 100
start_time = Sys.time()
print(paste0('Start time:', start_time))
get_sum_stats_base_total(base_folder_current, n_reps)
df_p_valuesBII = get_p_values_total(base_folder_current, n_reps)
end_time = Sys.time()
time_elapII = time_elapA + round(end_time - start_time, 2)
print(paste0('Time elapsed total:', time_elapII))

#PLOT
plot_p_vals(df_p_valuesBII)

############################################################
#RUN III - sse_inf_base_sim

#APPLY MCMC
model_type = 'base_inf_base_sim' 
flags_data_type = c(FALSE, FALSE, TRUE) #1)ss_events, 2) ss_individuals, 3) basline
iter = 1
base_folder_current = paste0(results_home, model_type, '/iter_', iter) 
print(base_folder_current)

#Repitions 
n = 5500
n_reps = 100
burn_in = 500
thinning_factor = 50 #(1/1000)*n;

#RUN MCMC
start_time = Sys.time()
print('Start time:')
print(start_time)
#run_mcmc_base_reps(n, n_reps, true_r0, sigma_base, flags_data_type, base_folder_current, burn_in)
end_time = Sys.time()
time_elapA = round(end_time - start_time, 2)
print(paste0('Time elapsed mcmc :', time_elapA))

###############
#APPLY SUMMARY STATS + p vals
n_reps = 100
start_time = Sys.time()
print(paste0('Start time:', start_time))
get_sum_stats_base_total(base_folder_current, n_reps)
df_p_valuesBIII = get_p_values_total(base_folder_current, n_reps)
end_time = Sys.time()
print(end_time)
time_elapIII = time_elapA + round(end_time - start_time, 2)
print(paste0('Time elapsed total:', time_elapIII))

#PLOT
plot_p_vals(df_p_valuesBIII)


###########################################################
#******************
#*RUN SSE-SSI
results_home = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism_II_iter_II/"

############################################################
#RUN II - sse_inf_ssi_sim

#APPLY MCMC
model_type = 'sse_inf_ssi_sim' #'sse_inf_sse_sim' 'sse_inf_base_sim'
flags_data_type = c(FALSE, TRUE, FALSE) #1)ss_events, 2) ss_individuals, 3) basline
iter = 1
base_folder_current = paste0(results_home, model_type, '/iter_', iter)
print(base_folder_current)

#Repitions
n = 5500
n_reps = 100
burn_in = 500
thinning_factor = 50 #(1/1000)*n;
start_time = Sys.time()
print('Start time:')
print(start_time)
run_mcmc_reps(n, n_reps, model_params, sigma, flags_data_type, base_folder_current, burn_in)
end_time = Sys.time()
time_elapA = round(end_time - start_time, 2)
print(paste0('Time elapsed mcmc :', time_elapA))

###############
#APPLY SUMMARY STATS + p vals
n_reps = 100
start_time = Sys.time()
print(paste0('Start time:', start_time))
get_sum_stats_total(base_folder_current, n_reps)
df_p_valuesII = get_p_values_total(base_folder_current, n_reps)
end_time = Sys.time()
#time_elapII = time_elapA + round(end_time - start_time, 2)
print(paste0('Time elapsed total:', time_elapIII))

#PLOT
plot_p_vals(df_p_valuesII)

#CHANGE run_mcmc_reps BACK TO 1:n_reps