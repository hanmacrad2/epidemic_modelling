#RUN ALL COMBINATIONS

#Setup
source("functions.R")
source("7A_Model_Criticism_SS_Events.R")
#results_home = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism_II/"
results_home = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism_II_iter_II/"

############################################################
# RUN I
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
#run_mcmc_reps(n, n_reps, model_params, sigma, flags_data_type, base_folder_current, burn_in)
end_time = Sys.time()
time_elapA = round(end_time - start_time, 2)
print(paste0('Time elapsed mcmc :', time_elapA))

###############
#APPLY SUMMARY STATS + p vals
n_reps = 69 #100
start_time = Sys.time()
print(paste0('Start time:', start_time))
get_sum_stats_total(base_folder_current, n_reps)
df_p_valuesII = get_p_values_total(base_folder_current, n_reps)
end_time = Sys.time()
#time_elapII = time_elapA + round(end_time - start_time, 2)
#print(paste0('Time elapsed total:', time_elapIII))

#PLOT
plot_p_vals(df_p_valuesII)

############################################################
#RUN III - sse_inf_base_sim

#APPLY MCMC
model_type = 'sse_inf_base_sim' #'sse_inf_sse_sim' 'sse_inf_base_sim'
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
# start_time = Sys.time()
# print('Start time:')
# print(start_time)
# #run_mcmc_reps(n, n_reps, model_params, sigma, flags_data_type, base_folder_current, burn_in)
# end_time = Sys.time()
# time_elapA = round(end_time - start_time, 2)
# print(paste0('Time elapsed mcmc :', time_elapA))

###############
#APPLY SUMMARY STATS + p vals
n_reps = 100
start_time = Sys.time()
print(paste0('Start time:', start_time))
get_sum_stats_total(base_folder_current, n_reps) 
df_p_valuesIII = get_p_values_total(base_folder_current, n_reps) 
end_time = Sys.time()
print(end_time)
#time_elapIII = time_elapA + round(end_time - start_time, 2)
#print(paste0('Time elapsed total:', time_elapIII))

#PLOT
plot_p_vals(df_p_valuesIII)
