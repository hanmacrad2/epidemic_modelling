#RUN BASE MODEL MODEL CRITICISM

#Setup
setwd("~/GitHub/epidemic_modelling")
source("epidemic_functions.R")
source("helper_functions.R")
source("Model_criticism/7A_Model_Criticism_SS_Events.R")
#results_home =  "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_criticism/model_criticism_1k_I/"

#RESULT REPITIONS
n_mcmc = 5500
n_reps = 500
burn_in = 500
thinning_factor = 5 #0 #(1/1000)*n;

#########################################################
# RUN I

###############
#APPLY MCMC
model_type = 'base_inf_sse_sim' #base_inf_base_sim' #'base_inf_ssi_sim'
flags_data_type = c(TRUE, FALSE, FALSE) #1)ss_events, 2) ss_individuals, 3) basline
iter = 1
base_folder_current = paste0(results_home, model_type, '/iter_', iter)
print(base_folder_current)
# 
# #START MCMC
start_time = Sys.time()
print('Start time:')
print(start_time)
run_mcmc_base_reps(n, n_reps, true_r0, sigma_base, flags_data_type, base_folder_current, burn_in)
end_time = Sys.time()
timei = get_time(start_time, end_time)


###############
#APPLY SUMMARY STATS + p vals
model_type = 'base_inf_sse_sim'
start_time = Sys.time()
print(paste0('Time elapsed:'), start_time)
get_sum_stats_base_total(base_folder_current, n_reps)
df_p_valuesBI = get_p_values_total(base_folder_current, n_reps)
end_time = Sys.time()
#timeI = get_timeII(start_time, end_time, timei)

#PLOT
plot_p_vals(df_p_valuesBI)

############################################################
#RUN II - base_inf_base_sim

#APPLY MCMC
model_type = 'base_inf_ssi_sim' #'sse_inf_sse_sim' 'sse_inf_base_sim'
flags_data_type = c(FALSE, TRUE, FALSE) #1)ss_events, 2) ss_individuals, 3) basline
base_folder_current = paste0(results_home, model_type, '/iter_', iter)
print(base_folder_current)

#Repitions
start_time = Sys.time()
print('Start time:')
print(start_time)
run_mcmc_base_reps(n, n_reps, true_r0, sigma_base, flags_data_type, base_folder_current, burn_in)
end_time = Sys.time()
timeii = get_time(start_time, end_time)


###############
#APPLY SUMMARY STATS + p vals
model_type = 'base_inf_ssi_sim'
start_time = Sys.time()
print(paste0('Start time:', start_time))
get_sum_stats_base_total(base_folder_current, n_reps)
df_p_valuesBII = get_p_values_total(base_folder_current, n_reps)
end_time = Sys.time()
#timeII = get_timeII(start_time, end_time, timeii)

#PLOT
plot_p_vals(df_p_valuesBII)

############################################################
#RUN III - base_inf_base_sim

#APPLY MCMC
model_type = 'base_inf_base_sim' #'sse_inf_sse_sim' 'sse_inf_base_sim'
flags_data_type = c(FALSE, FALSE, TRUE) #1)ss_events, 2) ss_individuals, 3) basline
base_folder_current = paste0(results_home, model_type, '/iter_', iter)
print(base_folder_current)
# 
# #RUN MCMC
start_time = Sys.time()
print('Start time:')
print(start_time)
run_mcmc_base_reps(n_mcmc, n_reps, model_params, sigma, flags_data_type, base_folder_current, burn_in)
end_time = Sys.time()
timeiii = get_time(start_time, end_time)

###############
#APPLY SUMMARY STATS + p vals
model_type = 'base_inf_base_sim'
start_time = Sys.time()
print(paste0('Start time:', start_time))
get_sum_stats_base_total(base_folder_current, thinning_factor, n_reps, n_mcmc) 
df_p_valuesBIII = get_p_values_total(base_folder_current, n_reps) 
end_time = Sys.time()
print(end_time)
#timeIII = get_timeII(start_time, end_time, timeiii)

#PLOT
plot_p_vals(df_p_valuesBIII)
