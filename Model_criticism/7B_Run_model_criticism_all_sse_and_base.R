#RUN ALL COMBINATIONS (SS Events & BASE)

#Setup
setwd("~/GitHub/epidemic_modelling")
source("epidemic_functions.R")
source("helper_functions.R")
source("Model_criticism/7A_Model_Criticism_SS_Events.R")
source("Model_criticism/8A_Model_criticism_base_model.R")
#results_folder =  "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_criticism/model_criticism_1k_I/"
results_folder =  "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_criticism/model_criticism_1k_II/"

#################
#RESULT REPITIONS
iter = 1
n_mcmc = 5500
n_reps = 1000
burn_in = 500
thinning_factor = 5 #0 #(1/1000)*n;

#### - MCMC params - ######
alphaX = 1.1 #0.8 #0.7 #0.8 #0.7 
betaX = 0.2 #0.05 #0.025 #0.2 #0.1 
gammaX = 2 #10 #8 #TRY WITH SMALLER GAMMA
true_r0 = alphaX + betaX*gammaX
true_r0
model_params = c(alphaX, betaX, gammaX, true_r0)

#MCMC - sigma
sigma_a = 0.4*alphaX
sigma_b = 1.0*betaX 
sigma_g = 0.85*gammaX
sigma_bg = 1.5*gammaX
sigma = c(sigma_a, sigma_b, sigma_g, sigma_bg)
sigma_base = 0.25 #0.5

#####################################################******************************************************
#RUN INFERENCE: SS EVENTS
inference_type = 'ss_events_infer/'
results_home =  paste0(results_folder, inference_type)
print(results_home)

###############
# RUN I

#APPLY MCMC
model_type = 'sse_inf_sse_sim' #base_sim_sse_inf' #'ssi_sim_sse_inf'
flags_data_type = c(TRUE, FALSE, FALSE) #1)ss_events, 2) ss_individuals, 3) basline
base_folder_current = paste0(results_home, model_type, '/iter_', iter)
print(base_folder_current)
 
#START MCMC
start_time = Sys.time()
print('Start time:')
print(start_time)
run_mcmc_reps_ss(n_mcmc, n_reps, model_params, sigma, flags_data_type, base_folder_current, burn_in)
end_time = Sys.time()
timei = get_time(start_time, end_time)


###############
#APPLY SUMMARY STATS + p vals
model_type = 'sse_inf_sse_sim'
start_time = Sys.time()
print(paste0('Time elapsed:'), start_time)
get_sum_stats_total(base_folder_current, thinning_factor, n_reps, n_mcmc) 
df_p_valuesI = get_p_values_total(base_folder_current, n_reps) 
end_time = Sys.time()
timeI = get_timeII(start_time, end_time, timei)

#PLOT
plot_p_vals(df_p_valuesI)

############################################################
#RUN II - sse_inf_ssi_sim

#APPLY MCMC
model_type = 'sse_inf_ssi_sim' #'sse_inf_sse_sim' 'sse_inf_base_sim'
flags_data_type = c(FALSE, TRUE, FALSE) #1)ss_events, 2) ss_individuals, 3) basline
base_folder_current = paste0(results_home, model_type, '/iter_', iter)
print(base_folder_current)

#Repitions
start_time = Sys.time()
print('Start time:')
print(start_time)
run_mcmc_reps_ss(n_mcmc, n_reps, model_params, sigma, flags_data_type, base_folder_current, burn_in)
end_time = Sys.time()
timeii = get_time(start_time, end_time)


###############
#APPLY SUMMARY STATS + p vals
model_type = 'sse_inf_ssi_sim'
start_time = Sys.time()
print(paste0('Start time:', start_time))
get_sum_stats_total(base_folder_current, thinning_factor, n_reps, n_mcmc) 
df_p_valuesII = get_p_values_total(base_folder_current, n_reps)
end_time = Sys.time()
timeII = get_timeII(start_time, end_time, timeii)

#PLOT
plot_p_vals(df_p_valuesII)

############################################################
#RUN III - sse_inf_base_sim

#APPLY MCMC
model_type = 'sse_inf_base_sim' #'sse_inf_sse_sim' 'sse_inf_base_sim'
flags_data_type = c(FALSE, FALSE, TRUE) #1)ss_events, 2) ss_individuals, 3) basline
base_folder_current = paste0(results_home, model_type, '/iter_', iter)
print(base_folder_current)

#RUN MCMC
start_time = Sys.time()
print('Start time:')
print(start_time)
run_mcmc_reps_ss(n_mcmc, n_reps, model_params, sigma, flags_data_type, base_folder_current, burn_in)
end_time = Sys.time()
timeiii = get_time(start_time, end_time)

###############
#APPLY SUMMARY STATS + p vals
model_type = 'sse_inf_base_sim'
start_time = Sys.time()
print(paste0('Start time:', start_time))
get_sum_stats_total(base_folder_current, thinning_factor, n_reps, n_mcmc) 
df_p_valuesIII = get_p_values_total(base_folder_current, n_reps) 
end_time = Sys.time()
print(end_time)
timeIII = get_timeII(start_time, end_time, timeiii)

#PLOT
plot_p_vals(df_p_valuesIII)


#*******************************************************************************************************
#*
#*
#*******************************************************************************************************
#RUN INFERENCE: BASE MODEL
inference_type = 'base_infer/'
results_home =  paste0(results_folder, inference_type)
print(results_home)

#########################################################
# RUN I - SSE/BASE
model_type = 'base_inf_sse_sim' #base_inf_base_sim' #'base_inf_ssi_sim'
flags_data_type = c(TRUE, FALSE, FALSE) #1)ss_events, 2) ss_individuals, 3) basline
iter = 1
base_folder_current = paste0(results_home, model_type, '/iter_', iter)
print(base_folder_current)

# APPLY MCMC
start_time = Sys.time()
print('Start time:')
print(start_time)
run_mcmc_base_reps(n_mcmc, n_reps, true_r0, sigma_base, flags_data_type, base_folder_current, burn_in)
end_time = Sys.time()
timeBi = get_time(start_time, end_time)

###############
#APPLY SUMMARY STATS + p vals
start_time = Sys.time()
print(paste0('Time elapsed:'), start_time)
get_sum_stats_base_total(base_folder_current, thinning_factor, n_reps, n_mcmc) 
df_p_valuesBI = get_p_values_total(base_folder_current, n_reps)
end_time = Sys.time()
timeBI = get_timeII(start_time, end_time, timeBi)

#PLOT
plot_p_vals(df_p_valuesBI)

############################################################
#RUN II - base_inf_base_sim

model_type = 'base_inf_ssi_sim' #'sse_inf_sse_sim' 'sse_inf_base_sim'
flags_data_type = c(FALSE, TRUE, FALSE) #1)ss_events, 2) ss_individuals, 3) basline
base_folder_current = paste0(results_home, model_type, '/iter_', iter)
print(base_folder_current)

#APPLY MCMC
start_time = Sys.time()
print('Start time:')
print(start_time)
run_mcmc_base_reps(n_mcmc, n_reps, true_r0, sigma_base, flags_data_type, base_folder_current, burn_in)
end_time = Sys.time()
timeBii = get_time(start_time, end_time)


###############
#APPLY SUMMARY STATS + p vals
start_time = Sys.time()
print(paste0('Start time:', start_time))
get_sum_stats_base_total(base_folder_current, thinning_factor, n_reps, n_mcmc) 
df_p_valuesBII = get_p_values_total(base_folder_current, n_reps)
end_time = Sys.time()
timeBII = get_timeII(start_time, end_time, timeBii)

#PLOT
plot_p_vals(df_p_valuesBII)

############################################################
#RUN III - base_inf_base_sim

model_type = 'base_inf_base_sim' #'sse_inf_sse_sim' 'sse_inf_base_sim'
flags_data_type = c(FALSE, FALSE, TRUE) #1)ss_events, 2) ss_individuals, 3) basline
base_folder_current = paste0(results_home, model_type, '/iter_', iter)
print(base_folder_current)


#APPLY MCMC
start_time = Sys.time()
print('Start time:')
print(start_time)
run_mcmc_base_reps(n_mcmc, n_reps, model_params, sigma_base, flags_data_type, base_folder_current, burn_in)
end_time = Sys.time()
timeBiii = get_time(start_time, end_time)

###############
#APPLY SUMMARY STATS + p vals
start_time = Sys.time()
print(paste0('Start time:', start_time))
get_sum_stats_base_total(base_folder_current, thinning_factor, n_reps, n_mcmc) 
df_p_valuesBIII = get_p_values_total(base_folder_current, n_reps) 
end_time = Sys.time()
print(end_time)
timeBIII = get_timeII(start_time, end_time, timeBiii)

#PLOT
plot_p_vals(df_p_valuesBIII)

#Test return r0
#r02 = getr0(base_folder_current, thinning_factor, n_reps, n_mcmc) 
