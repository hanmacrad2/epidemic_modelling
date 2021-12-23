#MODEL CRITICISM - BASE MODEL

#Setup
source("functions.R")
source("7A_Model_Criticism_SS_Events.R")
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
      sim_data = simulation_super_spreaders(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
      print('simulate ss events')
    } else if (flag2){
      sim_data = simulation_super_spreaders(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
      print('simulate ss individalss')
      saveRDS(sim_data, file = paste0(folder_rep, '/sim_data.rds'))
    } else if (flag3) {
      sim_data = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
      saveRDS(sim_data, file = paste0(folder_rep, '/sim_data.rds'))
      print('simulate base model')
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
