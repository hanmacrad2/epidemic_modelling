#INSPECT OUTPUT OF MODEL CRITICISM

#Set up params
setwd("~/GitHub/epidemic_modelling")
source("epidemic_functions.R")
source("helper_functions.R")
source("Model_criticism/7A_Model_Criticism_SS_Events.R")
source("Model_criticism/7B_Run_model_criticism_all_sse_and_base.R")
results_folder =  "~/PhD_Warwick/Project_Epidemic_Modelling/Results/model_criticism/model_criticism_1k_I/"

#RESULTS FOLDER: SS EVENTS
inference_type = 'ss_events_infer/'
results_home =  paste0(results_folder, inference_type)
print(results_home)
iter = 1

#*##########################################################
#1. GET & DISPLAY TOTAL REP RESULTS 
display_rep_results <- function(results_home, model_type, iter, rep, n_mcmc, true_r0,
                            upper_quant, trim_flag, list_i, time_elap, flags_data_type){
  #P values
  df_p_vals_tot = get_df_p_vals(results_home, model_type, iter)
  #plot_p_vals(df_p_vals_tot, model_type)
  
  #Specific rep results inspect
  results_inspect = paste0(results_home, model_type, "/iter_", iter, "/rep_", rep, '/')
  print(results_inspect)
  
  #Data
  sim_data_rep <- readRDS(paste0(results_inspect, 'sim_data.rds')) 
  print(paste0('Sim data rep ', rep, '='))
  print(sim_data_rep)
  mcmc_params <- readRDS(paste0(results_inspect, '/mcmc_params_rep_', rep, '.rds' ))
  df_sum_stats <- readRDS(paste0(results_inspect, 'df_summary_stats_', rep, '.rds'))
  df_true_sum_stats <- readRDS(paste0(results_inspect, 'df_true_sum_stats_rep_', rep, '.rds' ))
  list_p_vals <- readRDS(paste0(results_inspect, 'list_p_vals_rep', rep, '.rds'))
  print(df_true_sum_stats)
  
  #Plot simualated results (random sample)
  plot_sim_data_mcmc_runs(results_home, sim_data_rep, mcmc_params, list_i, flags_data_type)
  
  #Plot MCMC results 
  plot_mcmc_grid(n_mcmc, sim_data_rep, mcmc_params, true_r0, 'Neg Bin,', time_elap, rep, TRUE, TRUE)
  
  #Plot p vals & summary stats
  plot_rep_sum_stats(true_r0, model_type, sim_data_rep, df_sum_stats, df_true_sum_stats, list_p_vals, upper_quant, trim_flag) 
  print(list_p_vals)
}

#################
#P VALUE DATASET 
get_df_p_vals <- function(results_home, model_type, iter){
  
  #Results inspect
  results_inspect = paste0(results_home, model_type, "/iter_", iter)
  print(results_inspect)
  
  for(rep in 1:n_reps) {
    
    results_rep = paste0(results_inspect, "/rep_", rep, '/')
    
    list_p_vals <- readRDS(paste0(results_rep, 'list_p_vals_rep', rep, '.rds'))
    
    if (rep == 10){
      print(list_p_vals)
    }
    
    if (rep == 1) {
      df_p_values = data.frame(sum_inf_counts = list_p_vals[1],
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
      print(paste0('df_p_values', df_p_values))
      
    } else {
      df_p_values[nrow(df_p_values) + 1, ] = list_p_vals
    }
  }
  
  #Return p values
  df_p_values
}

#df_p_vals_si = get_df_p_vals(results_home, model_type, iter)

###################
#PLOT P VALUES
plot_p_vals <- function(df_p_vals, model_type){
  
  'Plot histograms of the p values'
  par(mfrow=c(4,5)) #c(3,4)
  
  #Prop lt than 0.05
  num_iters = length(df_p_vals[,1])
  
  for (i in c(1:20)){ #11
    
    #Prop lt 0.05
    val_05 = 0.05
    percent_lt_05 = (length(which(df_p_vals[,i] < val_05))/num_iters)*100
    
    if (i == 1){
      print(paste0('i = ', i))
      hist(df_p_vals[,i], breaks = 100, #freq = FALSE, 
           #xlim = c(xmin, xmax),
           xlab = paste0('p value, < 0.05: ', percent_lt_05, '%'),
           ylab = 'Num Samples', col = 'green',
           main = paste('', toupper(colnames(df_p_vals)[i]),'R0:', true_r0, '', model_type),
           cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      abline(v = 0.05, col = 'red', lwd = 2)
      
    } else {
      hist(df_p_vals[,i], breaks = 100, #freq = FALSE, 
           #xlim = c(xmin, xmax),
           xlab = paste0('p value, < 0.05: ', percent_lt_05, '%'),
           ylab = 'Num Samples', col = 'green',
           main = paste('', toupper(colnames(df_p_vals)[i]), ', n reps:', num_iters),
           cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      abline(v = 0.05, col = 'red', lwd = 2)
      
    }
  }
}

##################
#PLOT SUMMARY STATS
plot_rep_sum_stats <- function(true_r0, model_type, sim_data_rep, df_sum_stats, df_true_sum_stats,
                               list_p_vals, upper_quant, trim_flag){
  
  'Plot sim data, summary stats and true summary stat for a given mcmc rep' 
  #Setup
  par(mfrow = c(5, 5)) #c(3,4))
  len_data = length(list_p_vals)
  print(paste0('length of data:', len_data))
  colorsX <- rainbow(len_data+1)
  colors_line <- rainbow(c(15:15+len_data+1))
  
  #i. Sim_data - Infections
  plot.ts(sim_data_rep, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste(rep, "Infts,", model_type, "R0 = ", true_r0), #model_type
          col = colorsX[1],
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Columns
  for (i in c(1:len_data)){
    print(paste0('i = ', i))
    X = df_sum_stats[,i] # df_sum_stats[1:nrow(df_sum_stats),i]
    
    if (trim_flag){
      print('trimmed')
      X = upper_quantile(X, upper_quant)
    }
    
    # print(paste0('Col i = ', colnames(df_sum_stats)[i]))
    # print(paste0('True val i = ', round(df_true_sum_stats[nrow(df_true_sum_stats), i], 2)))
    # print(paste0('X = ', X))
    
    #Histogram
    hist(X, breaks = 100, #freq = FALSE, 
         xlab = paste0('', toupper(colnames(df_sum_stats)[i]), ', T:', round(df_true_sum_stats[nrow(df_true_sum_stats), i], 2)),
         ylab = 'Num Samples',
         col = colorsX[i+1],
         main = paste('', toupper(colnames(df_sum_stats)[i]),', p value:', round(list_p_vals[i],3)),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(v = df_true_sum_stats[nrow(df_true_sum_stats), i], col = colors_line[i], lwd = 2.5) #This should be true value, not p value

  }
  
}

###################################
#INSPECT MCMC
plot_sim_data_mcmc_runs <- function(results_home, sim_data, mcmc_params, list_idx, flags_data_type){
  
  #Results inspect
  sim_data_path = paste0(results_home, model_type, "/iter_", iter, "/rep_", rep, '/mcmc/')
  print(sim_data_path)
  #Flags
  
  #Plot
  colorsX <- rainbow(length(list_idx)+1)
  par(mfrow=c(3,4))
  
  #i. Sim infections data (True for the rep)
  
  #SSE or SSI
  if(!flags_data_type[3]) {
    plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
            main = paste0(rep, ", ", model_type, ", True a:", alphaX, ', b:', betaX, ', g:', gammaX, ', R0:', true_r0), #model_type
            col = colorsX[1],
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  } else if (flags_data_type[3]){
    plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
            main = paste0(rep, ", ", model_type, ", True R0:", true_r0), #model_type
            col = colorsX[1],
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    
  }

  
  #Loop - Get mcmc rep
  for (i in seq_along(list_idx)){
    get_sim_data_i(list_idx[i], mcmc_params, sim_data_path, colorsX[i])
  }
  
}

#Plot mcmc iteration data
get_sim_data_i <- function(i, mcmc_params, sim_data_path, colourX){
  
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

#GET MCMC DATA
get_mcmc_results <- function(results_home, model_type, iter, rep, true_r0, time_elap){
  
  'Get MCMC Results for the given rep'
  
  #Results inspect
  results_inspect = paste0(results_home, model_type, "/iter_", iter, "/rep_", rep, '/')
  print(results_inspect)
  
  #Data
  sim_data_rep <- readRDS(paste0(results_inspect, 'sim_data.rds')) 
  cat('Sum sim data = ', sum(sim_data_rep))
  mcmc_params <- readRDS(paste0(results_inspect, 'mcmc_params_rep_', rep, '.rds' ))
  
  #Plot MCMC results 
  plot_mcmc_grid(sim_data_rep, mcmc_params, true_r0, 'Neg Bin,', time_elap, rep, TRUE, TRUE)
  plot_sim_data_mcmc_runs(results_home, sim_data_rep, mcmc_params, list_i)
  
}

##############################################
#MODELs x3 APPLY - INSPECT SPECIFIC REPS
model_type = 'sse_inf_sse_sim' #sse_inf_base_sim' #sse_inf_ssi_sim' #sse_inf_sse_sim' #'' #  '' 
flags_data_type = c(TRUE, FALSE, FALSE)

#Time
time_elap = timeIII
time_elap = time_hours(time_elap)
print(paste0('timeI: ', time_elap))

#SSE
df_baseI = get_df_p_vals(results_home, model_type, iter)
plot_p_vals(df_baseI)
plot_p_vals(df_p_valuesIII)

#Rep specific
upper_quant = 0.99 #8 #8 #6 #7 #8 #0.99 #1.0
trim_flag = FALSE #TRUE #FALSE # # #
list_i = seq(from = 500, to = 5500, by = 500)

#Rep
rep = 6 #34 #10 #3 #12 #10 #33 #10 #16 12 #98 #98 #7 #72 ##38 #40 #34 #40 #39 #31 #28 #27 #14 #26 #19 #7 #1 #34 #26 #54
display_rep_results(results_home, model_type, iter, rep, n_mcmc, true_r0,
                upper_quant, trim_flag, list_i, time_elap, flags_data_type)

