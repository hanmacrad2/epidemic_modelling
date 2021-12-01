#INSPECT OUTPUT OF MODEL CRITICISM

#Set up params
results_home = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism/"

#FUNCTIONS - GET DATA
get_rep_results <- function(results_home, model_type, rep, true_r0){
  
  #Results inspect
  results_inspect = paste0(results_home, model_type, "/iter_", iter, "/rep_", rep, '/')
  print(results_inspect)
  
  #Data
  sim_data_rep <- readRDS(paste0(results_inspect, 'ss_data.rds'))
  cat('Sum sim data = ', sum(sim_data_rep))
  df_sum_stats <- readRDS(paste0(results_inspect, 'df_summary_stats_', rep, '.rds'))
  list_p_vals <- readRDS(paste0(results_inspect, 'list_p_vals_', rep, '.rds'))
  
  #Plot results
  plot_rep_results(true_r0, model_type, sim_data_rep, df_sum_stats, list_p_vals) 
  
}

#Plot results 
plot_rep_results <- function(true_r0, model_type, sim_data_rep, df_sum_stats, list_p_vals){
  
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
    sum_stat_trim = upper_quantile(X, 0.99)

    hist(sum_stat_trim, breaks = 100, #freq = FALSE, 
         #xlim = c(xmin, xmax),
         xlab = paste('', toupper(colnames(df_sum_stats)[i]), '< 99th quantile'),
         ylab = 'Num Samples',
         col = colorsX[i+1],
         main = paste('', toupper(colnames(df_sum_stats)[i]),', p value:', round(list_p_vals[i],3)),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(v = df_sum_stats[nrow(df_sum_stats),i], col = colors_line[i], lwd = 2.5) #This should be true value, not p value
  }
  
}


#Apply
rep = 7 #14, 73
get_rep_results(results_home, model_type, rep, true_r0)

