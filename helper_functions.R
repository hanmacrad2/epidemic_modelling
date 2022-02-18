#R - helper functions
#library(data.table)

################
#TIME

print_time <- function(start_time, end_time){
  'Print difference between end & start time'
  
  time_elap = round(end_time - start_time, 2)
  print('Time elapsed:') 
  print(time_elap)
  
}

get_time <- function(start_time, end_time, show = TRUE){
  'Print difference between end & start time'
  
  time_elap = round(end_time - start_time, 2)
  if(show){
    print('Time elapsed:') 
    print(time_elap)
  }
  time_elap
}

get_timeII <- function(start_time, end_time, timeI){
  'Print difference between end & start time'
  
  time_elap = round(end_time - start_time, 2) + timeI
  time_elap = time_hours(time_elap)
  print('Time elapsed:') 
  print(time_elap)
  
}

time_hours <- function(time_secs){
  
  time_hours = as.numeric(time_secs, units = "hours")
  time_hours = round(time_hours, 2)
  
}

#TIME
#time_elapI = round(difftime(end_time, start_time, units='hours'), 2) #round(end_time - start_time, 2)
#round(difftime(timeI, units='hours'), 2) #round(end_time - start_time, 2)
#as.numeric(timeI, units = "hours")

################
# DATA

#Trim data for a histogram between the lb quartile and upper bound quartile
trim_q <- function(x, lb, ub){
  'Trim data for a histogram between the lb quartile and upper bound quartile'
  
  x[(x > quantile(x, lb)) & (x < quantile(x, ub))]
}

upper_quantile <- function(x, ub){
  'Trim data for a histogram less then the upper bound quartile'
  
  x[x < quantile(x, ub)]
}

trim_t <- function(x){
  'Trim to remove outliers.
  Tukey defined it at 1.5× the interquartile range above and below the first and third quartile and not the mean'
  x[(x > quantile(x, 0.25)-1.5*IQR(x)) & (x < quantile(x, 0.75)+1.5*IQR(x))]
}

#****************************************************
#* DATAFRAMES

#Rename list of cols
rename_cols <- function(df, list_old, list_new) {
  df_new = df
  setnames(df_new, old = list_old, new = list_new)
  #cat('new_cols', names(df_new))
  df_new 
}

#Rename colums
rename_col <- function(df_new, old_name, new_name){
  #Rename column
  names(df_new)[names(df_new) == old_name] <- new_name
}


#PRINT

# print(paste0('seed_Count = ', seed_count))
# print(paste0('alphaX = ', alphaX))
# print(paste0('a_mcmc_mean = ', a_mcmc_mean))
# print(paste0('betaX = ', betaX))
# print(paste0('b_mcmc_mean = ', b_mcmc_mean))
# print(paste0('gammaX = ', gammaX))
# print(paste0('g_mcmc_mean = ', g_mcmc_mean))
# print(paste0('true_r0 = ', true_r0))
# print(paste0('accept_rate_a = ', round(mcmc_params[[5]],2)))
# print(paste0('a_rte_b = ', round(mcmc_params[[6]], 2)))
# 
# print(paste0('a_rte_g = ', round(mcmc_params[[7]],2)))
# print(paste0('a_rte_b_g = ', round(mcmc_params[[8]],2)))
# print(paste0('a_rte_rj0 = ', round(mcmc_params[[9]],2)))
# print(paste0('a_rte_rj1 = ', round(mcmc_params[[10]],2)))
# print(paste0('base_pc = ', base_pc))
# print(paste0('bayes_factor = ', bayes_factor))
# print(paste0('total_time = ', total_time))

