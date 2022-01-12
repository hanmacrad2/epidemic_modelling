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

get_time <- function(start_time, end_time){
  'Print difference between end & start time'
  
  time_elap = round(end_time - start_time, 2)
  print('Time elapsed:') 
  print(time_elap)
  
}

get_timeII <- function(start_time, end_time, timeI){
  'Print difference between end & start time'
  
  time_elap = round(end_time - start_time, 2) + timeI
  print('Time elapsed:') 
  print(time_elap)
  
}

#TIME
#time_elapI = round(difftime(end_time, start_time, units='hours'), 2) #round(end_time - start_time, 2)

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
