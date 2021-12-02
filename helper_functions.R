#R - helper functions
library(data.table)

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
  df_new 
}

#Rename colums
rename_col <- function(df_new, old_name, new_name){
  #Rename column
  names(df_new)[names(df_new) == old_name] <- new_name
}