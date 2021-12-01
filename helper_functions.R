#R - helper functions

#Trim data for a histogram between the lb quartile and upper bound quartile
trim_q <- function(x, lb, ub){
  'Trim data for a histogram between the lb quartile and upper bound quartile'
  
  x[(x > quantile(x, lb)) & (x < quantile(x, ub))]
}

upper_quartile <- function(x, ub){
  'Trim data for a histogram less then the upper bound quartile'
  
  x[x < quantile(x, ub)]
}