#Old versions

#################
#2. SUMMARY STATS
get_summary_stats <- function(data, flag_create){
  
  'Calculate summary statisitcs of the simulated data'
  #Summary stats params
  start_half2 = (length(data)/2)+1
  stop_half2 = length(data)
  
  if (flag_create){
    
    #Df
    summary_stats_results = data.frame(
      sum_infects = sum(data),
      median_infect = median(data),
      max_infect = max(data),
      sd_infect = sd(data),
      
      infect_q_75 = quantile(data)[4][1][1],
      infect_q_87_5 = quantile(data, probs = seq(0, 1, 0.125))[8][1][1], #mean(quantile(data)[4][1][1], quantile(data)[5][1][1]),
      max_dif = max((diff(data))), #Change from absolute difference
      med_dif = median(diff(data)),
      dif_dif = diff(diff(data)),
      
      #mean_upper_dif = mean(c(quantile(abs(diff(data)))[4][1][1], quantile(abs(diff(data)))[5][1][1])),
      sum_1st_half  = sum(data[1:(length(data)/2)]), #sum(which(data < quantile(data)[3][1][1]))
      sum_2nd_half =  sum(data[start_half2:stop_half2]) #sum(which(data > quantile(data)[3][1][1]))
    )
    
  } else {
    #List
    summary_stats_results = list(sum(data), median(data), max(data),
                                 sd(data), quantile(data)[4][1][1], 
                                 mean(quantile(data)[4][1][1], quantile(data)[5][1][1]),
                                 max(abs(diff(data))), median(abs(diff(data))),
                                 mean(c(quantile(abs(diff(data)))[4][1][1], quantile(abs(diff(data)))[5][1][1])),
                                 sum(data[1:(length(data)/2)]),
                                 sum(data[start_half2:stop_half2]) #sum(which(data > quantile(data)[3][1][1]))
    )
  }
  
  summary_stats_results
  
}

#EXAMPLE
list_p_vals = sapply(1:ncol(df_summary_stats_rep), function(x) get_p_values(df_summary_stats_rep[,x], df_true_ss[,x]))

#Miscellaneous - Divide by point in time
vecX = c(1,3,4,0,1,3,5,11,24,50,10,25,70,60,85,100,125)

unlist(sapply(2:length(vecX), function(i) vecX[i-1]/vecX[i]))

max_diff = max(lapply(1:length(vecX), function(i) diff(vecX)/mean(c(vecX[i],vecX[i+1]))))
max_diff

max(lapply(1:length(vecX)-1, function(i) diff(vecX)/mean(c(vecX[i],vecX[i+1])))) #get_p_values(df_summary_stats_rep[], df_true_ss[,x])))
#
#
dif_vec = diff(vecX)
dif_vec/mean(c(vecX[i],vecX[i+1]))

diff(vecX)/mean(c(vecX[i],vecX[i+1]))
max(lapply(vecX, function(i) diff(vecX)/mean(c(vecX[i],vecX[i+1]))))

mean(c(vecX[1:length(vecX)-1], vecX[2:length(vecX)]))

#Roll apply
library(zoo)
a <- 1:9
rollapply(a, 2, mean, by = 1, align = "left", partial = TRUE)
help(rollapply)
