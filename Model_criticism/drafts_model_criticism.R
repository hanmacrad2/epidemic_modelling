#Old versions

#################
#2. SUMMARY STATS
get_summary_stats <- function(data, flag_create){
  
  'Calculate summary statisitcs of the simulated data'
  #Summary stats params
  start_half2 = (length(data)/2)+1
  stop_half2 = length(data)
  
  if (flag_create){
    
    #Df of summary stats (20)
    summary_stats_results = data.frame( 
      
      sum_infects = sum(data),
      sum_1st_half  = sum(data[1:(length(data)/2)]),
      sum_2nd_half =  sum(data[start_half2:stop_half2]),
      
      median_infect = median(data),
      max_infect = max(data),
      sd_infect = sd(data),
      
      infect_q_75 = quantile(data)[4][1][1],
      infect_q_87_5 = quantile(data, probs = seq(0, 1, 0.125))[8][1][1],
      
      #Differences
      max_dif = max((diff(data))), #Change from absolute difference
      med_dif = median(diff(data)),
      max_dif2nd = max(diff(diff(data))),
      med_dif2nd = median(diff(diff(data))),
      
      #Norm Differences
      max_dif_normI = max(diff(data)/data[1:length(data)-1]),
      max_dif_normII = max(diff(data)/rollapply(data, 2, mean, by = 1)),
      max_dif2nd_I = max(diff(diff(data))/data[1:length(data)-1]),
      max_dif_2ndII = max(diff(diff(data))/rollapply(data, 2, mean, by = 1)),
      
      med_dif_normI = median(diff(data)/data[1:length(data)-1]),
      med_dif_normII = median(diff(data)/rollapply(data, 2, mean, by = 1)),
      med_dif2nd_I = median((diff(diff(data))/data[1:length(data)-1])),
      med_dif_2ndII = median(diff(diff(data))/rollapply(data, 2, mean, by = 1))
      
    )
    
    
    
  } else {
    #List
    summary_stats_results = list(sum(data), sum(data[1:(length(data)/2)]), sum(data[start_half2:stop_half2]),
                                 median(data), max(data), sd(data),
                                 quantile(data)[4][1][1], quantile(data, probs = seq(0, 1, 0.125))[8][1][1],
                                 max(diff(data)), median(diff(data)), max(diff(diff(data))), median(diff(diff(data))),
                                 max(diff(data)/data[1:length(data)-1]), max(diff(data)/rollapply(data, 2, mean, by = 1)),
                                 max(diff(diff(data))/data[1:length(data)-1]),  max(diff(diff(data))/rollapply(data, 2, mean, by = 1)),
                                 median(diff(data)/data[1:length(data)-1]), median(diff(data)/rollapply(data, 2, mean, by = 1)),
                                 median((diff(diff(data))/data[1:length(data)-1])), median(diff(diff(data))/rollapply(data, 2, mean, by = 1))
                                 

    )
  }
  
  summary_stats_results
  
}

#Total p values 
get_p_values_total <- function(base_folder_current, n_reps){
  
  for(rep in 1:n_reps) {
    cat('rep = ', rep)
    #Get results
    folder_rep = paste0(base_folder_current, "/rep_", rep, '/')
    cat('folder_rep', folder_rep)
    true_rep_sim = readRDS(paste0(folder_rep, '/sim_data.rds'))
    #Get true summary statistics 
    df_true_ss = get_summary_stats(true_rep_sim, TRUE)
    
    #Data
    df_summary_stats_rep <- readRDS(paste0(folder_rep, '/df_summary_stats_', rep, '.rds' ))
    
    #Get p values
    list_p_vals = sapply(1:ncol(df_summary_stats_rep), function(x) get_p_values(df_summary_stats_rep[,x], df_true_ss[,x]))
    saveRDS(list_p_vals, file = paste0(folder_rep, '/list_p_vals_rep', rep, ".rds"))
    
    list_all_p_vals = sapply(1:ncol(df_summary_stats_rep), function(x) get_p_values_list(df_summary_stats_rep[,x], df_true_ss[,x]))
    saveRDS(list_all_p_vals, file = paste0(folder_rep, '/list_all_p_vals_rep_', rep, ".rds"))
    
    #Save all 
    if (!exists("df_p_values")) {
      df_p_values = data.frame(sum_infects = list_p_vals[1],
                               sum_1st_half = list_p_vals[2],
                               sum_2nd_half = list_p_vals[3],
                               median_infect = list_p_vals[4],
                               max_infect = list_p_vals[5],
                               sd_infect = list_p_vals[6],
                               infect_q_75 = list_p_vals[7],
                               infect_q_87_5 = list_p_vals[8],
                               max_dif = list_p_vals[9],
                               med_dif  = list_p_vals[10],
                               max_dif2nd =  list_p_vals[11],
                               med_dif2nd =  list_p_vals[12],
                               max_dif_normI =  list_p_vals[13],
                               max_dif_normII =  list_p_vals[14],
                               max_dif2nd_I =  list_p_vals[15],
                               max_dif_2ndII =  list_p_vals[16],
                               med_dif_normI =  list_p_vals[17],
                               med_dif_normII =  list_p_vals[18],
                               med_dif2nd_I =  list_p_vals[19],
                               med_dif_2ndII =  list_p_vals[20]
                               
      )
      print(paste0('df_p_values', df_p_values))
      
    } else {
      df_p_values[nrow(df_p_values) + 1, ] = list_p_vals
    }
    
  }
  
  #Ensure its a df
  df_p_values = as.data.frame(df_p_values)
  
  #SaveRDS
  saveRDS(df_p_values, file = paste0(base_folder_current, '/total_p_values_iter_', iter, '.rds' ))
  
  #Return p values
  df_p_values
  
}

#EXAMPLE
list_p_vals = sapply(1:ncol(df_summary_stats_rep), function(x) get_p_values(df_summary_stats_rep[,x], df_true_ss[,x]))

#Miscellaneous - Divide by point in time
vecX = c(1,3,4,0,0,1,3,5,11,24,50,10,25,70,60,85,100,125)

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


rollapply(a, 2, mean, by = 1)

round(diff(a)/rollapply(a, 2, mean, by = 1),2)

#DIf/mean
round(diff(vecX)/rollapply(vecX, 2, mean, by = 1),2)

#Dif/num 1
k = round(diff(vecX)/vecX[1:length(vecX)-1], 3)

#******************************************************************************
#OPTION 1 dif_normI
dif_normI = round(diff(a)/a[1:length(a)-1], 3)

#OPTION 2 dif_normI
dif_normII = round(diff(vecX)/(rollapply(vecX, 2, mean, by = 1)+1),2)
dif_normII

#dif2nd_II
#dif2nd_I

#OPTION III Ratio
dif_normIII = round(diff(a)/((a[1:length(a)-1]/a[2:length(a)])+1), 3)
dif_normIII

dif_normIII = round(diff(a)/(rollapply(a, 2, FUN = function(x, y) x/y, by = 1)+1),2)
dif_normIII
