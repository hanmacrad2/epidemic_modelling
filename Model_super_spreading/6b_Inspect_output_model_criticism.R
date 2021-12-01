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
  colors_line <- rainbow(c(20:20+len_data+1))
  
  #Sim_data
  #i.Infections
  plot.ts(sim_data_rep, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste(rep, "Day Infts SS Evnts, SS model, ", "R0 = ", true_r0), #model_type
          col = colorsX[1],
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Columns
  for (i in c(1:len_data)){
    print(colnames(df_sum_stats)[i])
    cat('num zeros =', length(which(df_sum_stats[,i] == 0)))
    sum_stat_trim = upper_quartile(df_sum_stats[,i], 0.95)
    cat('  After trim num zeros = ', length(which(df_sum_stats[,i] == 0)))
    cat('\n')
    hist(df_sum_stats[,i], breaks = 100, #freq = FALSE, 
         #xlim = c(xmin, xmax),
         xlab = paste('', toupper(colnames(df_sum_stats)[i])),
         ylab = 'Num Samples',
         col = colorsX[i+1],
         main = paste('', toupper(colnames(df_sum_stats)[i]),', p value:', round(list_p_vals[i],2)),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    #abline(v = list_p_vals[i], col = colors_line[i], lwd = 2.5) #This should be true value, not p value
  }
  
}

#Apply
rep = 14
get_rep_results(results_home, model_type, rep, true_r0)

#1. SS EVENTS  - MAKE FUNCTION
rep_interest = 14
results_inspect = paste0(results_home, model_type, "/iter_", iter, "/rep_", rep_interest, '/')
results_inspect
setwd(results_inspect)
#Data
d14 <- readRDS('df_summary_stats_14.rds')
ss14 <- readRDS('ss_data.rds')
list_p_vals14 <- readRDS(paste0(results_inspect, 'list_p_vals_73.rds'))

#Data
rep_interest = 73
results_inspect = paste0(results_home, model_type, "/iter_", iter, "/rep_", rep_interest, '/')
results_inspect
setwd(results_inspect)
d73 <- readRDS('df_summary_stats_73.rds')
ss73 <- readRDS('ss_data.rds')
plot.ts(ss73)
#p vals
list_p_vals <- readRDS('list_p_vals_73.rds')








#OLD!!!!!!!!!!!!!
#Read in data
sim_dataX =  readRDS(paste0(results_inspect, '/ss_data.Rdata'))
sim_dataX

sim_dataX =  read.table(paste0(results_inspect, '/ss_data.Rdata'))
sim_dataX

#Summary stats
(load(paste0(results_inspect, '/df_summary_stats_8.Rdata')))

df_sum_stats =  read.table(paste0(results_inspect, '/df_summary_stats_8.Rdata'))
df_sum_stats


load(paste0(results_inspect, '/df_summary_stats_8.RData'))

df_sum_stats <- readRDS("stuff.RDS")

#WORKING
#Test saving
newdf3 = newdf1
save(newdf3, file = paste0(results_home, 'check1.RData'))
#Load
load(paste0(results_home, 'check1.RData'))
load('check1.RData')

#CHECK
#save(newdf3, file = paste0(results_home, 'check1.RData'))
#Load
load('df_summary_stats_8.RData')
#load('check1.RData')

#Option II - WORKING :d 
saveRDS(newdf3, 'newdf3.rds')
#remove(newdf3)
d3 <- readRDS('newdf3.rds')
