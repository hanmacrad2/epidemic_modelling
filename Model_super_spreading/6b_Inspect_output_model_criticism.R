#INSPECT OUTPUT OF MODEL CRITICISM

#Look at specific iterations
results_home = "~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism/"

#1. SS EVENTS  - MAKE FUNCTION
rep_interest = 14
results_inspect = paste0(results_home, model_type, "/iter_", iter, "/rep_", rep_interest, '/')
results_inspect
setwd(results_inspect)
#Data
d14 <- readRDS('df_summary_stats_14.rds')
ss14 <- readRDS('ss_data.rds')

#Data
rep_interest = 73
results_inspect = paste0(results_home, model_type, "/iter_", iter, "/rep_", rep_interest, '/')
results_inspect
setwd(results_inspect)
d73 <- readRDS('df_summary_stats_73.rds')
ss73 <- readRDS('ss_data.rds')





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
