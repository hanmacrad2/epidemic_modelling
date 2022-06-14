#****************************************************************
#1. APPLY SSI MODEL MCMC + DATA AUGMENTATION
#****************************************************************

#SETUP 
setwd("~/GitHub/epidemic_modelling") 
source("epidemic_functions.R") 
source("plot_functions.R") 
source("helper_functions.R") 

#DATA SIMULATION PARAMS
num_days = 50
shape_g = 6; scale_g = 1 #Infectious pressure (lambda) - gamma params

#SSI specific (*TO DO: DESIGN OF EXPERIMENTS FOR PARAM COMBINATIONS)
model_typeX = 'SSI'; 
aX = 0.8; bX = 0.1; cX = 10 
true_r0 = aX + bX*cX
true_r0
model_params = list(m1 = aX, m2 = bX, m3 = cX, true_r0 = true_r0)

#MCMC PARAMS  
n_mcmc = 10000 #1000
#SIGMA
sigma_a = 0.4*aX; sigma_b = 1.0*bX #0.1 #SHOULD SIGMA BE DEFINED MORE RIGOROUS
sigma_c = 0.85*cX; sigma_bc = 1.5*cX
sigma = list(sigma_a = sigma_a, sigma_b = sigma_b, #Acc rate too big -> Make sigma bigger. 
             sigma_c = sigma_c, sigma_bc = sigma_bc) #Acc rate too small -> make sigma smaller
time_elap = 0
mcmc_inputs = list(n_mcmc = n_mcmc, sigma = sigma, 
                   model_params = model_params, x0 = 1, seed_count = seed_count)


#****************************************************************
#DATASET - GENERATED USING SIMULATION FUNCTIONS
#****************************************************************
seed_count = 3 #seed_count = seed_count + 1 #print(paste0('i mcmc = ', i))
set.seed(seed_count)
sim_data = simulation_super_spreaders(num_days, shape_g, scale_g, aX, bX, cX)

#PLOTS
par(mfrow=c(2,1))
non_ss = sim_data[[1]]
plot.ts(non_ss, ylab = 'Daily Infections count', main = 'Non Super-Spreaders' )
ss = sim_data[[2]]
plot.ts(ss, ylab = 'Daily Infections count', main = 'Super-Spreaders')

#Total
sim_dataX = non_ss + ss
plot.ts(sim_dataX, ylab = 'Daily Infections count', main = 'Total - Super Spreaders Model, Daily Infections count')

#****************************************************************
#I. APPLY MCMC SSI MODEL   
#***************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_out = MCMC_SSI(sim_data, mcmc_inputs = mcmc_inputs)
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_out$time_elap = time_elap

#PLOT RESULTS
PLOT_MCMC_GRID(sim_dataX, mcmc_ssi_out,
               mcmc_inputs = mcmc_inputs,
               FLAGS_LIST = list(DATA_AUG = FALSE, BC_TRANSFORM = TRUE,
                                 PRIOR = TRUE, JOINT = TRUE,
                                 B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                 RJMCMC = FALSE))

#****************************************************************
# II APPLY MCMC SSI MODEL + DATA AUGMENTATION  
#***************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))
mcmc_ssi_da = MCMC_SSI(sim_data, mcmc_inputs = mcmc_inputs)

end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_da$time_elap = time_elap

#PLOT RESULTS
PLOT_MCMC_GRID(sim_dataX, mcmc_ssi_da,
               mcmc_inputs = mcmc_inputs,
               FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                 PRIOR = TRUE, JOINT = TRUE,
                                 B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                 RJMCMC = FALSE))

#****************************************************************
# III APPLY MCMC SSI MODEL + NON-SS EXTREME CASE
#***************************************************************
n_mcmc = 1000
#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))

mcmc_ssi_da2 = MCMC_SSI(sim_data, mcmc_inputs = mcmc_inputs,
                           FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                             PRIOR = TRUE, JOINT = TRUE,
                                             B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                             FLAG_NS_DATA_AUG = TRUE, FLAG_SS_DATA_AUG = FALSE))
end_time = Sys.time()
time_elap = get_time(start_time, end_time)
mcmc_ssi_da2$time_elap = time_elap

#PLOT RESULTS
PLOT_MCMC_GRID(sim_dataX, mcmc_ssi_da2,
               mcmc_inputs = mcmc_inputs,
               FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                 PRIOR = TRUE,  JOINT = TRUE,
                                 B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,  RJMCMC = FALSE,
                                 FLAG_NS_DATA_AUG = TRUE, FLAG_SS_DATA_AUG = FALSE))

# plot_mcmc_grid(n_mcmc, sim_dataX, mcmc_params_da3, true_r0, time_elap, seed_count,
#                model_params,
#                model_type = model_typeX,
#                FLAG_G_PRIOR_B = TRUE, gam_priors_on_b = c(10, 1/100),
#                rjmcmc = RJMCMCX, data_aug = TRUE,
#                mod_par_names = c('a', 'b', 'c'))

#****************************************************************
# APPLY MCMC SSI MODEL + SS EXTREME CASE
#***************************************************************

#START MCMC
start_time = Sys.time()
print(paste0('start_time:', start_time))

mcmc_ssi_da3 = MCMC_SSI(sim_data, mcmc_inputs = mcmc_inputs,
                           FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                             PRIOR = TRUE, JOINT = TRUE,
                                             B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,  RJMCMC = FALSE,
                                             FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = TRUE))

end_time = Sys.time()
mcmc_ssi_da3$time_elap = time_elap

#PLOT RESULTS
PLOT_MCMC_GRID(sim_dataX, mcmc_ssi_da3,
               mcmc_inputs = mcmc_inputs,
               FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                 PRIOR = TRUE,  JOINT = TRUE,
                                 B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,  RJMCMC = FALSE,
                                 FLAG_NS_DATA_AUG = FALSE, FLAG_SS_DATA_AUG = TRUE))