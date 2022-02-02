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


#Loops
for (t in 2:num_days) {
  
  if(sim_data[t] == 0){ 
    print('inf = 0')
  } else {
    inner_sum_vec <- vector('numeric', sim_data[t])
    
    for (y_t in 0:sim_data[t]){
      print(paste0('yt = ', y_t))
  }
  
  }
}

##############################
#1. MODEL COMPARISON
##############################

rjmcmc_sse_base <- function(data, n, sigma, x0 = 1, prior = TRUE) { #thinning_factor, burn_in
  
  'Returns mcmc samples of alpha & acceptance rate'
  print('MCMC SUPERSPREADING')
  
  #Initialise params
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  alpha_vec[1] <- x0; beta_vec[1] <- x0;
  gamma_vec[1] <- x0; r0_vec[1] <- x0;
  
  #Extract params
  sigma_a = sigma[1]; sigma_b = sigma[2]
  sigma_g = sigma[3]; sigma_bg = sigma[4];
  
  #Result vectors
  count_accept1 = 0; count_accept2 = 0;
  count_accept3 = 0; count_accept4 = 0; 
  count_accept5 = 0; count_accept6 = 0;
  #flag_true = FALSE
  
  #Create folder for mcmc results 
  #folder_mcmc = paste0(folder_results, '/mcmc')
  #ifelse(!dir.exists(file.path(folder_mcmc)), dir.create(file.path(folder_mcmc), recursive = TRUE), FALSE)
  
  #MCMC chain
  for(i in 2:n) { #1:n
    
    #******************************************************
    #ALPHA
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma_a) 
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    #log alpha
    logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev  #+ prior1 - prior
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - alpha_dash + alpha_vec[i-1]
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
    }
    
    #************************************************************************
    #BETA (Only if B > 0)
    if (beta_vec[i-1] > 0){ #& (i > 1)
      
      beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma_b) 
      if(beta_dash < 0){
        beta_dash = abs(beta_dash)
      }
      #loglikelihood
      logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1])
      logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
      log_accept_prob = logl_new - logl_prev
      #Priors
      if (prior){
        log_accept_prob = log_accept_prob - beta_dash + beta_vec[i-1]
      }
      #Metropolis Acceptance Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        count_accept2 = count_accept2 + 1
      } else {
        beta_vec[i] <- beta_vec[i-1]
      }
      
      #************************************************************************
      #GAMMA
      gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma_g) 
      if(gamma_dash < 1){
        gamma_dash = 2 - gamma_dash #abs(gamma_dash)
      }
      #Acceptance Probability
      logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash) 
      logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
      log_accept_prob = logl_new - logl_prev 
      #Priors
      if (prior){
        log_accept_prob = log_accept_prob - gamma_dash + gamma_vec[i-1]
      }
      #Metropolis Acceptance Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        gamma_vec[i] <- gamma_dash
        count_accept3 = count_accept3 + 1
      } else {
        gamma_vec[i] <- gamma_vec[i-1]
      }
      
      #R0
      r0_vec[i] = alpha_vec[i] + beta_vec[i]*gamma_vec[i]
      
      #*****************************************************
      #GAMMA-BETA
      gamma_dash <- gamma_vec[i] + rnorm(1, sd = sigma_bg)#Alter sigma_bg depending on acceptance rate. 
      #Acc rate too big -> Make sigma bigger. Acc rate too small -> make sigma smaller
      if(gamma_dash < 1){ #If less then 1
        gamma_dash = 2 - gamma_dash #abs(gamma_dash)
      }
      #New beta 
      beta_new = (r0_vec[i] - alpha_vec[i])/gamma_dash #Proposing new Gamma AND Beta. Beta_dash = f(R0 & gamma_dash)
      
      if(beta_new >= 0){ #Only accept values of beta > 0
        
        logl_new = log_like_ss_lse(data, alpha_vec[i], beta_new, gamma_dash)
        logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i])
        log_accept_prob = logl_new - logl_prev  
        #Priors
        if (prior){
          log_accept_prob = log_accept_prob - beta_dash + beta_vec[i]
        }
        #Metropolis Step
        if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
          beta_vec[i] <- beta_new
          count_accept4 = count_accept4 + 1
        } 
      }
    }
    
    #RJMCMC Step 
    #Reverse of proposals. Prob of proposing 0 when not 0, = 1. Therefore one of qs is 1. 
    if ((beta_vec[i-1] > 0) | (gamma_vec[i-1] > 0)){ #Proposal. 
      print('B 0 proposal')
      #Not a random walk metropolis - as not using current values to decide the next. As current values are zero - choosing non zero.
      #Independence Sampler 
      beta_dash = 0
      gamma_dash = 0
      #Everything cancels
      logl_new = log_like_ss_lse_B0(data, alpha_vec[i], beta_dash, gamma_dash)
      print(paste0('B0 log_l = ', logl_new))
      logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
      print(paste0('Prev log_l = ', logl_prev))
      log_accept_prob = logl_new - logl_prev 
      
      #Metropolis Step
      unif_var = runif(1)
      print(paste0('unif_var = ', unif_var))
      print(paste0('log_accept_prob = logl_new - logl_prev:', log_accept_prob))
      
      if(!(is.na(log_accept_prob)) && log(unif_var) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        gamma_vec[i] <- gamma_dash
        count_accept5 = count_accept5 + 1
      } 
      
    } else { #This acceptance prob will be the reverse of the first version
      print('B Independ. proposal')
      
      #Independence sampler - Propose from prior. If VERY lucky value is accepted to be able to jump between models. 
      beta_dash = rexp(1) #q - the proposal distribution is equal to the prior disribution. Reason: Acc prob = like*prior*q/(like*prior*q)
      gamma_dash = rexp(1) + 1
      
      #Everything cancels
      logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_dash)
      logl_prev = log_like_ss_lse_B0(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
      log_accept_prob = logl_new - logl_prev 
      
      unif_var = runif(1)
      print(paste0('unif_var = ', unif_var))
      print(paste0('log_accept_prob = ', log_accept_prob))
      
      #Metropolis Step
      if(!(is.na(log_accept_prob)) && log(unif_var) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        gamma_vec[i] <- gamma_dash
        count_accept6 = count_accept6 + 1
      } 
    }
    
  }
  
  #Final stats
  accept_rate1 = 100*count_accept1/n
  accept_rate2 = 100*count_accept2/n
  accept_rate3 = 100*count_accept3/n
  accept_rate4 = 100*count_accept4/n
  accept_rate5 = 100*count_accept5/n
  accept_rate6 = 100*count_accept6/n
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3, accept_rate4, accept_rate5, accept_rate6))
}