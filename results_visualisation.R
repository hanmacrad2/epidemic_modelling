#Results - Dispaly results of MCMC

#Final MCMC mean vs R0 mean

#Results
vec_r0 = c(0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0)
vec_final_mean = c(0.58, 0.52, 1.01, 0.82, 1.1, 1.32, 0.925, 2.66, 2.89, 3.39, 2.6, 3.65, 3.41, 4.31, 3.725, 4.66, 7.5, 9.06) 
df_accuracy = data.frame(vec_r0, vec_final_mean)
group = c(rep.int(1, length(vec_r0)), rep.int(2, length(vec_r0)))

ggplot(df_accuracy, aes(x = vec_r0, y = resp, color = grp) ) + #3.5 == 3.41, 3 = 3.65, 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

#
vec_final_mean = c(0.58, 0.52, 1.01, 0.82, 1.1, 1.32, 0.925, 2.66, 2.89, 2.6, 3.63, 3.39, 2.6, 3.625, 3.41, 4.31, 4.66, 7.5, 9.06) 
