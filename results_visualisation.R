#Results - Dispaly results of MCMC

#Final MCMC mean vs R0 mean

#Results


ggplot(dat, aes(x = x1, y = resp, color = grp) ) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)