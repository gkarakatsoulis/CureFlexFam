library(dplyr)

myfunction = function(N, lambda, exp_rate, logit_par, weib_par, seed){
  
  # We start with determining the population ------------------------------
  set.seed(seed)
  
  
  # -------------------------------------------------------------------
  #                       1. Model Covariates generation
  # -------------------------------------------------------------------
  
  
  z1 = runif(N, -1, 1)
  z2 = rbinom(n = N, size = 1, prob = 0.5)
  
  
  
  # -------------------------------------------------------------------
  #                       2. Logistic regression part
  # -------------------------------------------------------------------
  
  
  
  # a. Generate the logistic regression
  logit = 0.5 + logit_par[1] * z1 + logit_par[2] * z2
  
  
  # b. Calculate the probabilities
  prob = plogis(logit)
  
  
  # c. Generate the vector of the cure status (Y_i: 1 if uncured, 0 if cured)
  y_i = rbinom(n = N, size = 1, prob = prob)
  
  
  
  
  # -------------------------------------------------------------------
  #                       3. Censoring Mechanism
  # -------------------------------------------------------------------
  
  
  # a. Generate the censoring times (coming from exponential distribution)
  cens_time_i = rexp(n = N, rate = exp_rate)
  
  
  
  # -------------------------------------------------------------------
  #                       4. Susceptible times to event
  # -------------------------------------------------------------------
  
  # We will use the inversion method.
  
  # a. Generate a uniform random sample (to use it for the S_u)
  u_cdf = runif(N, 0, 1)
  
  # b. Transform u_cdf to be for the S_p cdf (which we will inverse)
  r = prob + (1-prob) * (1-u_cdf)
  
  
  # c. Use the cure probabilities to derive the first formula (=F(t))
  x = (1-r^lambda) / (1-prob^lambda)
  
  # d. Define the Weibull regression parameters
  a1 = weib_par[1]
  a2 = weib_par[2]
  
  # e. Calculate the S_u random sample times using the inverse formula
  t_i = (-log(1-x))^(1/a2) / a1
  
  
  # -------------------------------------------------------------------
  #                       5. Merge in a dataset
  # -------------------------------------------------------------------
  
  df = data.frame('Cure_status' = y_i,
                  'z1' = z1,
                  'z2' = z2,
                  'Censoring_time' = cens_time_i,
                  'Event_time' = t_i)
  
  df$Observed_time = apply(df[,c('Censoring_time', 'Event_time')], 1, min)
  
  df$Censor_status = ifelse(df$Censoring_time < df$Event_time | df$Cure_status == 0, 0, 1)

  
  file_name = paste0('_Censor_Exp(', round(exp_rate,2), ')',
                     '_logit_(', paste0(round(exp_rate,2), collapse = '_'),
                     ')_weibul_(', paste0(round(exp_rate,2), collapse = '_'),
                     ')_lambda_', round(exp_rate,2))
  
  myfile_save = gsub('Rscripts', paste0('Simulated_Datasets/', file_name,'//'), getwd())
  
  ifelse(!dir.exists(myfile_save), dir.create(myfile_save), FALSE)
  
  
  openxlsx::write.xlsx(df, paste0(myfile_save, 'Dataset_', seed, '.xlsx'))
  
  return(df)
  
  
}

myseed = 1:100

for (j in myseed){
  
  myfunction(N = 500,
             exp_rate = 1/3,
             lambda = 1,
             seed = j,
             logit_par = c(0.5, 1),
             weib_par = c(0.8, 1))
  
}

