# This R script only optimizes the loglikelihood over the incidence and prevalence parts (glm and Weibull)

library(dplyr)

myfile_wd = gsub('Rscripts', 'Simulated_Datasets/', getwd())

myfolders = list.files(myfile_wd)

mydatasets = list.files(paste0(myfile_wd, myfolders[1]))

df = list()

for (i in mydatasets){
  
  df[[i]] = readxl::read_xlsx(paste0(myfile_wd, myfolders[1], '/', i))
  
}



# -------------------------------------------------------------
# Derive the functions for the log likelihood
# ------------------------------------------------------------- 




# ---------------------------
# -------------------------------------------------------------
# Derive the functions for the incidence part
# ------------------------------------------------------------- 
# ---------------------------

pi_fun = function(X, beta){ # Function for the cure probabilities
  
  z = 1/( 1 + exp(-X %*% beta) ) |> as.vector()
  
  return(z)
  
}


loglike_glm = function(parameter, X, y){ # Loglikelihood
  
  
  beta = parameter[1:3]
  
  pi = pi_fun(X = X, beta = beta)
  
  z = (1-y) * log(1-pi) + y * log(pi)
  
  result = sum(z)
  
  return(result)
  
}


# ---------------------------
# -------------------------------------------------------------
# Derive the functions for the latency part
# ------------------------------------------------------------- 
# ---------------------------


# Population Survival Function ---------------------------------
sp = function(X, beta, gamma, t, alpha){
  
  a1 = alpha[1]
  a2 = alpha[2]
  
  pi = pi_fun(X, beta)
  F_t = 1 - exp(-(a1*t)^a2)
  
  z = (1-(1-pi^gamma)*F_t)^(1/gamma) |> as.vector()
  
  return(z)
  
}


# Uncured Survival Function -----------------------
su = function(X, beta, gamma, t, alpha){
  
  sp = sp(X, beta, gamma, t, alpha)
  pi = pi_fun(X, beta)
  
  z = (sp - pi)/(1-pi)
  
  z = ifelse(z < 10^-9, 10^-9, z)
  
  return(z)
}

# Uncured pdf ----------------------------------
fu = function(X, beta, gamma, t, alpha){
  
  a1 = alpha[1]
  a2 = alpha[2]
  
  sp = sp(X, beta, gamma, t, alpha)
  pi = pi_fun(X, beta)
  f = a1*a2 * (t*a1)^(a2-1) * exp(-(t*a1)^a2)
  
  z = ((1-pi^gamma) / (1-pi)) * (sp^(1-gamma) * f) / gamma
  
  return(z) # For gamma=1 it is just the Weibull pdf!
}


# Uncured hazard function --------------------------------------
hu = function(X, beta, gamma, t, alpha){ 
  
  fu = fu(X, beta, gamma, t, alpha)
  su = su(X, beta, gamma, t, alpha)
  
  z = fu/su
  
  return(z) # For gamma=1 it is just the Weibull hazard!
  
}

loglike_surv = function(alpha, beta, gamma, X, y, t, delta){ # Loglikelihood
  
  pi = pi_fun(X, beta)
  
  su = su(X, beta, gamma, t, alpha)
  hu = hu(X, beta, gamma, t, alpha)
  
  
  hu = ifelse(hu < 10^(-9), 10^(-9), hu)
  
  z = (1-y) * log(1-pi) + y * log(pi) + delta * log(hu) + y * log(su)
  
  result = sum(z)
  
  return(result)
  
}


loglike = function(alpha, beta, gamma, X, y, t, delta){ # Loglikelihood
  
  
  z1 = loglike_surv(alpha = alpha, beta = beta, gamma = gamma, X = X, y = y, t = t, delta = delta)
  z2 = loglike_glm(parameter = beta, X = X, y = y)
  
  z = z1 + z2
  
  result = sum(z)
  
  return(result)
  
}


# -------------------------------------------------------------
# Optimization using optim function
# ------------------------------------------------------------- 

# Since optim minimizes by default, multiply the function by −1 for maximization

loglike_neg = function(parameter, X, y, t, delta){
  
  alpha = parameter[1:2]
  beta = parameter[3:5]
  gamma = parameter[6]
  
  
  z = -loglike(alpha, beta, gamma, X, y, t, delta)
  
  return(z)
  
}


tmp = lapply(df, function(mydf){
  
  X = as.matrix(cbind(rep(1, nrow(mydf)), mydf[, c('z1', 'z2')])) # Design matrix
  
  
  res = optim(
    par = c(1, 1, 1, 1, 1, 0.5), # Initialization
    fn = loglike_neg, # Function for minimization
    X = X,
    y = mydf$Cure_status,
    t = mydf$Observed_time,
    delta = mydf$Censor_status)
  
  z = cbind(c('a1', 'a2', 'intercept', 'z1', 'z2', 'gamma'), res$par) %>% as.data.frame()
  
  colnames(z) = c('Coef', 'Estimate')
  
  return(z)
  
})

tmp = data.table::rbindlist(tmp, idcol = 'Dataset')

tmp$Estimate = as.numeric(tmp$Estimate)

tmp %>%
  
  group_by(Coef) %>%
  
  rstatix::get_summary_stats(Estimate)


