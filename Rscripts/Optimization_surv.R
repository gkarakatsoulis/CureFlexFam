# This R script only optimizes the loglikelihood over the prevalence part (Weibull)

# It actually compares the "optim" function with the "WeibullR::MLEw2p"

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

# Function for the cure probabilities
pi_fun = function(X, beta){
  
  z = 1/( 1 + exp(-X %*% beta) ) |> as.vector()
  
  return(z)
  
}

# Population Survival Function ---------------------------------
sp = function(X, beta, gamma = 1, t, alpha){
  
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
fu = function(X, beta, gamma = 1, t, alpha){
  
  a1 = alpha[1]
  a2 = alpha[2]
  
  sp = sp(X, beta, gamma, t, alpha)
  pi = pi_fun(X, beta)
  f = a1*a2 * (t*a1)^(a2-1) * exp(-(t*a1)^a2)
  
  z = ((1-pi^gamma) / (1-pi)) * (sp^(1-gamma) * f) / gamma
  
  return(z) # For gamma=1 it is just the Weibull pdf!
}


# Uncured hazard function --------------------------------------
hu = function(X, beta, gamma = 1, t, alpha){ 
  
  fu = fu(X, beta, gamma, t, alpha)
  su = su(X, beta, gamma, t, alpha)
  
  z = fu/su
  
  return(z) # For gamma=1 it is just the Weibull hazard!
  
}

loglike_surv = function(alpha, beta, gamma = 1, X, y, t, delta){ # Loglikelihood
  
  pi = pi_fun(X, beta)
  
  su = su(X, beta, gamma, t, alpha)
  hu = hu(X, beta, gamma, t, alpha)
  
  
  hu = ifelse(hu < 10^(-9), 10^(-9), hu)
  
  z = y * log(1-pi) + (1-y) * log(pi) + delta * log(hu) + y * log(su)
  
  result = sum(z)
  
  return(result)
  
}



# -------------------------------------------------------------
# Optimization using optim function
# ------------------------------------------------------------- 

# Since optim minimizes by default, multiply the function by −1 for maximization

loglike_neg = function(alpha, beta, gamma = 1, X, y, t, delta){
  
  z = -loglike_surv(alpha, beta, gamma, X, y, t, delta)
  
  return(z)
  
}


tmp = lapply(df, function(mydf){
  
  X = as.matrix(cbind(rep(1, nrow(mydf)), mydf[, c('z1', 'z2')])) # Design matrix
  
  
  res = optim(
    par = c(1,1), # Initialization
    fn = loglike_neg, # Function for minimization
    beta = c(1,1,1),
    X = X,
    y = mydf$Cure_status,
    t = mydf$Observed_time,
    delta = mydf$Censor_status)
  
  z = cbind(c('a1', 'a2'), res$par) %>% as.data.frame()
  
  colnames(z) = c('Coef', 'Estimate')
  
  return(z)
  
})

tmp = data.table::rbindlist(tmp, idcol = 'Dataset')

tmp$Estimate = as.numeric(tmp$Estimate)

# -------------------------------------------------------------
# Optimization using WeibullR::MLEln2p function
# ------------------------------------------------------------- 


tmp2 = lapply(df, function(mydf){
  
  
  mod = WeibullR::MLEw2p(x = mydf$Event_time[which(mydf$Censor_status == 1)],
                          s = mydf$Observed_time[which(mydf$Censor_status == 0 & mydf$Cure_status == 1)])
  
  z = mod[c('Eta', 'Beta')]
  
  z = c(1/z[1], z[2])
  
  names(z) = c('a1','a2')
  
  z = cbind('Coef' = names(z), 'Estimate2' = z) |> as.data.frame()
  
  
  return(z)
  
})

tmp2 = data.table::rbindlist(tmp2, idcol = 'Dataset')




# -------------------------------------------------------------
# Results comparison
# ------------------------------------------------------------- 


results = merge(tmp, tmp2)

results$Estimate = as.numeric(results$Estimate)
results$Estimate2 = as.numeric(results$Estimate2)
results$diff = results$Estimate - results$Estimate2

summary(results$diff)







