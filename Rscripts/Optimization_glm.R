# This R script only optimizes the loglikelihood over the incidence part (glm)

# It actually compares the "optim" function with the "glm"

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

# -------------------------------------------------------------
# Optimization using optim function
# ------------------------------------------------------------- 

# Since optim minimizes by default, multiply the function by −1 for maximization

loglike_neg = function(parameter, X, y){
  
  z = -loglike_glm(parameter, X, y)
  
  return(z)
  
}


tmp = lapply(df, function(mydf){
  
  X = as.matrix(cbind(rep(1, nrow(mydf)), mydf[, c('z1', 'z2')])) # Design matrix
  
  
  res = optim(
    par = c(1,1,1), # Initialization
    fn = loglike_neg, # Function for minimization
    X = X,
    y = mydf$Cure_status)
  
  z = cbind(c('Intercept', 'z1', 'z2'), res$par) %>% as.data.frame()
  
  colnames(z) = c('Coef', 'Estimate')
  
  return(z)
  
})

tmp = data.table::rbindlist(tmp, idcol = 'Dataset')

# -------------------------------------------------------------
# Optimization using glm function
# ------------------------------------------------------------- 


tmp2 = lapply(df, function(mydf){
  
  X = as.matrix(cbind(rep(1, nrow(mydf)), mydf[, c('z1', 'z2')])) # Design matrix
  
  mod = glm(mydf$Cure_status ~ X, data = mydf, family = 'binomial')
  
  z = coef(mod)
  
  z = z[which(!is.na(z))]
  
  z = cbind('Coef' = names(z), 'Estimate2' = z) |> as.data.frame()
  
  z$Coef = gsub('X', '', z$Coef)
  
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







