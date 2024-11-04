library(dplyr)
library(survival)

myfile_wd = gsub('Rscripts', '/Simulated_Datasets', getwd()) # The folder with the output of Simulations.R

myfolders = list.files(myfile_wd)

mydata = list.files(paste0(myfile_wd, '/', myfolders[1]))


df = list()

for (i in mydata){
  
  df[[i]] = readxl::read_xlsx(paste0(myfile_wd, '/', myfolders[1], '/', i))
  
}

# ------------------------------------------------------------------------------
# Apply logistic regression to see the results of the simulated datasets
# ------------------------------------------------------------------------------

tmp = lapply(df, function(x){
  
  model = glm(Cure_status ~ z1 + z2, data = x, family = 'binomial')
  
  z = cbind(names(model$coefficients), model$coefficients) |> as.data.frame()
  
  colnames(z) = c('Coef', 'Estimate')
  
  return(z)
  
})


tmp = data.table::rbindlist(tmp, idcol = 'Dataset')

results = tmp %>%
  
  mutate(Estimate = as.numeric(Estimate)) %>%
  
  group_by(Coef) %>%
  
  summarise(m = mean(Estimate), se = sd(Estimate)/sqrt(length(Estimate)),
            lower = quantile(Estimate, probs = 0.025) |> round(2),
            upper = quantile(Estimate, probs = 0.975) |> round(2)) %>%
  
  ungroup() %>%
  
  mutate(low = (m - 1.96*se)  |> round(2),
         up = (1+1.96*se)  |> round(2)) %>%
  
  transmute(Coef,
            M=m,
            'Bootstrap 95%CI' = paste0('(', lower, ', ', upper, ')'),
            'Wald 95%CI' = paste0('(', low, ', ', up, ')'))



# ------------------------------------------------------------------------------
#    True Proportion of Cured and Censored Cases in the simulated datasets
# ------------------------------------------------------------------------------

tmp2 = lapply(df, function(x){
  
  z = c(((1 - sum(x$Cure_status)/ nrow(x)) *100 )|> round(2),
        ((1 - sum(x$Censor_status)/ nrow(x)) *100 )|> round(2))
  
  z = matrix(z, nrow = 1) |> as.data.frame()
  
  colnames(z) = c('Cure_rate', 'Censoring_proportion')
  
  return(z)
  
}) %>%
  
  data.table::rbindlist(., idcol = 'Dataset')

tmp2 %>%
  
  select(-Dataset) %>%
  
  rstatix::get_summary_stats()


# ------------------------------------------------------------------------------
#    Time to event for the Susceptible in the simulated datasets
# ------------------------------------------------------------------------------

tmp_weibull = lapply(df, function(x){
  
  x = x %>% filter(Cure_status == 1)
  
  failures = x$Observed_time[which(x$Censor_status == 1)]
  censored = x$Observed_time[which(x$Censor_status == 0)]
  
  fit = WeibullR::MLEw2p(x = failures, s = censored)
  
  z = cbind(names(fit), fit) |> as.data.frame()
  
  colnames(z) = c('Coef', 'Estimate')
  
  return(z)
  
})


tmp_weibull = data.table::rbindlist(tmp_weibull, idcol = 'Dataset')

results = tmp_weibull %>%
  
  filter(Coef != 'LL') %>%
  
  mutate(Estimate = as.numeric(Estimate),
         Estimate = ifelse(Coef == 'Eta', 1/Estimate, Estimate)) %>%
  
  group_by(Coef) %>%
  
  summarise(m = mean(Estimate), se = sd(Estimate)/sqrt(length(Estimate)),
            lower = quantile(Estimate, probs = 0.025) |> round(2),
            upper = quantile(Estimate, probs = 0.975) |> round(2)) %>%
  
  ungroup() %>%
  
  mutate(low = (m - 1.96*se)  |> round(2),
         up = (1+1.96*se)  |> round(2)) %>%
  
  transmute(Coef,
            M=m,
            'Bootstrap 95%CI' = paste0('(', lower, ', ', upper, ')'),
            'Wald 95%CI' = paste0('(', low, ', ', up, ')'))


