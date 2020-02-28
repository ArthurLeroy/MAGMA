library(tidyverse)
library(MASS)
library(Matrix)
library(mvtnorm)

setwd(dir = 'C:/Users/user/CloudStation/Maths/These/Processus Gaussiens/Code R/Algo multitask GP')
source('Computing_functions.R')


training_VEM = function(db, prior_mean, ini_hp, kern_0 = kernel_mu, kern_i = kernel)
{ ## db : database with all individuals in training set. Column required : 'ID', Timestamp', 'Input', 'Output'
  ## prior_mean : prior mean parameter of the mean GP (mu_0)
  ## ini_param : initial values of HP for the kernels
  ## kern_0 : kernel associated to covariance functions of the mean GP
  ## kern_i : kernel associated to common covariance functions of all individuals GPs
  ####
  ## return : list of trained HP, boolean to indicate convergence
  n_loop_max = 25
  list_ID = unique(db$ID)
  hp = list('theta_0' = ini_hp$theta_0, 
            'theta_i' = ini_hp$theta_i %>% list() %>% rep(length(list_ID))  %>% setNames(nm = list_ID))
  cv = 'FALSE'
  
  for(i in 1:n_loop_max)
  { 
    print(i)
    param = e_step(db, prior_mean, kern_0, kern_i, hp)   
    new_hp = m_step(db, hp, mean = param$mean, cov = param$cov, kern_0, kern_i, prior_mean)
    
    logL_new = logL_multi_GP(new_hp, db, kern_i, kern_0, param$mean, param$cov, m_0 = prior_mean)
    eps = (logL_new - logL_multi_GP(hp, db, kern_i, kern_0, param$mean, param$cov, m_0 = prior_mean)) / 
      abs(logL_new)
    
    c(logL_new, logL_multi_GP(hp, db, kern_i, kern_0, param$mean, param$cov, m_0 = prior_mean), eps) %>% print()
    #print(eps)
    #print(param$cov %>% det())
    if(eps <= 0){stop('Likelihood descreased')}
    if(eps > 0 & eps < 1e-3)
    {
      cv = 'TRUE'
      break
    }
    hp = new_hp
  }
  return(list('theta_0' = new_hp$theta_0, 'theta_i' = new_hp$theta_i, 'convergence' = cv, 'param' = param))
}

