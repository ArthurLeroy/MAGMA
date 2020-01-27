library(tidyverse)
library(optimr)
library(mvtnorm)

##### USEFUL FUNCTIONS ########
'%notin%' <- Negate('%in%')

##### KERNELS DEFINITION ######
kernel = function(t1, t2, theta = list(1, 0.5, 0.5))
{ ## t1 : first observational timestamp
  ## t2 : second observatinal timestamp
  ## theta : list of hyperparamaters of the kernel
  ## return : the value of dot production <f(t1),f(t2)> computed from the kernel
  a = theta[[1]] 
  b = theta[[2]]
  c = 0
  
  if(all(t1 == t2)){c = exp(theta[[3]])}  #t is the value of pointwise variance (i.e. sigma)
  kern = exp(a) * exp((-1/(2 * exp(b))) * t((t1 - t2)) * (t1 - t2)) + c
  return(kern)
}

#kernel2 = function(t1 , t2, theta)
#{ ## t1 : first observational timestamp
  ## t2 : second observatinal timestamp
  ## theta : list of hyperparamaters of the kernel
  ## return : the value of dot production <f(t1),f(t2)> computed from the kernel
#  a = theta[[1]] 
#  b = theta[[2]]
#  kern = a^2 * exp((-1/(2 * b^2)) * t((t1 - t2)) %*% (t1 - t2))
#  return(kern)
#  }

kern_to_cov = function(x, kern = kernel, theta = list(1, 0.2, 0.4))
{ ## x : vector or tibble of all timestamps for each individuals (input vector). 2 columns, 'ID' and 'Timestamp' required
  ## kern : indicates which kernel function to use to compute the covariance matrix
  ## theta : list of the required hyperparameters
  ####
  ## return : list of inverse covariance matrices (1 by individual) of the input vectors according to kernel() function
  
  ## If user uses kern_to_cov() directly with an observation's vector for one individual
  if(is.vector(x))
  { mat = sapply(x, function(t) sapply(x, function(s) kern(t, s, theta = theta)))
    if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
    rownames(mat) = paste0('X', x)
    colnames(mat) = paste0('X', x)
    return(mat)
  }
  
  ## If the tibble is composed of observations from only one individual
  if(x %>% pull(ID) %>% unique() %>% length() == 1)
  {
    mat = sapply(x$Timestamp, function(t) sapply(x$Timestamp, function(s) kern(t, s, theta = theta)))
    if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
    rownames(mat) = paste0('X', x$Timestamp)
    colnames(mat) = paste0('X', x$Timestamp)
    return(mat)
  }
  
  ## If the tibble is composed of observations from different individuals, identified by a variable 'ID'
  else
  { list_mat = list()
    for(i in unique(x$ID))
    {
      indiv = x %>% filter(ID == i) %>% pull(Timestamp)
      mat = sapply(indiv, function(t) sapply(indiv, function(s) kern(t, s, theta = theta)))
      if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
      rownames(mat) = paste0('X', indiv)
      colnames(mat) = paste0('X', indiv)
      list_mat[[i]] = mat
    }
  return(list_mat)
  } 
}

kern_to_inv = function(x, kern = kernel, theta = list(1, 0.2, 0.4))
{ ## x : vector or tibble of all timestamps for each individuals (input vector). 2 columns, 'ID' and 'Age' required
  ## kern : indicates which kernel function to use to compute the covariance matrix
  ## theta : list of the required hyperparameters
  ####
  ## return : list of inverse covariance matrices (1 by individual) of the input vectors according to kernel()
  
  ## If user uses kern_to_inv() directly with an observation's vector for one individual
  if(is.vector(x))
  { 
    mat = sapply(x, function(t) sapply(x, function(s) kern(t, s, theta = theta))) %>% solve()
    if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
    rownames(mat) = paste0('X', x)
    colnames(mat) = paste0('X', x)
    return(mat)
  }
  
  ## If the tibble is composed of observations from only one individual
  if(x %>% pull(ID) %>% unique() %>% length() == 1)
  {
    mat = sapply(x$Timestamp, function(t) sapply(x$Timestamp, function(s) kern(t, s, theta = theta))) %>% solve()
    if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
    rownames(mat) = paste0('X', x$Timestamp)
    colnames(mat) = paste0('X', x$Timestamp)
    return(mat)
  }
  
  ## If the tibble is composed of observations from different individuals, identified by a variable 'ID'
  else
  { list_mat = list()
    for(i in unique(x$ID))
    {
      indiv = x %>% filter(ID == i) %>% pull(Timestamp)
      mat = sapply(indiv, function(t) sapply(indiv, function(s) kern(t, s, theta = theta))) %>% solve()
      if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
      rownames(mat) = paste0('X', indiv)
      colnames(mat) = paste0('X', indiv)
      list_mat[[i]] = mat
    }
  return(list_mat)
  } 
}

##### LOGLIKELIHOOD FUNCTIONS ####
logL = function(y , t , theta_i, sigma , kern = kernel)
{
  cov = kern_to_cov(t, kern = kern, theta = theta_i)
  return(1/2 * t(y) %*% solve(cov) %*% y + 1/2 %*% det(cov))
}

##### EM FUNCTIONS ############
e_step_GP = function(t, kern, theta)
{
   return(kern_to_inv(t, kern = kern, theta))
}

m_step_GP = function(mean, inv)
{
  
}

#### EM MEAN GP FUNCTIONS ####
update_inv = function(prior_inv, list_inv_i)
{
  new_inv = prior_inv
  
  for(x in list_inv_i)
  {
    inv_i = x
    common_times = intersect(row.names(x), row.names(new_inv))
    new_inv[common_times, common_times] = new_inv[common_times, common_times] + inv_i[common_times, common_times]
  }
  return(new_inv)
}

update_mean = function(prior_mean, prior_inv, list_inv_i, list_value_i)
{
  if(length(prior_mean) == 1){prior_mean = rep(prior_mean, ncol(prior_inv))}
  weighted_mean = prior_inv %*% prior_mean
  #row.names(weithed_mean) = row.names(prior_inv)

  for(j in 1:length(list_value_i)) 
  {
    weighted_i = list_inv_i[[j]] %*% list_value_i[[j]]
    #row.names(weithed_i) = row.names(list_inv_i[[j]])
    
    common_times = intersect(row.names(weighted_i), row.names(weighted_mean))
    
    weighted_mean[common_times,] = weighted_mean[common_times,] + weighted_i[common_times,]
  }
  return(weighted_mean)
}

m_step_meanGP = function(mean, cov)
{
  return(new_theta_0)
}
