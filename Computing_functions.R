library(tidyverse)

##### USEFUL FUNCTIONS ########
'%notin%' <- Negate('%in%')
###############################


##### KERNELS DEFINITION ######
kernel = function(t1, t2, theta = list(1, 0.2))
{ ## t1 : first observational timestamp
  ## t2 : second observatinal timestamp
  ## theta : list of hyperparamaters of the kernel
  ## return : the value of dot production <f(t1),f(t2)> computed from the kernel
  a = theta[[1]] 
  b = theta[[2]]
  
  t = 0
  ##if(t1 == t2){t = 0.1}  t is the value of pointwise variance (i.e. sigma)
  kern = a * exp((-1/(2 * b)) * t((t1 - t2)) * (t1 - t2))  + t
  return(kern)
}

#kernel2 = function(t1 , t2, theta)
{ ## t1 : first observational timestamp
  ## t2 : second observatinal timestamp
  ## theta : list of hyperparamaters of the kernel
  ## return : the value of dot production <f(t1),f(t2)> computed from the kernel
  a = theta[[1]] 
  b = theta[[2]]
  
  t = 0
  if(t1 == t2){t = 0.1}
  kern = a * exp((-1/(2 * b)) * t((t1 - t2)) * (t1 - t2)) + t
  return(kern)
}

kern_to_cov = function(x, kern = kernel, theta = list(1, 0.2))
{ ## x : vector or tibble of all timestamps for each indivicuals (input vector). 2 columns, 'ID' and 'Timestamp' required
  ## kern : indicates which kernel function to use to compute the covariance matrix
  ## theta : list of the required hyperparameters
  ####
  ## return : list of inverse covariance matrices (1 by individual) of the input vectors according to kernel() function
  
  ## If user use ker_to_cov() directly with a observations vector for one individual
  if(is.vector(x))
  { mat = sapply(x, function(t) sapply(x, function(s) kernel(t, s, theta = theta)))
    if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
    rownames(mat) = paste0('X', x)
    colnames(mat) = paste0('X', x)
    return(mat)
  }
  
  ## If the tibble is composed of observations from only one individual
  if(x %>% pull(ID) %>% unique() %>% length() == 1)
  {
    mat = sapply(x$Timestamp, function(t) sapply(x$Timestamp, function(s) kernel(t, s, theta = theta)))
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
      mat = sapply(indiv, function(t) sapply(indiv, function(s) kernel(t, s, theta = theta)))
      if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
      rownames(mat) = paste0('X', indiv)
      colnames(mat) = paste0('X', indiv)
      list_mat[[i]] = mat
    }
  return(list_mat)
  } 
}

kern_to_inv = function(x, kern = kernel, theta = list(1, 0.2))
{ ## x : vctor or tibble of all timestamps for each indivicuals (input vector). 2 columns, 'ID' and 'Age' required
  ## kern : indicates which kernel function to use to compute the covariance matrix
  ## theta : list of the required hyperparameters
  ####
  ## return : list of inverse covariance matrices (1 by individual) of the input vectors according to kernel()
  
  ## If user use ker_to_inv() directly with an observations vector for one individual
  if(is.vector(x))
  { 
    mat = sapply(x, function(t) sapply(x, function(s) kernel(t, s, theta = theta))) %>% solve()
    if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
    rownames(mat) = paste0('X', x)
    colnames(mat) = paste0('X', x)
    return(mat)
  }
  
  ## If the tibble is composed of observations from only one individual
  if(x %>% pull(ID) %>% unique() %>% length() == 1)
  {
    mat = sapply(x$Timestamp, function(t) sapply(x$Timestamp, function(s) kernel(t, s, theta = theta))) %>% solve()
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
      mat = sapply(indiv, function(t) sapply(indiv, function(s) kernel(t, s, theta = theta))) %>% solve()
      if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
      rownames(mat) = paste0('X', indiv)
      colnames(mat) = paste0('X', indiv)
      list_mat[[i]] = mat
    }
  return(list_mat)
  } 
}
###############################

##### LOGLIKELIHOOD FUNCTIONS ####
logL_i = function(y = , t = , theta_i = , sigma = , kern = kernel)
{
  cov = kern_to_cov(t, kern = kern, theta = theta_i)
  return(-1/2 * t(y) %*% solve(cov) %*% y - 1/2 %*% det(cov) )
}

logL_0 = function(obs = , theta_0)
{
  
}
###############################

##### EM FUNCTIONS ###########
e_step_GP = function(t, kern, theta)
{
   return(kern_to_inv(t, kern = kern, theta = theta))
}

m_step_GP = function(mean, inv)
{
  
}
##############################


#### EM MEAN GP FUNCTIONS ####
update_inv = function(old_inv, list_inv_i)
{
  
  return(new_inv)
}

update_mean = function(old_mean, old_inv, list_value_i)
{
  weighted_m = old_inv %*% old_mean
    
  for(i in list_value_i)
  {
   # sum de tous les indices communs
  }
  return(weighted_m)
}

m_step_meanGP = function(mean, cov)
{
  return(new_theta_0)
}
##############################
