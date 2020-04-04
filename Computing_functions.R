library(tidyverse)
library(optimr)
library(Matrix)
library(MASS)

##### USEFUL FUNCTIONS ########
'%notin%' <- Negate('%in%')

dmvnorm <- function (x, mu, inv_Sigma, log = FALSE, tol = 1e-06) 
{
  if (is.vector(x)) 
    x = t(as.matrix(x))
  n = length(mu)
  if (is.vector(mu)) {
    p <- length(mu)
    if (is.matrix(x)) {
      mu <- matrix(rep(mu, nrow(x)), ncol = p, byrow = TRUE)
    }
  }
  else {
    p <- ncol(mu)
  }
  if (!all(dim(inv_Sigma) == c(p, p)) || nrow(x) != nrow(mu)) 
    stop("incompatible arguments")
  
  z <- t(x - mu)
  logdetS <- try(- determinant(inv_Sigma, logarithm = TRUE)$modulus,
                 silent=TRUE)
  attributes(logdetS) <- NULL
  
  ssq <- t(z) %*% inv_Sigma %*% z
  loglik <- -(n * (log(2*pi)) +  logdetS + ssq)/2
  if (log) return(loglik) else return(exp(loglik))
}

mat_dist = function(x,y)
{ ## return the matrice of distances between all pairs of vectors in x and y
  outer(x, y, Vectorize(function(p,q) sum((p - q)^2)) ) %>% return()
}

##### KERNELS DEFINITION ######
kernel = function(mat, theta = c(1, 0.5))
{ ## mat : the matrix M[i,j] =  t(i - j) %*% (i - j), for all i,j in the vector of timestamps
  ## theta : list of hyperparamaters of the kernel
  ####
  ## return : the value of dot production <f(t1),f(t2)> computed from the kernel
  
  exp(theta[[1]] - exp(-theta[[2]])/2 * mat) %>% return()
}

kernel_mu = function(mat, theta = c(1,0.5))
{ ## mat : the matrix M[i,j] =  t(i - j) %*% (i - j), for all i,j in the vector of timestamps
  ## theta : list of hyperparamaters of the kernel
  ####
  ## return : the value of dot production <f(t1),f(t2)> computed from the kernel
  
  exp(theta[[1]] - exp(-theta[[2]])/2 * mat) %>% return()
}

kern_to_cov = function(x, kern = kernel, theta = c(1, 0.5), sigma = 0.2)
{ ## x : vector or tibble of all timestamps for each individuals (input vector). 2 columns, 'ID' and 'Timestamp' required
  ## kern : indicates which kernel function to use to compute the covariance matrix
  ## theta : list of the required hyperparameters
  ####
  ## return : list of inverse covariance matrices (1 by individual) of the input vectors according to kernel() function
  
  if(is.vector(x))
  { 
    if(length(x) == 1)
    {
      mat = sigma^2 + exp(theta[[1]])
    } 
    else 
    {
      x = x %>% sort()
      M = dist(x)^2
      mat = as.matrix(kern(M, theta)) + diag(sigma^2 + exp(theta[[1]]) ,length(x))
    }
    
    
    if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
    rownames(mat) = paste0('X', x)
    colnames(mat) = paste0('X', x)
    return(mat)
  }
  
  ## If x is a tibble composed of observations from different individuals, identified by a variable 'ID'
  ## In this case, the list of HP must provide a set of HP for each individuals
  else
  { 
    floop = function(i)
    {
      if(theta %>% is.vector())
      {
        theta = list(c(theta, sigma))
        names(theta) = i
      }
      
      indiv = x %>% filter(ID == i) %>% arrange(Timestamp) %>% pull(Timestamp) %>% sort()
      if(length(indiv) == 1)
      {
        mat = theta[[i]][3]^2 + exp(theta[[i]][1])
      } 
      else 
      {
        M = dist(indiv)^2
        mat = as.matrix(kern(M, theta[[i]][1:2])) + diag(theta[[i]][3]^2 + exp(theta[[i]][1]), length(indiv))
      }

      if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
      rownames(mat) = paste0('X', indiv)
      colnames(mat) = paste0('X', indiv)
      return(mat)
    }
    list_mat = sapply(unique(x$ID), floop, simplify = FALSE, USE.NAMES = TRUE)
    if(length(list_mat) == 1){return(list_mat[[1]])}
    
    return(list_mat)
  } 
}

kern_to_inv = function(x, kern = kernel, theta = c(1, 0.5), sigma = 0.2)
{ ## x : vector or tibble of all timestamps for each individuals (input vector). 2 columns, 'ID' and 'Timestamp' required
  ## kern : indicates which kernel function to use to compute the covariance matrix
  ## theta : list of the required hyperparameters
  ####
  ## return : list of inverse covariance matrices (1 by individual) of the input vectors according to kernel() function
  if(is.vector(x))
  { 
    if(length(x) == 1)
    {
      mat = sigma^2 + exp(theta[[1]])
    } 
    else 
    {
      x = x %>% sort()
      M = dist(x)^2
      mat = as.matrix(kern(M, theta)) + diag(sigma^2 + exp(theta[[1]]) ,length(x))
    }
    #inv = tryCatch(solve(mat), error = function(e){MASS::ginv(mat)})
    inv = solve(mat)
    
    if(inv %>% dim() %>% is.null()){inv = as.matrix(inv)}
    rownames(inv) = paste0('X', x)
    colnames(inv) = paste0('X', x)
    return(inv)
  }
  
  ## If x is a tibble composed of observations from different individuals, identified by a variable 'ID'
  ## In this case, the list of HP must provide a set of HP for each individuals
  else
  { 
    floop = function(i)
    { 
      if(theta %>% class() == 'numeric')
      {
        theta = list(c(theta, sigma))
        names(theta) = i
      }
      
      indiv = x %>% filter(ID == i) %>% arrange(Timestamp) %>% pull(Timestamp) %>% sort()
      if(length(indiv) == 1)
      {
        mat = theta[[i]][3]^2 + exp(theta[[i]][1])
      } 
      else 
      {
        M = dist(indiv)^2
        mat = as.matrix(kern(M, theta[[i]][1:2])) + diag(theta[[i]][3]^2 + exp(theta[[i]][1]), length(indiv))
      }
      #inv = tryCatch(solve(mat), error = function(e){MASS::ginv(mat)})
      inv = solve(mat)
      
      if(inv %>% dim() %>% is.null()){inv = as.matrix(inv)}
      rownames(inv) = paste0('X', indiv)
      colnames(inv) = paste0('X', indiv)
      return(inv)
    }
    list_mat = sapply(unique(x$ID), floop, simplify = FALSE, USE.NAMES = TRUE)
    
    if(length(list_mat) == 1){return(list_mat[[1]])}
    
    return(list_mat)
  } 
}


##### LOGLIKELIHOOD FUNCTIONS ####
logL_GP_mod = function(hp, db, mean, kern, new_cov)
{ ## hp : vector or list of parameters of the kernel with format (a, b, sigma)
  ## db : tibble containing values we want to compute logL on. Required columns : Timestamp, Output
  ## mean : mean of the GP at corresponding timestamps
  ## kern : kernel used to compute the covariance matrix at corresponding timestamps
  ## new_cov : posterior covariance matrix of the mean GP (mu_0). Used to compute correction term (cor_term)
  ####
  ##return : value of the modified Gaussian log-likelihood for one GP as it appears in the model

  if(length(mean) == 1){mean = rep(mean, nrow(db))} ## mean is equal for all timestamps
  sigma = ifelse((length(hp) == 3), hp[[3]], 0.1) ## mean GP (mu_0) is noiseless and thus has only 2 hp

  inv =  kern_to_inv(db$Timestamp, kern, theta = hp[1:2], sigma) 
  
  LL_norm = - dmvnorm(db$Output, mean, inv, log = T) ## classic gaussian loglikelihood
  cor_term =  0.5 * (inv * new_cov) %>% sum() ## correction term (0.5 * Trace(inv %*% new_cov))
  
  return(LL_norm + cor_term)
}

logL_GP_mod_common_hp = function(hp, db, mean, kern, new_cov)
{ ## hp : vector of common hyperparameters for all individuals. Format : c(a, b, sigma)
  ## db : tibble of data. Required columns : ID, Timestamp, Output
  ## mean : mean of the GP at union of observed timestamps
  ## kern : kernel used to compute the covariance matrix at corresponding timestamps
  ## new_cov : posterior covariance matrix of the mean GP (mu_0). Used to compute correction term (cor_term)
  ####
  ##return : value of the modified Gaussian log-likelihood for the sum of all indiv with same HPs
  
  sigma = ifelse((length(hp) == 3), hp[[3]], 0.1) ## mean GP is noiseless (0.1 is for computational issues) has only 2 hp

  LL_norm = 0
  cor_term = 0
  t_i_old = NULL

  for(i in unique(db$ID))
  {
    t_i = db %>% filter(ID == i) %>% pull(Timestamp)
    input_i = paste0('X', t_i)
    y_i = db %>% filter(ID == i) %>% pull(Output)
    
    if( !identical(t_i, t_i_old) )
    { ## We update the inverse cov matrix only if necessary (if different timestamps)
      inv =  kern_to_inv(t_i, kern, theta = hp[1:2], sigma)
    }

    LL_norm = LL_norm - dmvnorm(y_i, mean %>% filter(Timestamp %in% t_i) %>% pull(Output), inv, log = T) 
    cor_term = cor_term + 0.5 * (inv * new_cov[input_i, input_i]) %>% sum()  ##(0.5 * Trace(inv %*% new_cov))
    
    t_i_old = t_i
  }
  return(LL_norm + cor_term)
}

logL_monitoring = function(hp, db, kern_i, kern_0, mean_mu, cov_mu, m_0)
{ ## hp : list of parameters of the kernel for each individuals. Format : list(theta_0, list(theta_i)_i)
  ## db : tibble containing values we want to compute logL on. Required columns : Timestamp, Output
  ## kern_i : kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i)
  ## kern_0 : kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0)
  ## mean_mu : posterior mean of the mean GP (mu_0). Needed to compute the log-likelihood
  ## cov_mu : posterior covariance matrix of the mean GP (mu_0). Needed to compute correction term
  ## m_0 : prior value of the mean parameter of the mean GP (mu_0). Length = 1 or nrow(db)
  ####
  ##return : value of expectation of joint log-likelihood of the model. The function to be maximised in step M
  
  ## The full likelihood is composed of M+1 independent parts, depending on only theta_0, or theta_i respectively
  ## for each i. The following code computes and sums these M+1 (modified) gaussian likelihoods.

  ll_0 = logL_GP_mod(hp$theta_0, db = mean_mu, mean = m_0, kern_0, cov_mu)

  funloop = function(i)
  { 
    t_i = db %>% filter(ID == i) %>% pull(Timestamp)
    logL_GP_mod(hp$theta_i[[i]], db %>% filter(ID == i),
                mean = mean_mu %>% filter(Timestamp %in% t_i) %>% pull(Output),
                kern_i, cov_mu[paste0('X', t_i), paste0('X', t_i)]) %>% return()
  }
  sum_ll_i = sapply(unique(db$ID), funloop) %>% sum()
  
  return(-ll_0 - sum_ll_i)
}

#### GRADIENT OF LogL FUNCTION #####
deriv_hp1 = function(mat, theta)
{
  exp(theta[[1]] - exp(-theta[[2]])/2 * mat) %>% return()
}

deriv_hp2 = function(mat, theta)
{
  grad = 0.5 * exp(- theta[[2]]) * mat
  
  return( exp(theta[[1]]) * grad * exp(- grad) )
}

gr_GP_mod = function(hp, db, mean, kern, new_cov)
{ 
  y = db$Output
  t = db$Timestamp
  
  sigma = ifelse((length(hp) == 3), hp[[3]], 0.1)
  inv = kern_to_inv(t, kern, theta = hp[1:2], sigma)
  prod_inv = inv %*% (y - mean) 
  cste_term = prod_inv %*% t(prod_inv) + inv %*% ( new_cov %*% inv - diag(1, length(t)) )
  
  
  g_1 = 1/2 * (cste_term %*% kern_to_cov(t, deriv_hp1, theta = hp[1:2], sigma = 0)) %>%  diag() %>% sum()
  g_2 = 1/2 * (cste_term %*% as.matrix(deriv_hp2(dist(t)^2, theta = hp[1:2]) ))  %>%  diag() %>% sum()
  
  if(length(hp) == 3)
  {
    g_3 = hp[[3]] * (cste_term %>% diag() %>% sum() )
    (- c(g_1, g_2, g_3)) %>% return()
  }
  else   (- c(g_1, g_2)) %>% return()
}

gr_GP_mod_common_hp = function(hp, db, mean, kern, new_cov)
{ 
  sigma = ifelse((length(hp) == 3), hp[[3]], 0.1)
  g_1 = 0
  g_2 = 0
  g_3 = 0
  t_i_old = NULL
  
  for(i in unique(db$ID))
  {
    t_i = db %>% filter(ID == i) %>% pull(Timestamp)
    input_i = paste0('X', t_i)
    y_i = db %>% filter(ID == i) %>% pull(Output)
    
    if( !identical(t_i, t_i_old) )
    { ## We update the inverse cov matrix only if necessary (if different timestamps)
      inv = kern_to_inv(t_i, kern, theta = hp[1:2], sigma)
    }
    prod_inv = inv %*% (y_i - mean %>% filter(Timestamp %in% t_i) %>% pull(Output)) 
    cste_term = prod_inv %*% t(prod_inv) + inv %*% 
                ( new_cov[input_i,input_i] %*% inv - diag(1, length(t_i)) )
    
    g_1 = g_1 + 1/2 * (cste_term %*% kern_to_cov(t_i, deriv_hp1, theta = hp[1:2], sigma = 0)) %>%  diag() %>% sum()
    g_2 = g_2 + 1/2 * (cste_term %*% as.matrix(deriv_hp2(dist(t_i)^2, theta = hp[1:2]) )) %>%  diag() %>% sum()
    
    t_i_old = t_i
    if(length(hp) == 3)
    {
      g_3 = g_3 + hp[[3]] * (cste_term %>% diag() %>% sum() )
    }
  }
  if(length(hp) == 3) return(- c(g_1, g_2, g_3)) else  return(- c(g_1, g_2))
}

##### EM FUNCTIONS #########
e_step = function(db, m_0, kern_0, kern_i, hp)
{ ## db : full database with all individuals. Columns required : ID, Timestamp, Input, Output
  ## kern_i : kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i)
  ## kern_0 : kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0) 
  ## hp : set of hyper-parameters optimised during the M step
  ####
  ## return : mean and covariance parameters of the mean GP (mu_0)
  all_t = unique(db$Timestamp) %>% sort()
  inv_0 = kern_to_inv(all_t, kern_0, hp$theta_0, sigma = 0.1)
  inv_i = kern_to_inv(db, kern_i, hp$theta_i, sigma = 0)
  value_i = base::split(db$Output, list(db$ID))

  new_inv = update_inv(prior_inv = inv_0, list_inv_i = inv_i)
  #new_cov = tryCatch(solve(new_inv), error = function(e){MASS::ginv(new_inv)}) ## fast or slow matrix inversion if singular
  new_cov = solve(new_inv)
  
  weighted_mean = update_mean(prior_mean = m_0, prior_inv = inv_0, list_inv_i = inv_i, list_value_i = value_i)
  new_mean = new_cov %*% weighted_mean %>% as.vector()
  
  list('mean' = tibble('Timestamp' = all_t, 'Output' = new_mean), 'cov' = new_cov) %>% return()
}

m_step = function(db, old_hp, mean, cov, kern_0, kern_i, m_0, common_hp)
{ ## db : db : full database with all individuals. Columns required : ID, Timestamp, Input, Output
  ## old_hp : the set of hyper-parameters from the previous step of the EM
  ## mean : mean parameter of the mean GP (mu_0), computed during the E step
  ## cov : covariance parameter of the mean GP (mu_0), computed during the E step
  ## kern_i : kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i)
  ## kern_0 : kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0)
  ## m_0 : prior value of the mean parameter of the mean GP (mu_0). Length = 1 or nrow(db)
  ####
  ## return : set of optimised hyper parameters for the different kernels of the model
  list_ID = unique(db$ID)

  t1 = Sys.time()
  new_theta_0 = opm(old_hp$theta_0, logL_GP_mod, gr = gr_GP_mod, db = mean, mean = m_0, kern = kern_0, new_cov = cov,
                    method = "L-BFGS-B", control = list(kkt = FALSE))[1,1:2]

  if(common_hp)
  {
    param = opm(old_hp$theta_i[[1]], logL_GP_mod_common_hp, gr = gr_GP_mod_common_hp , db = db, mean = mean,
                kern = kern_i, new_cov = cov, method = "L-BFGS-B", control = list(kkt = F))[1,1:3]
    new_theta_i = param %>% list() %>% rep(length(list_ID))  %>% setNames(nm = list_ID) 
  }
  else
  {
    floop = function(i)
    {
      t_i = db %>% filter(ID == i) %>% pull(Timestamp)
      return(opm(old_hp$theta_i[[i]] %>% unlist(), logL_GP_mod, gr = gr_GP_mod , db = db %>% filter(ID == i),
                 mean = mean %>% filter(Timestamp %in% t_i) %>% pull(Output), kern = kern_i,
                 new_cov = cov[paste0('X', t_i), paste0('X', t_i)], method = "L-BFGS-B", 
                 control = list(kkt = F))[1,1:3])
    }
    new_theta_i = sapply(list_ID, floop, simplify = FALSE, USE.NAMES = TRUE)
  }
                
  t2 = Sys.time()
  print(t2-t1)
  list('theta_0' = new_theta_0, 'theta_i' = new_theta_i) %>% return()
}

#### UPDATE FUNCTIONS ####
update_inv = function(prior_inv, list_inv_i)
{ ## prior_inv : inverse of the covariance matrix of the prior mean GP (mu_0). dim = all timestamps 
  ## list_inv_i : list of inverse of the covariance matrices of each individuals. dim = timestamps of i  
  ####
  ## return : inverse of the covariance of the posterior mean GP (mu_0 | (y_i)_i). dim = (all timestamps)^2 
  
  new_inv = prior_inv
  
  for(x in list_inv_i)
  {
    inv_i = x
    common_times = intersect(row.names(inv_i), row.names(new_inv))
    new_inv[common_times, common_times] = new_inv[common_times, common_times] + inv_i[common_times, common_times]
  }
  return(new_inv)
}

update_mean = function(prior_mean, prior_inv, list_inv_i, list_value_i)
{ ## prior_mean : mean parameter of the prior mean GP (mu_0)
  ## prior_inv : inverse of the covariance matrix of the prior mean GP (mu_0). dim = (all timestamps)^2 
  ## list_inv_i : list of inverse of the covariance matrices of each individuals. dim = (timestamps of i)^2  
  ## list_value_i : list of outputs (y_i) for each individuals. dim = (timestamps of i) x 1
  ####
  ## return : mean parameter of the posterior mean GP (mu_0 | (y_i)_i). dim = (all timestamps) x 1
  
  if(length(prior_mean) == 1){prior_mean = rep(prior_mean, ncol(prior_inv))}
  weighted_mean = prior_inv %*% prior_mean
  #row.names(weithed_mean) = row.names(prior_inv)

  for(i in list_inv_i %>% names()) 
  {
    weighted_i = list_inv_i[[i]] %*% list_value_i[[i]]
    #row.names(weithed_i) = row.names(list_inv_i[[i]])
    
    common_times = intersect(row.names(weighted_i), row.names(weighted_mean))
    weighted_mean[common_times,] = weighted_mean[common_times,] + weighted_i[common_times,]
  }
  return(weighted_mean)
}
