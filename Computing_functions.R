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
  eS <- eigen(inv_Sigma, symmetric = TRUE) 
  ev <- (eS$values)^(-1)
  #if (!all(ev >= -tol * abs(ev[1]))) 
  #  stop("inv_Sigma is not positive definite")
  z <- t(x - mu)
  logdetS <- try(- determinant(inv_Sigma, logarithm = TRUE)$modulus,
                 silent=TRUE)
  attributes(logdetS) <- NULL
  #iS <- MASS::ginv(Sigma)
  ssq <- diag(t(z) %*% inv_Sigma %*% z)
  loglik <- -(n * (log(2*pi)) +  logdetS + ssq)/2
  if (log) return(loglik) else return(exp(loglik))
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

mat_dist = function(x,y)
{ ## return the matrice of distances between all pairs of vectors in x and y
  outer(x, y, Vectorize(function(p,q) sum((p - q)^2)) ) %>% return()
}

kern_to_cov = function(x, kern = kernel, theta = c(1, 0.5), sigma = 0.2)
{ ## x : vector or tibble of all timestamps for each individuals (input vector). 2 columns, 'ID' and 'Timestamp' required
  ## kern : indicates which kernel function to use to compute the covariance matrix
  ## theta : list of the required hyperparameters
  ####
  ## return : list of inverse covariance matrices (1 by individual) of the input vectors according to kernel() function
  
  if(is.vector(x))
  {    
    x = x %>% sort()
    M = mat_dist(x,x)
    mat = kern(M, theta) + diag(sigma^2, length(x))
    
    if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
    rownames(mat) = paste0('X', x)
    colnames(mat) = paste0('X', x)
    return(mat)
  }

  ## If x is a tibble composed of observations from different individuals, identified by a variable 'ID'
  ## In this case, the list of HP must provide a set of HP for each individuals
  else
  { 
    loop = function(i)
    {
      if(theta %>% is.vector())
      {
        theta = list(c(theta, sigma))
        names(theta) = i
      }
      
      indiv = x %>% filter(ID == i) %>% arrange(Timestamp) %>% pull(Timestamp) %>% sort()
      M = mat_dist(indiv,indiv)
      mat = kern(M, theta[[i]][1:2]) + diag((theta[[i]][[3]])^2, length(indiv))
      
      if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
      rownames(mat) = paste0('X', indiv)
      colnames(mat) = paste0('X', indiv)
      return(mat)
    }
    list_mat = lapply(unique(x$ID), loop)
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
      x = x %>% sort()
      M = mat_dist(x,x)
      mat = kern(M, theta) + diag(sigma^2, length(x))
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
      loop = function(i)
      {
        if(theta %>% class() == 'numeric')
        {
          theta = list(c(theta, sigma))
          names(theta) = i
        }
        
        indiv = x %>% filter(ID == i) %>% arrange(Timestamp) %>% pull(Timestamp) %>% sort()
        M = mat_dist(indiv,indiv)
        mat = kern(M, theta[[i]][1:2]) + diag((theta[[i]][[3]])^2, length(indiv))
        #inv = tryCatch(solve(mat), error = function(e){MASS::ginv(mat)})
        inv = solve(mat)
        
        if(inv %>% dim() %>% is.null()){inv = as.matrix(inv)}
        rownames(inv) = paste0('X', indiv)
        colnames(inv) = paste0('X', indiv)
        return(inv)
      }
      list_mat = lapply(unique(x$ID), loop)
      if(length(list_mat) == 1){return(list_mat[[1]])}
      
      return(list_mat)
    } 
}

##### LOGLIKELIHOOD FUNCTIONS ####
logL_GP_mod = function(param, db, mean, kern, new_cov)
{ ## param : vector or list of parameters of the kernel with format (a, b, sigma)
  ## db : tibble containing values we want to compute logL on. Required columns : Timestamp, Output
  ## mean : mean of the GP at corresponding timestamps
  ## kern : kernel used to compute the covariance matrix at corresponding timestamps
  ## new_cov : posterior covariance matrix of the mean GP (mu_0). Used to compute correction term (cor_term)
  ####
  ##return : value of the modified Gaussian log-likelihood for one GP as it appears in the model
  if(length(mean) == 1){mean = rep(mean, nrow(db))} ## mean is equal for all timestamps
  if(length(param) == 2){param[3] = 0.01} ## mean GP (mu_0) is noiseless and thus has only 2 hp
  
  cov =  kern_to_cov(db$Timestamp, kern, theta = param[1:2], sigma = param[3]) 
  #inv = tryCatch(solve(cov), error = function(e){MASS::ginv(cov)}) ## fast or slow matrix inversion if singular
  inv = solve(cov)
  LL_norm = - dmvnorm(db$Output, mean, inv, log = T) ## classic gaussian loglikelihood
  cor_term =  0.5 * (inv * new_cov) %>% sum() ## correction term (0.5 * Trace(inv %*% new_cov))
  
  return(LL_norm + cor_term)
}

logL_multi_GP = function(param, db, kern_i, kern_0, m_hat, K_hat, m_0 = 0)
{ ## param : list of parameters of the kernel for each individuals. Format : list(theta_0, list(theta_i)_i)
  ## db : tibble containing values we want to compute logL on. Required columns : Timestamp, Output
  ## kern_i : kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i)
  ## kern_0 : kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0)
  ## m_hat : posterior mean of the mean GP (mu_0). Needed to compute the log-likelihood
  ## K_hat : posterior covariance matrix of the mean GP (mu_0). Needed to compute correction term
  ## m_0 : prior value of the mean parameter of the mean GP (mu_0). Length = 1 or nrow(db)
  ####
  ##return : value of expectation of joint log-likelihood of the model. The function to be maximised in step M
  
  ## The full likelihood is composed of M+1 independent parts, depending on only theta_0, or theta_i respectively
  ## for each i. The following code computes and sums these M+1 (modified) gaussian likelihoods.

  ll_0 = logL_GP_mod(param$theta_0, db = m_hat, mean = m_0, kern_0, K_hat)

  funloop = function(i)
  { 
    t_i = db %>% filter(ID == i) %>% pull(Timestamp)
    logL_GP_mod(param$theta_i[[i]], db %>% filter(ID == i), mean = m_hat %>% filter(Timestamp %in% t_i) %>% pull(Output),
                kern_i, K_hat[paste0('X', t_i), paste0('X', t_i)]) %>% return()
  }
  sum_ll_i = sapply(unique(db$ID), funloop) %>% sum()
  return(-ll_0 - sum_ll_i)
}


#### GRADIENT OF LogL FUNCTION #####
deriv_hp1 = function(t1, t2, theta = list(1, 0.5), sigma = 0.2)
{
  a = theta[[1]]
  b = theta[[2]]
  
  return( exp(a) * exp((-1/(2 * exp(b))) * t((t1 - t2)) %*% (t1 - t2)) )
}

deriv_hp2 = function(t1, t2, theta = list(1, 0.5), sigma = 0.2)
{
  a = theta[[1]]
  b = theta[[2]]
  
  resi = t(t1 - t2) %*% (t1 - t2)
  norm = 0.5 * exp(-b)
  
  return( exp(a) * norm * resi * exp(- norm * resi) )
}

deriv_sigma = function(t1, t2, theta = list(1, 0.5), sigma = 0.2)
{
  if(all(t1 == t2))
  {
    return(2 * sigma)
  }
  else{return(0)}
}

gr_fn = function(x, db = db_obs)
{
  y = db$Output
  t = db$Timestamp
  inv = kern_to_inv(t, kernel, theta = x)
  
  cste_term = inv %*% y %*% t(inv %*% y) - inv
  
  g_1 = (-1/2 * cste_term %*% kern_to_cov(t, deriv1, theta = x[1:2], sigma = x[3])) %>%  diag() %>% sum()
  g_2 = (-1/2 * cste_term %*% kern_to_cov(t, deriv2, theta = x[1:2], sigma = x[3])) %>%  diag() %>% sum()
  g_3 = (-1/2 * cste_term %*% kern_to_cov(t, deriv3, theta = x[1:2], sigma = x[3])) %>%  diag() %>% sum()
  grad = c(g_1, g_2, g_3)
  
  return(grad)
}

##### TEST #######
# m_hat_test = tibble('Output' = seq(50, 40, length.out = length(unique(db_train$Timestamp))),
#                     'Timestamp' = unique(db_train$Timestamp))
# K_hat_test = kern_to_cov(unique(db_train$Timestamp), kernel, theta = c(1,1,0.2))
 # list_ID = unique(db_train$ID)
 # param_test = list('theta_0' = ini_hp$theta_0, 'theta_i' = ini_hp$theta_i%>% list() %>% rep(length(list_ID)) %>% setNames(nm = list_ID))

# logL_multi_GP(param = param_test, db = db_train, kern_i = kernel,  kern_mu = kernel_mu, m_hat = m_hat_test, 
#               K_hat = K_hat_test) 
# result0 <- opm(param_test, logL_multi_GP, method= c("Nelder-Mead", "L-BFGS-B"), control = list(kkt =FALSE))
# result0 = summary(result0, order=value)
# result0
# 
# result1 <- opm(param_test[c('a_3', 'b_3', 'sigma')], logL_GP, method= meth0, gr = gr_fn, control = list(kkt =FALSE))
# result1 = summary(result1, order=value)
# result1
# 
# pred_gp(db_train %>% filter(ID == 3), timestamps = seq(9,21, 0.03), mean = 45,
#         theta = c(result1$a_3[1], result1$b_3[1]), sigma = result1$sigma[1]) %>%
#         plot_gp(data = db_train %>% filter(ID == 3))

##### EM FUNCTIONS #########
e_step = function(db, m_0, kern_0, kern_i, hp)
{ ## db : full database with all individuals. Columns required : ID, Timestamp, Input, Output
  ## kern_i : kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i)
  ## kern_0 : kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0) 
  ## hp : set of hyper-parameters optimised during the M step
  ####
  ## return : mean and covariance parameters of the mean GP (mu_0)
  all_t = unique(db$Timestamp) %>% sort()
  inv_0 = kern_to_inv(all_t, kern_0, hp$theta_0, sigma = 0.01)
  inv_i = kern_to_inv(db, kern_i, hp$theta_i, sigma = 0)
  value_i = base::split(db$Output, list(db$ID))

  new_inv = update_inv(prior_inv = inv_0, list_inv_i = inv_i)
  #new_cov = tryCatch(solve(new_inv), error = function(e){MASS::ginv(new_inv)}) ## fast or slow matrix inversion if singular
  new_cov = solve(new_inv)
  
  weighted_mean = update_mean(prior_mean = m_0, prior_inv = inv_0, list_inv_i = inv_i, list_value_i = value_i)
  new_mean = new_cov %*% weighted_mean %>% as.vector()
  
  list('mean' = tibble('Timestamp' = all_t, 'Output' = new_mean), 'cov' = new_cov) %>% return()
}

m_step = function(db, old_hp, mean, cov, kern_0, kern_i, m_0)
{ ## db : db : full database with all individuals. Columns required : ID, Timestamp, Input, Output
  ## old_hp : the set of hyper-parameters from the previous step of the EM
  ## mean : mean parameter of the mean GP (mu_0), computed during the E step
  ## cov : covariance parameter of the mean GP (mu_0), computed during the E step
  ## kern_i : kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i)
  ## kern_0 : kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0)
  ## m_0 : prior value of the mean parameter of the mean GP (mu_0). Length = 1 or nrow(db)
  ####
  ## return : set of optimised hyper parameters for the different kernels of the model
  
  t1 = Sys.time()
  new_theta_0 = opm(old_hp$theta_0, logL_GP_mod, db = mean, mean = m_0, kern = kern_0, new_cov = cov,
                    method = "Nelder-Mead", control = list(kkt = FALSE, maxit = 25))[1,1:2]
  
  funloop = function(i)
  {
    t_i = db %>% filter(ID == i) %>% pull(Timestamp)
    return(opm(old_hp$theta_i[[i]] %>% unlist(), logL_GP_mod, db = db %>% filter(ID == i),
                           mean = mean %>% filter(Timestamp %in% t_i) %>% pull(Output), kern = kern_i,
                           new_cov = cov[paste0('X', t_i), paste0('X', t_i)], method = "Nelder-Mead", 
                           control = list(kkt = F, maxit = 10))[1,1:3])
  }
  new_theta_i = sapply(unique(db$ID), funloop, simplify = FALSE, USE.NAMES = TRUE)
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
    common_times = intersect(row.names(x), row.names(new_inv))
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

  for(j in length(list_value_i) %>% seq_len()) 
  {
    weighted_i = list_inv_i[[j]] %*% list_value_i[[j]]
    #row.names(weithed_i) = row.names(list_inv_i[[j]])
    
    common_times = intersect(row.names(weighted_i), row.names(weighted_mean))
    
    weighted_mean[common_times,] = weighted_mean[common_times,] + weighted_i[common_times,]
  }
  return(weighted_mean)
}


##### Old functions #####
# kern_to_cov = function(x, kern = kernel, theta = list(1, 0.5), sigma = 0.2)
# { ## x : vector or tibble of all timestamps for each individuals (input vector). 2 columns, 'ID' and 'Timestamp' required
# ## kern : indicates which kernel function to use to compute the covariance matrix
# ## theta : list of the required hyperparameters
# ## return : list of inverse covariance matrices (1 by individual) of the input vectors according to kernel() function
# 
# ## If user uses kern_to_cov() directly with an observation's vector for one individual
# if(is.vector(x))
# { 
#   x = x %>% sort()
#   mat = sapply(x, function(t) sapply(x, function(s) kern(t, s, theta = theta))) + diag(sigma^2, length(x))
#   if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
#   rownames(mat) = paste0('X', x)
#   colnames(mat) = paste0('X', x)
#   return(mat)
# }
# 
# ## If the tibble is composed of observations from only one individual
# if(x %>% pull(ID) %>% unique() %>% length() == 1)
# {
#   x = x %>% arrange(Timestamp)
#   mat = sapply(x$Timestamp, function(t) sapply(x$Timestamp, function(s) kern(t, s, theta = theta))) + 
#     diag(sigma^2, nrow(x))
#   if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
#   rownames(mat) = paste0('X', x$Timestamp)
#   colnames(mat) = paste0('X', x$Timestamp)
#   return(mat)
# }
# 
# ## If the tibble is composed of observations from different individuals, identified by a variable 'ID'
# ## In this case, the list of HP must provide a set of HP for each individuals
# else
# { 
#   loop = function(i)
#   {
#     indiv = x %>% filter(ID == i) %>% arrange(Timestamp) %>% pull(Timestamp)
#     
#     mat = sapply(indiv, function(t) sapply(indiv, function(s) kern(t, s, theta[[i]][1:2]))) +
#       diag((theta[[i]][3])^2, length(indiv))
#     
#     if(mat %>% dim() %>% is.null()){mat = as.matrix(mat)}
#     rownames(mat) = paste0('X', indiv)
#     colnames(mat) = paste0('X', indiv)
#     return(mat)
#   }
#   list_mat = lapply(unique(x$ID), loop)
#   return(list_mat)
# } 
# }
