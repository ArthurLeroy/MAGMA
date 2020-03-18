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
    floop = function(i)
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
    floop = function(i)
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
  if(length(hp) == 2){hp[3] = 0.01} ## mean GP (mu_0) is noiseless and thus has only 2 hp
  
  cov =  kern_to_cov(db$Timestamp, kern, theta = hp[1:2], sigma = hp[3]) 
  #inv = tryCatch(solve(cov), error = function(e){MASS::ginv(cov)}) ## fast or slow matrix inversion if singular
  inv = solve(cov)
  LL_norm = - dmvnorm(db$Output, mean, inv, log = T) ## classic gaussian loglikelihood
  cor_term =  0.5 * (inv * new_cov) %>% sum() ## correction term (0.5 * Trace(inv %*% new_cov))
  
  return(LL_norm + cor_term)
}

logL_multi_GP = function(hp, db, kern_i, kern_0, mu_k_param, m_k)
{ ## hp : list of parameters of the kernel for each individuals. Format : list(theta_0, list(theta_i)_i)
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
  
  floop = function(k)
  {
    logL_GP_mod(hp$theta_k[[k]], db = mu_k_param$mean[[k]], mean = m_k[[k]] , kern_0, mu_k_param$cov[[k]]) %>%
    return()
  }
  sum_ll_k = sapply(names(m_k), floop) %>% sum()
  
  floop2 = function(i)
  { 
    t_i = db %>% filter(ID == i) %>% pull(Timestamp)
    logL_clust_multi_GP(hp$theta_i[[i]], db %>% filter(ID == i), mu_k_param, kern_0) %>% return()
  }
  sum_ll_i = sapply(unique(db$ID), floop2) %>% sum()
  return(-sum_ll_k - sum_ll_i)
}

logL_clust_multi_GP = function(hp, db, mu_k_param, kern)
{ ## hp : vector or list of parameters of the kernel with format (a, b, sigma)
  ## db : tibble containing values we want to compute logL on. Required columns : Timestamp, Output
  ## mu_k_param : list of parameters for the K mean Gaussian processes
  ## kern : kernel used to compute the covariance matrix at corresponding timestamps
  ####
  ##return : value of the modified Gaussian log-likelihood for one GP as it appears in the model
  
  if(length(hp) == 2){hp[3] = 0.01} ## mean GP (mu_0) is noiseless and thus has only 2 hp
  names_k = mu_k_param$mean %>% names()
  t_i = db$Timestamp
  y_i = db$Output
  i = unique(db$ID)
  
  inv =  kern_to_inv(db$Timestamp, kern, theta = hp[1:2], sigma = hp[3]) 
  #inv = tryCatch(solve(cov), error = function(e){MASS::ginv(cov)}) ## fast or slow matrix inversion if singular
  
  LL_norm = - dmvnorm(y_i, rep(0, length(y_i)), inv, log = T) ## classic gaussian centered loglikelihood

  corr1 = 0
  corr2 = 0
  
  for(k in seq_len(length(names_k)))
  {
    tau_i_k = mu_k_param$tau_i_k[[k]][[i]]
    mean_mu_k = mu_k_param$mean[[k]] %>% filter(Timestamp %in% t_i) %>% pull(Output)
    corr1 = corr1 + tau_i_k * mean_mu_k
    corr2 = corr2 + tau_i_k * ( mean_mu_k %*% t(mean_mu_k) + mu_k_param$cov[[k]][paste0('X', t_i), paste0('X', t_i)] )
  }
  return( LL_norm - y_i %*% inv %*% corr1 + 0.5 * sum(diag(inv %*% corr2)) )
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
e_step_VEM = function(db, m_k, kern_0, kern_i, hp, old_tau_i_k)
{ ## db : full database with all individuals. Columns required : ID, Timestamp, Input, Output
  ## m_k : prior means of the mu_k processes. List of vector of size K
  ## kern_i : kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i)
  ## kern_0 : kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0) 
  ## hp : set of hyper-parameters theta_k, theta_i, and pi_k (the proportion of individuals in cluster k)
  ## old_tau_i_k : values au tau_i_k from previous iterations. List(list(tau)_i)_k
  ## 
  ####
  ## return : mean and covariance parameters of the mean GP (mu_0)
  pi_k = hp$pi_k
  all_t = unique(db$Timestamp) %>% sort()
  t_clust = tibble('ID' = rep(names(m_k), each = length(all_t)) , 'Timestamp' = rep(all_t, length(m_k)),
                   'Input' = rep(paste0('X',all_t), length(m_k)))
  inv_k = kern_to_inv(t_clust, kern_0, hp$theta_k, sigma = 0)
  inv_i = kern_to_inv(db, kern_i, hp$theta_i, sigma = 0)
  value_i = base::split(db$Output, list(db$ID))
  
  ## Update each mu_k parameters
  floop = function(k)
  {
    new_inv = update_inv_VEM(prior_inv = inv_k[[k]], list_inv_i = inv_i, old_tau_i_k[[k]])
    #new_cov = tryCatch(solve(new_inv), error = function(e){MASS::ginv(new_inv)}) ## fast or slow matrix inversion if singular
    solve(new_inv) %>% return()
  }
  cov_k = sapply(names(m_k), floop, simplify = FALSE, USE.NAMES = TRUE)
  
  floop2 = function(k)
  {
    weighted_mean = update_mean_VEM(m_k[[k]], inv_k[[k]], inv_i, value_i, old_tau_i_k[[k]])
    new_mean = cov_k[[k]] %*% weighted_mean %>% as.vector()
    tibble('Timestamp' = all_t, 'Output' = new_mean) %>% return()
  }
  mean_k = sapply(names(m_k), floop2, simplify = FALSE, USE.NAMES = TRUE)
  
  ## Update tau_i_k
  tau_i_k = update_tau_i_k_VEM(db, m_k, mean_k, cov_k, kern_i, hp, pi_k)
  
  list('mean' = mean_k, 'cov' = cov_k, 'tau_i_k' = tau_i_k) %>% return()
}

m_step_VEM = function(db, old_hp, list_mu_param, kern_0, kern_i, m_k)
{ ## db : db : full database with all individuals. Columns required : ID, Timestamp, Input, Output
  ## old_hp : the set of hyper-parameters from the previous step of the EM
  ## list_mu_param : List of parameters of the K mean GPs. Format list('mean', 'cov', 'tau_i_k')
  ## kern_i : kernel used to compute the covariance matrix of individuals GP at corresponding timestamps (Psi_i)
  ## kern_0 : kernel used to compute the covariance matrix of the mean GP at corresponding timestamps (K_0)
  ## m_k : prior value of the mean parameter of the mean GPs (mu_k). Length = 1 or nrow(db)
  ####
  ## return : set of optimised hyper parameters for the different kernels of the model
  
  t1 = Sys.time()
  funloop = function(k)
  {
    print(c('k',k, Sys.time() - t1))
    new_theta_k = c(opm(old_hp$theta_k[[k]][1:2], logL_GP_mod, db = list_mu_param$mean[[k]], mean = m_k[[k]], 
                      kern = kern_0, new_cov = list_mu_param$cov[[k]] , method = "Nelder-Mead",
                      control = list(kkt = FALSE))[1,1:2] %>% unlist(), 0.1)
  }
  new_theta_k = sapply(names(m_k), funloop, simplify = FALSE, USE.NAMES = TRUE)
  
  funloop2 = function(i)
  {
    print(c('i',i, Sys.time() - t1))
    t_i = db %>% filter(ID == i) %>% pull(Timestamp)
    return(opm(old_hp$theta_i[[i]] %>% unlist(), logL_clust_multi_GP, db = db %>% filter(ID == i),
               mu_k_param = list_mu_param, kern = kern_i, method = "Nelder-Mead", 
               control = list(kkt = F, maxit = 30))[1,1:3])
  }
  new_theta_i = sapply(unique(db$ID), funloop2, simplify = FALSE, USE.NAMES = TRUE)
  
  pi_k = sapply( list_mu_param$tau_i_k, function(x) x %>% unlist() %>% mean() ) 

  t2 = Sys.time()
  print(t2 - t1)
  list('theta_k' = new_theta_k, 'theta_i' = new_theta_i, 'pi_k' = pi_k) %>% return()
}

##### UPDATE FUNCTIONS ####
update_inv_VEM = function(prior_inv, list_inv_i, tau_i_k)
{ ## prior_inv : inverse of the covariance matrix of the prior mean GP (mu_0). dim = all timestamps 
  ## list_inv_i : list of inverse of the covariance matrices of each individuals. dim = timestamps of i 
  ## tau_i_k : list of probability for each indiv 'i' to belong to cluster k
  ####
  ## return : inverse of the covariance of the posterior mean GP (mu_0 | (y_i)_i). dim = (all timestamps)^2 
  
  new_inv = prior_inv
  
  for(x in list_inv_i %>% names())
  {
    inv_i = list_inv_i[[x]]
    common_times = intersect(row.names(inv_i), row.names(new_inv))
    new_inv[common_times, common_times] = new_inv[common_times, common_times] + 
                                          tau_i_k[[x]] * inv_i[common_times, common_times]
  }
  return(new_inv)
}

update_mean_VEM = function(prior_mean, prior_inv, list_inv_i, list_value_i, tau_i_k)
{ ## prior_mean : mean parameter of the prior mean GP (mu_0)
  ## prior_inv : inverse of the covariance matrix of the prior mean GP (mu_0). dim = (all timestamps)^2 
  ## list_inv_i : list of inverse of the covariance matrices of each individuals. dim = (timestamps of i)^2  
  ## list_value_i : list of outputs (y_i) for each individuals. dim = (timestamps of i) x 1
  ## tau_i_k : list of probability for each indiv 'i' to belong to cluster k 
  ####
  ## return : mean parameter of the posterior mean GP (mu_0 | (y_i)_i). dim = (all timestamps) x 1
  
  if(length(prior_mean) == 1){prior_mean = rep(prior_mean, ncol(prior_inv))}
  weighted_mean = prior_inv %*% prior_mean
  #row.names(weithed_mean) = row.names(prior_inv)
  
  for(x in list_inv_i %>% names()) 
  {
    weighted_i = tau_i_k[[x]] * list_inv_i[[x]] %*% list_value_i[[x]]
    #row.names(weithed_i) = row.names(list_inv_i[[j]])
    
    common_times = intersect(row.names(weighted_i), row.names(weighted_mean))
    weighted_mean[common_times,] = weighted_mean[common_times,] + weighted_i[common_times,]
  }
  return(weighted_mean)
}

update_tau_i_k_VEM = function(db, m_k, mean_k, cov_k, kern_i, hp, pi_k)
{
  c_i = 0
  c_k = 0
  mat_logL = matrix(NA, nrow = length(names(m_k)), ncol = length(unique(db$ID)) )
  
  for(i in unique(db$ID))
  { c_i = c_i + 1
  
    for(k in names(m_k))
    { c_k = c_k + 1
      t_i = db %>% filter(ID == i) %>% pull(Timestamp)
      mat_logL[c_k,c_i] = - logL_GP_mod(hp$theta_i[[i]], db %>% filter(ID == i),
                                       mean_k[[k]] %>% filter(Timestamp %in% t_i) %>% pull(Output), kern_i,
                                       cov_k[[k]][paste0('X',t_i),paste0('X',t_i)]) %>% exp() 
    } 
    c_k = 0
  } 
 
  (pi_k * mat_logL) %>% apply(2,function(x) x / sum(x)) %>%
                        `rownames<-`(names(m_k)) %>%
                        `colnames<-`(unique(db$ID)) %>% 
                        apply(1, as.list) %>%
                        return()
}


### Tests ####
### Test update functions
# list_inv_i_test = kern_to_inv(db_train)
# prior_inv_test = kern_to_inv(seq(10,20,0.5), theta = c(2, 0.5), sigma = 0.1)
# tau_k_test = list('1' = 0.1, '2' = 0.6 , '3' = 0.3, '4' = 0.9 , '5' = 0.4, '6' = 0.8, '7' = 0.5,
#                     '8' = 0.3, '9' = 0.1, '10' = 0.2)
# list_value_i_test = base::split(db_train$Output, list(db_train$ID))
# 
# bla = update_inv_VEM(prior_inv_test, list_inv_i_test, tau_k_test)
# fu = update_mean_VEM(100, prior_inv_test, list_inv_i_test, list_value_i_test, tau_i_k_test)
# solve(bla) %*% fu
### Test EM functions
# ini_hp_test = list('theta_k' = c(2,0.5,0), 'theta_i' = c(1, 1, 0.2))
# k = seq_len(5)
# m_k_test = as.vector(rep(45, length(k)), 'list')
# names(m_k_test) = paste0('K', k)
# pi_k_test = runif(length(k))
# pi_k_test = pi_k_test / sum(pi_k_test)
# hp_test = list('theta_k' = ini_hp_test$theta_k %>% list() %>% rep(length(names(m_k_test)))  %>% setNames(nm = names(m_k_test)),
#                'theta_i' = ini_hp_test$theta_i %>% list() %>% rep(length(unique(db_train$ID)))  %>% setNames(nm = unique(db_train$ID)))
# 
# tau_i_k_test = replicate(length(k), runif(length(unique(db_train$ID)))) %>%
#                apply(2,function(x) x / sum(x)) %>%
#               `colnames<-`(paste0('K', k)) %>%
#               `rownames<-`(unique(db_train$ID)) %>%
#               apply(2, as.list)

# e_step_VEM(db_train, m_k_test, kernel_mu, kernel, hp_test, tau_i_k_test, pi_k_test)
#
# list_mu_param_test = e_step_VEM(db_train, m_k_test, kernel_mu, kernel, hp_test, tau_i_k_test, pi_k_test)
# m_step_VEM(db_train, hp_test, list_mu_param_test, kernel_mu, kernel, m_k_test, tau_i_k_test)
