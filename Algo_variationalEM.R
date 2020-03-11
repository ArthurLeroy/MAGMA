library(tidyverse)
library(MASS)
library(Matrix)
library(mvtnorm)

#setwd(dir = 'C:/Users/user/CloudStation/Maths/These/Processus Gaussiens/Code R/Algo multitask GP')
source('Computing_functions_VEM.R')

##### TRAINING FUNCTIONS ####
training_VEM = function(db, prior_mean, ini_hp, kern_0 = kernel_mu, kern_i = kernel, ini_tau_i_k)
{ ## db : database with all individuals in training set. Column required : 'ID', Timestamp', 'Input', 'Output'
  ## prior_mean : prior mean parameter of the K mean GPs (mu_k)
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
  
  tau_i_k = ini_tau_i_k
  pi_k = ### compute them with the ini_tau_i_k
  
  for(i in 1:n_loop_max)
  { 
    print(i)
    ## E-Step
    param = e_step(db, prior_mean, kern_0, kern_i, hp, tau_i_k, pi_k)  
    
    ## Monitoring of the LL
    logLL_complete = logL_multi_GP(hp, db, kern_i, kern_0, param$mean, param$cov, m_0 = prior_mean) - 
      0.5 * ( nrow(param$cov) + log(det(param$cov)) ) 
    c(logL_new, logLL_complete, eps) %>% print()
    ## M-Step
    new_hp = m_step(db, hp, mean = param$mean, cov = param$cov, kern_0, kern_i, prior_mean, param$tau)
    
    ## Testing the stoping condition
    logL_new = logL_multi_GP(new_hp, db, kern_i, kern_0, param$mean, param$cov, m_k = prior_mean)
    eps = (logL_new - logL_multi_GP(hp, db, kern_i, kern_0, param$mean, param$cov, m_k = prior_mean)) / 
          abs(logL_new)

    if(eps <= 0){stop('Likelihood descreased')}
    if(eps > 0 & eps < 1e-3)
    {
      cv = 'TRUE'
      break
    }
    hp = new_hp
    tau_i_k = param$tau_i_k
    pi_k = new_hp$pi_k
  }
  return(list('theta_0' = new_hp$theta_0, 'theta_i' = new_hp$theta_i, 'convergence' = cv, 'param' = param))
}

train_new_gp_VEM = function(db, mean_mu, cov_mu, ini_hp, kern_i)
{
  LL_GP<- function(hp, db, kern) 
  {
    return(-dmvnorm(db$Output, mean_mu, solve(kern_to_cov(db$Timestamp, kern, theta = hp[1:2], sigma = hp[3]) + cov_mu),
                    log = T))
  }
  new_hp = opm(ini_hp, LL_GP, db = db, kern = kern_i, method = "Nelder-Mead", control = list(kkt = FALSE))[1,1:3] 
  
  tau_k = ### Voir M-step
  
  list('theta' = new_hp[1:2], 'sigma' = new_hp[3], 'tau_k' = tau_k) %>% return()
}

################ PREDICTION FUNCTIONS ################
posterior_mu_k = function(db, timestamps, m_k, kern_0, kern_i, hp, tau_i_k)
{ ## db : matrix of data columns required ('Timestamp', 'Input', 'Output')(Input format : paste0('X', Timestamp))
  ## timestamps : timestamps on which we want a prediction
  ## prior_mean : prior mean value of the mean GP (scalar value or vector of same length as 'timestamps')
  ## kern_0 : kernel associated to the covariance function of the mean GP
  ####
  ## return : pamameters of the mean GP at timestamps
  inv_0 = kern_to_inv(timestamps, kern_0, hp$theta_0, sigma = 0.001)
  inv_i = kern_to_inv(db, kern_i, hp$theta_i, sigma = 0)
  value_i = base::split(db$Output, list(db$ID))
  
  new_inv = update_inv(prior_inv = inv_0, list_inv_i = inv_i)
  new_cov = tryCatch(solve(new_inv), error = function(e){MASS::ginv(new_inv)}) ## fast or slow matrix inversion if singular
  rownames(new_cov) = rownames(new_inv)
  colnames(new_cov) = colnames(new_inv)
  
  weighted_mean = update_mean(prior_mean = m_0, prior_inv = inv_0, list_inv_i = inv_i, list_value_i = value_i)
  new_mean = new_cov %*% weighted_mean
  
  #names(mean_mu) = paste0('X', t_mu)
  list('Timestamp' = timestamps, 'Mean' = new_mean, 'Cov' = new_cov) %>% return()
}

pred_gp_clust = function(db, timestamps, list_mu, kern = kernel, theta = list(1,0.2), sigma = 0.2, tau_k)
{ ## db : tibble of data columns required ('Timestamp', 'Input', 'Output')(Input format : paste0('X', Timestamp))
  ## timestamps : timestamps on which we want a prediction
  ## mean_mu : mean value of mean GP at timestamps (obs + pred) (matrix dim: timestamps x 1, with Input rownames)
  ## cov_mu : covariance of mean GP at timestamps (obs + pred) (square matrix, with Input row/colnames)
  ## kernel : kernel associated to the covariance function of the GP
  ## theta : list of hyperparameters for the kernel of the GP
  ## sigma : variance of the error term of the models
  ## tau_k : vector of probability to belong to each cluster for the new individual to predict
  ####
  ## return : pamameters of the gaussian density predicted at timestamps 
  tn = db %>% pull(Timestamp)
  input = db %>% pull(Input)
  yn = db %>% pull(Output)
  input_t = paste0('X', timestamps)
  all_times = union(tn,timestamps)
  
  if(is.null(cov_mu))
  {
    cov_mu = matrix(0, length(all_times), length(all_times),
                    dimnames = list(paste0('X', all_times), paste0('X', all_times)))
  }
  if(is.null(mean_mu)){mean_mu = matrix(0, length(all_times), 1, dimnames = list(paste0('X', all_times)))}
  if(length(mean_mu) == 1){mean_mu = matrix(mean_mu, length(all_times), 1, dimnames = list(paste0('X', all_times)))}
  
  inv_mat = (kern_to_cov(tn, kern, theta, sigma) + cov_mu[input, input]) %>% solve()
  cov_tn_t = kern(mat_dist(tn, timestamps), theta) + cov_mu[input,input_t]
  cov_t_t = kern(mat_dist(timestamps, timestamps), theta) + cov_mu[input_t ,input_t]
  
  tibble('Timestamp' = timestamps,
         'Mean' = mean_mu[input_t,] + t(cov_tn_t) %*% inv_mat %*% (yn - mean_mu[input,]) %>% as.vector(),
         'Var' = (cov_t_t - t(cov_tn_t) %*% inv_mat %*% cov_tn_t) %>% diag()) %>% return()
}

################ PLOT FUNCTIONS ######################
plot_gp = function(pred_gp, data = NULL)
{ ## pred_gp : tibble coming out of the pred_gp() function, columns required : 'Timestamp', 'Mean', 'Var'
  ## data : tibble of observational data, columns required : 'Timestamp', 'Output'
  ####
  ## return : plot of the predicted curve of the GP with the 0.95 confidence interval (optional : data points)
  gg = ggplot() +
    geom_line(data = pred_gp, aes(x = Timestamp, y = Mean), color = 'blue') +
    geom_ribbon(data = pred_gp, aes(x = Timestamp, ymin = Mean - 1.96* sqrt(Var), 
                                    ymax = Mean +  1.96* sqrt(Var)), alpha = 0.2)
  
  if(!is.null(data)){gg = gg + geom_point(data = data, aes(x = Timestamp, y = Output, col = ID), shape = 4)}
  return(gg)
}

################ APPLICATION #########################

full_algo_clust = function(db, new_db, timestamps, kern_i, pi_k, plot = T, 
                     prior_mean = NULL, kern_0 = NULL, list_hp = NULL, mu = NULL, ini_hp = NULL, hp_new_i = NULL)
{ ## db : Database containing all training data from all individuals. Column: ID - Timestamp - Output.
  ## new_db : Database containing data for a new individual we want a prediction on.
  ## timestamps : Timestamps we want to predict at.
  ## kern_i : Kernel associated to individual GPs.
  ## pi_k : Initial probability to belong to each cluster for individuals.
  ## plot : Boolean indicating whether we want to display a graph at the end of computations.
  ## prior_mean : Prior arbitrary value for the mean processes. Optional, not needed if 'mu' is given.
  ## kern_0 : Kernel associated to the mean GPs. Optional, not needed if 'mu' is given.
  ## list_hp : Hyper-parameters for all individuals in training set. Optional, computed if NULL.
  ## mu : Database containing parameters of the mean GPs at all prediction timestamps. Optional, computed if NULL.
  ## ini_hp : Initial values of the HP to start the training. Optional, not needed if 'list_hp' is given.
  ## hp_new_i : Hyper-pameters for the new individual to predict. Optional, computed if NULL. 
  ####
  ## return : predicted GP parameters | posterior K mean processes | all trained hyperparameters
  if(is.null(list_hp)){list_hp = training(db, prior_mean, ini_hp, kern_0, kern_i)[c('theta_k', 'theta_i')]}
  
  t_pred = timestamps %>% union(unique(db$Timestamp)) %>% union(unique(new_db$Timestamp)) %>% sort()
  if(is.null(mu)){mu = posterior_mu(db, t_pred, prior_mean, kern_0, kern_i, list_hp)}
  mean_mu = mu$Mean
  t_mu = mu$Timestamp
  obs_input = new_db$Input
  
  if(is.null(hp_new_i)){hp_new_i = train_new_gp(new_db, mean_mu = mean_mu[obs_input,], 
                                                cov_mu = mu$Cov[obs_input,obs_input], ini_hp$theta_i, kern_i)}
  
  pred = pred_gp(new_db, timestamps, mean_mu, mu$Cov, kern_i, hp_new_i$theta, hp_new_i$sigma)
  
  if(plot)
  {
    (pred %>% plot_gp(data = rbind(new_db, db)) + geom_point(mu, aes(Timestamp, Output, col = ID))) %>% print()
  }
  
  list('Prediction' = pred, 'Mean_process' = mu_k, 'Hyperparameters' = list('other_i' = list_hp, 'new_i' = hp_new_i)) %>% 
    return()
}
