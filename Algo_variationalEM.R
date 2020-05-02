library(tidyverse)
library(MASS)
library(Matrix)
library(mvtnorm)
library(optimr)

#setwd(dir = 'C:/Users/user/Google Drive/Travail/GitHub/gpclust/')
source('Computing_functions_VEM.R')

##### TRAINING FUNCTIONS ####
training_VEM = function(db, prior_mean_k, ini_hp, kern_0, kern_i, ini_tau_i_k, common_hp_k = T, common_hp_i = T)
{ ## db : database with all individuals in training set. Column required : 'ID', Timestamp', 'Input', 'Output'
  ## prior_mean : prior mean parameter of the K mean GPs (mu_k)
  ## ini_hp : initial values of HP for the kernels
  ## kern_0 : kernel associated to covariance functions of the mean GP
  ## kern_i : kernel associated to common covariance functions of all individuals GPs
  ## ini_tau_i_k : initial values of probabiliy to belong to each cluster for each individuals. 
  ####
  ## return : list of trained HP, boolean to indicate convergence
  n_loop_max = 25
  list_ID = unique(db_train$ID)
  
  ID_k = names(prior_mean_k)
  hp = list('theta_k' = ini_hp$theta_k %>% list() %>% rep(length(ID_k))  %>% setNames(nm = ID_k), 
            'theta_i' = ini_hp$theta_i %>% list() %>% rep(length(list_ID))  %>% setNames(nm = list_ID))
  cv = 'FALSE'
  tau_i_k = ini_tau_i_k
  hp[['pi_k']] = sapply( tau_i_k, function(x) x %>% unlist() %>% mean() ) 
  logLL_monitoring = - Inf

  for(i in 1:n_loop_max)
  { 
    print(i)
    ## E-Step
    param = e_step_VEM(db, prior_mean_k, kern_0, kern_i, hp, tau_i_k)  

    ## Monitoring of the LL
    (new_logLL_monitoring = logL_monitoring(hp, db, kern_i, kern_0, mu_k_param = param , m_k = prior_mean_k) + 
                            0.5 * (length(param$cov) * nrow(param$cov[[1]]) +
                                  Reduce('+', lapply(param$cov, function(x) log(det(x)))) )) %>% print()
    diff_moni = new_logLL_monitoring - logLL_monitoring
    
    if(diff_moni < 0){stop('Likelihood descreased')}
    
    ## M-Step
    new_hp = m_step_VEM(db, hp, list_mu_param = param, kern_0, kern_i, prior_mean_k, common_hp_k, common_hp_i)
    
    ## Testing the stoping condition
    logL_new = logL_monitoring(new_hp, db, kern_i, kern_0, mu_k_param = param, m_k = prior_mean_k)
    eps = (logL_new - logL_monitoring(hp, db, kern_i, kern_0, mu_k_param = param, m_k = prior_mean_k)) / 
          abs(logL_new)
    
    print(c('eps', eps))
    if(eps>0 & eps < 1e-3)
    {
      cv = 'TRUE'
      break
    }
    hp = new_hp
    tau_i_k = param$tau_i_k
    logLL_monitoring = new_logLL_monitoring
  }
  return(list('theta_k' = new_hp$theta_k, 'theta_i' = new_hp$theta_i, 'convergence' = cv,  'param' = param))
}

### Sûrement besoin d'un VEM pour le nouvel individu
train_new_gp_VEM = function(db, param_mu_k, ini_hp, kern_i, hp)
{
  mean_mu_k = param_mu_k$mean
  cov_mu_k = param_mu_k$cov
  pi_k = sapply(param_mu_k$tau_i_k, function(x) Reduce("+", x)/ length(x)) 
  
  tau_k = update_tau_i_k_VEM = function(db, m_k, mean_mu_k, cov_mu_k, kern_i, hp, pi_k)
  
  LL_GP<- function(hp, db, kern) 
  {

    
    floop = function(k)
    {
      ## To avoid pathological behaviour of the opm optimization function in rare cases
      if((hp %>% abs %>% sum) > 20){return(10^10)}
      
      - tau_k[[k]] * dmvnorm(db$Output, mean_mu[[k]], 
                             solve(kern_to_cov(db$Timestamp, kern, theta = hp[1:2], sigma = hp[3]) + cov_mu_k[[k]]),
                             log = T) %>% return() 
    }
    sapply(names(mean_mu_k), floop) %>% sum() %>% return()
  }
  new_hp = opm(ini_hp, LL_GP, db = db, kern = kern_i, method = "L-BFGS-B", control = list(kkt = FALSE))[1,1:3] 
  
  
  list('theta' = new_hp[1:2], 'sigma' = new_hp[3], 'tau_k' = tau_k) %>% return()
}

################ PREDICTION FUNCTIONS ################
posterior_mu_k = function(db, timestamps, m_k, kern_0, kern_i, hp)
{ ## db : tibble of data columns required ('ID', 'Timestamp',  'Output')
  ## timestamps : timestamps on which we want a prediction
  ## prior_mean : prior mean value of the mean GP (scalar value or vector of same length as 'timestamps')
  ## kern_0 : kernel associated to the covariance function of the mean GP
  ####
  ## return : pamameters of the mean GP at timestamps
  
  t_clust = tibble('ID' = rep(names(m_k), each = length(timestamps)) , 'Timestamp' = rep(timestamps, length(m_k)),
                   'Input' = rep(paste0('X', timestamps), length(m_k)))
  inv_k = kern_to_inv(t_clust, kern_0, hp$theta_k, sigma = 0)
  inv_i = kern_to_inv(db, kern_i, hp$theta_i, sigma = 0)
  value_i = base::split(db$Output, list(db$ID))
  tau_i_k = hp$param$tau_i_k
  
  ## Update each mu_k parameters
  floop = function(k)
  {
    new_inv = update_inv_VEM(prior_inv = inv_k[[k]], list_inv_i = inv_i, tau_i_k[[k]])
    #new_cov = tryCatch(solve(new_inv), error = function(e){MASS::ginv(new_inv)}) ## fast or slow matrix inversion if singular
    solve(new_inv) %>% return()
  }
  cov_k = sapply(names(m_k), floop, simplify = FALSE, USE.NAMES = TRUE)
  
  floop2 = function(k)
  {
    weighted_mean = update_mean_VEM(m_k[[k]], inv_k[[k]], inv_i, value_i, tau_i_k[[k]])
    new_mean = cov_k[[k]] %*% weighted_mean %>% as.vector()
    tibble('Timestamp' = timestamps, 'Output' = new_mean) %>% return()
  }
  mean_k = sapply(names(m_k), floop2, simplify = FALSE, USE.NAMES = TRUE)
  
  #names(mean_mu) = paste0('X', t_mu)
  list('mean' = mean_k, 'cov' = cov_k) %>% return()
}

pred_gp_clust = function(db, timestamps = NULL, list_mu, kern, hp)
{ ## db : tibble of data columns required ('Timestamp', 'Input', 'Output')(Input format : paste0('X', Timestamp))
  ## timestamps : timestamps on which we want a prediction
  ## mean_mu : mean value of mean GP at timestamps (obs + pred) (matrix dim: timestamps x 1, with Input rownames)
  ## cov_mu : covariance of mean GP at timestamps (obs + pred) (square matrix, with Input row/colnames)
  ## kern : kernel associated to the covariance function of the GP
  ## hp : list of hyperparameters and tau_k for the new individual
  ####
  ## return : pamameters of the gaussian density predicted at timestamps 
  tn = db %>% pull(Timestamp)
  #input = db %>% pull(Input)
  input = paste0('X', db$Timestamp)
  yn = db %>% pull(Output)
  
  if(is.null(timestamps)){timestamps = seq(min(tn), max(tn), length.out = 500)}
  input_t = paste0('X', timestamps)
  
  theta = hp$theta
  sigma = hp$sigma
  tau_k = hp$tau_k
  
  mean = 0
  cov = 0
  for(k in list_mu$mean %>% names())
  {
    mean_mu_obs = list_mu$mean[[k]] %>% filter(Timestamp %in% tn) %>% pull(Output)
    mean_mu_pred = list_mu$mean[[k]] %>% filter(Timestamp %in% timestamps) %>% pull(Output)
    cov_mu = list_mu$cov[[k]]
    inv_mat = (kern_to_cov(tn, kern, theta, sigma) + cov_mu[input, input]) %>% solve() 
    cov_tn_t = kern(mat_dist(tn, timestamps), theta) + cov_mu[input,input_t]
    cov_t_t = kern_to_cov(timestamps, kern, theta, sigma) + cov_mu[input_t ,input_t] 
    
    mean = mean + tau_k[[k]] * (mean_mu_pred + t(cov_tn_t) %*% inv_mat %*% (yn - mean_mu_obs) %>% as.vector())
    cov = cov + tau_k[[k]] * (cov_t_t - t(cov_tn_t) %*% inv_mat %*% cov_tn_t) %>% diag() 
  }
  tibble('Timestamp' = timestamps, 'Mean' = mean, 'Var' = cov) %>% return()
}

pred_gp_clust_animate = function(db, timestamps = NULL, list_mu, kern, hp)
{
  ## Inputs : same as for a classic GP prediction
  ####
  ## return : tibble of classic GP predictions but with an inscreasing number of data points considered as 'observed'
  db %>% arrange(Timestamp)
  all_pred = tibble()
  
  if(is.null(timestamps)){timestamps = seq(min(db$Timestamp), max(db$Timestamp), length.out = 500)}
  
  for(j in 1:nrow(db))
  {
    pred_j = pred_gp_clust(db[1:j,], timestamps, list_mu, kern, hp) %>% mutate(Nb_data = j)
    all_pred = all_pred %>% rbind(pred_j)
  }
  return(all_pred)
}

################ PLOT FUNCTIONS ######################
plot_gp_clust = function(pred_gp, data = NULL, data_train = NULL, mean = NULL, mean_CI = F)
{ ## pred_gp : tibble coming out of the pred_gp() function, columns required : 'Timestamp', 'Mean', 'Var'
  ## data : tibble of observational data, columns required : 'Timestamp', 'Output'
  ####
  ## return : plot of the predicted curve of the GP with the 0.95 confidence interval (optional : data points)
  gg = ggplot() +
    geom_line(data = pred_gp, aes(x = Timestamp, y = Mean), color = 'blue') +
    geom_ribbon(data = pred_gp, aes(x = Timestamp, ymin = Mean - 1.96* sqrt(Var), 
                                    ymax = Mean +  1.96* sqrt(Var)), alpha = 0.2)
  
  if(!is.null(data_train)){gg = gg + geom_point(data = data_train, aes(x = Timestamp, y = Output, col = ID), shape = 4)}
  if(!is.null(data)){gg = gg + geom_point(data = data, aes(x = Timestamp, y = Output), size = 2, shape = 18)}
  if(!is.null(mean))
  {
    for(i in 1:length(mean))
    {
      gg = gg + geom_line(data = mean[[i]], aes(x = Timestamp, y = Output), linetype = 'dashed') 
    }
  }
  if(mean_CI){}## Code someday possibility to diplay mean process' CI
    
  return(gg)
}

plot_animate = function(pred_gp, data = NULL, data_train = NULL, mean = NULL, mean_CI = F, file = "gganim.gif")
{ ## pred_gp : tibble coming out of the pred_gp_animate() function, columns required : 'Timestamp', 'Mean', 'Var'
  ## data : tibble of observational data, columns required : 'Timestamp', 'Output' (Optional)
  ####
  ## return : plot the animated curves of the GP with the 0.95 confidence interval (optional display raw data)

  gg = plot_gp_clust(pred_gp, data, data_train, mean, mean_CI) +
    transition_states(Nb_data, transition_length = 2, state_length = 1)
  animate(gg, renderer = gifski_renderer(file)) %>% return()
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


############### SIMULATED DATA ######################
simu_indiv = function(ID, t, kern = kernel_mu, theta, mean, var)
{ ## Return Age and Performance of a simulated individual
  ## nb : Number of timestamps
  ## Id : identifier of the individual
  ## tmin : Value of the first timestamp
  ## tmax : Value of the last timestamp
  ## a : Multiplying coefficient
  ## b : Intercept
  ## var : Variance of the additive noise
  # t = seq(tmin, tmax, length.out = nb)
  # indiv = tibble('ID' = rep(as.character(ID), nb), 'Timestamp' = t, 'Input' = paste0('X', t),
  #                'Output' = a*t + b + rnorm(nb,0,var),'a' = rep(runif(1,0,5) %>% round(2), nb),
  #                'b' = rep(runif(1,0,0.5) %>% round(2), nb))
  
  db = tibble('ID' = ID,
              'Timestamp' = t,
              'Input' = paste0('X', t),
              'Output' = rmvnorm(1, rep(mean,length(t)), kern_to_cov(t, kern, theta, sigma = var)) %>% as.vector())
  return(db)
}

#set.seed(42)
M = 10
N = 10
t = matrix(0, ncol = N, nrow = M)
for(i in 1:M){t[i,] = sample(seq(10, 20, 0.01),N, replace = F) %>% sort()}

db_train = simu_indiv(ID = '1', t[1,], kernel_mu, theta = c(2,1), mean = 45, var = 0.2)
for(i in 2:M)
{
  k = i %% 5
  if(k == 0){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(2,2), mean = 45, var = 0.6))}
  if(k == 1){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(2,1), mean = 40, var = 0.2))}
  if(k == 2){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(3,2), mean = 50, var = 0.3))}
  if(k == 3){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(1,2), mean = 35, var = 0.4))}
  if(k == 4){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(1,1), mean = 55, var = 0.5))}
}


db_obs = simu_indiv(ID = (M+1) %>% as.character(), sample(seq(10, 20, 0.01), N, replace = F) %>%
                                   sort(), kernel_mu, theta = c(2,1), mean = 45, var = 0.2)

# ################ INITIALISATION ######################
# ini_hp = list('theta_k' = c(1,1), 'theta_i' = c(1, 1, 0.2))
# 
#### TEST ####
# k = seq_len(2)
# tau_i_k_test = replicate(length(k), rep(1,length(unique(db_train$ID)))) %>%
#   apply(1,function(x) x / sum(x)) %>%
#   `rownames<-`(paste0('K', k)) %>%
#   `colnames<-`(unique(db_train$ID)) %>%
#   apply(1, as.list)
# 
# prior_mean_k = list('K1' = 42, 'K2' = 50)
# ini_hp_test = list('theta_k' = c(2, 0.5, 0.1), 'theta_i' = c(1, 1, 0.2))
# common_hp_k = T
# common_hp_i = T
# # Training
# t1 = Sys.time()
# training_test = training_VEM(db_train, prior_mean_k, ini_hp_test, kernel_mu, kernel,
#                              tau_i_k_test, common_hp_k, common_hp_i)
# t2 = Sys.time()
# c_time = (t2-t1) %>% print()
# pi_k_test = lapply(training_test$param$tau_i_k, function(x) Reduce("+", x)/ length(x)) 
# ## Posterior mu_k
# timestamps = seq(10, 20, 0.01)
# post_test = posterior_mu_k(db_train, timestamps, prior_mean_k, kernel_mu, kernel, training_test)
# ## Pred GP
# pred_gp_test = pred_gp_clust(db_obs[1:5,], timestamps, post_test, kern = kernel,
#              list('theta' = training_test$theta_i$`1`[1:2], 'sigma' = training_test$theta_i$`1`[[3]],
#                   'tau_k' = pi_k_test ))
# # Plot GP
# plot_gp_clust(pred_gp_test, data = db_obs, data_train = db_train, mean = post_test$mean)
# # Plot GIF
# pred_gif = pred_gp_clust_animate(db_obs, timestamps, post_test, kern = kernel,
#                       list('theta' = training_test$theta_i$`1`[1:2], 'sigma' = training_test$theta_i$`1`[[3]],
#                            'tau_k' = pi_k_test ))
# plot_animate(pred_gif, data = db_obs, data_train = db_train, mean = post_test$mean, file = "bla.gif")
