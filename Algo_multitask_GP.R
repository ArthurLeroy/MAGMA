library(tidyverse)
library(MASS)
library(Matrix)
library(mvtnorm)
library(optimr)

#setwd(dir = 'C:/Users/user/CloudStation/Maths/These/Processus Gaussiens/Code R/Algo multitask GP')
source('Computing_functions.R')

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
for(i in 1:M){t[i,] = sample(seq(10, 20, 0.5),N, replace = F) %>% sort()}

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


db_obs = simu_indiv(ID = (M+1) %>% as.character(), sample(seq(10, 20, 0.5), N, replace = F) %>% sort(),
                    kernel_mu, theta = c(2,1), mean = 40, var = 0.2)

# ################ INITIALISATION ######################
ini_hp = list('theta_0' = c(1,1), 'theta_i' = c(1, 1, 0.2))
m_0 = 50

################ TRAINING FUNCTIONS ################## 
training = function(db, prior_mean, ini_hp, kern_0, kern_i, common_hp = T)
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
  logLL_monitoring = - Inf
  
  for(i in 1:n_loop_max)
  { 
    print(i)
    ## E-Step
    param = e_step(db, prior_mean, kern_0, kern_i, hp)   
    
    ## Monitoring of the LL
    new_logLL_monitoring = logL_monitoring(hp, db, kern_i, kern_0, param$mean, param$cov, m_0 = prior_mean) - 
                            0.5 * ( nrow(param$cov) + log(det(param$cov)) ) 
    
    c(new_logLL_monitoring) %>% print()
    
    diff_moni = new_logLL_monitoring - logLL_monitoring
    #if(diff_moni < 0){stop('Likelihood descreased')}
    
    ## M-Step
    new_hp = m_step(db, hp, mean = param$mean, cov = param$cov, kern_0, kern_i, prior_mean, common_hp)

    ## Testing the stoping condition
    logL_new = logL_monitoring(new_hp, db, kern_i, kern_0, param$mean, param$cov, m_0 = prior_mean)
    eps = (logL_new - logL_monitoring(hp, db, kern_i, kern_0, param$mean, param$cov, m_0 = prior_mean)) / 
          abs(logL_new)

    if(eps > 0 & eps < 1e-10)
    {
      cv = 'TRUE'
      break
    }
    hp = new_hp
    logLL_monitoring = new_logLL_monitoring
  }
  return(list('theta_0' = new_hp$theta_0, 'theta_i' = new_hp$theta_i, 'convergence' = cv, 'param' = param))
}


train_new_gp = function(db, mean_mu, cov_mu, ini_hp, kern_i)
{
  mean = mean_mu %>% filter(Timestamp %in% db$Timestamp) %>% pull(Output)
  new_cov = cov_mu[paste0('X', db$Timestamp), paste0('X', db$Timestamp)]
  LL_GP<- function(hp, db, kern) 
  {
    return(-dmvnorm(db$Output, mean, solve(kern_to_cov(db$Timestamp, kern, theta = hp[1:2], sigma = hp[3]) + new_cov),
                    log = T))
  }
  new_hp = opm(ini_hp, fn = LL_GP, gr = gr_one_GP, kern = kern_i, db = db, method = "L-BFGS-B", control = list(kkt = FALSE)) 
  list('theta' = new_hp[1:2], 'sigma' = new_hp[3]) %>% return()
}

# t1 = Sys.time()
# bla = train_new_gp(db_obs, rep(0, nrow(db_obs)), 
#                    cov_mu = kern_to_cov(db_obs$Timestamp, theta = c(2,1), sigma = 0.5), 
#                    ini_hp$theta_i, kernel)
# t2 = Sys.time()
# t = t2 - t1


################ PREDICTION FUNCTIONS ################
posterior_mu = function(db, timestamps, m_0, kern_0, kern_i, hp)
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

pred_gp = function(db, timestamps, mean_mu = NULL , cov_mu = NULL, 
                   kern = kernel, theta = list(1,0.2), sigma = 0.2)
{ ## db: tibble of data columns required ('Timestamp', 'Input', 'Output')(Input format : paste0('X', Timestamp))
  ## timestamps : timestamps on which we want a prediction
  ## mean_mu : mean value of mean GP at timestamps (obs + pred) (matrix dim: timestamps x 1, with Input rownames)
  ## cov_mu : covariance of mean GP at timestamps (obs + pred) (square matrix, with Input row/colnames)
  ## kernel : kernel associated to the covariance function of the GP
  ## theta : list of hyperparameters for the kernel of the GP
  ## sigma : variance of the error term of the models
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

full_algo = function(db, new_db, timestamps, kern_i, plot = T, prior_mean,
                     kern_0 = NULL, list_hp = NULL, mu = NULL, ini_hp = NULL, hp_new_i = NULL)
{ ## db : Database containing all training data from all individuals. Column: ID - Timestamp - Output.
  ## new_db : Database containing data for a new individual we want a prediction on.
  ## timestamps : Timestamps we want to predict at.
  ## kern_i : Kernel associated to individual GPs.
  ## plot : Boolean indicating whether we want to display a graph at the end of computations.
  ## prior_mean : Prior arbitrary value for the mean process. Optional, not needed if 'mu' is given.
  ## kern_0 : Kernel associated to the mean GPs. Optional, not needed if 'mu' is given.
  ## list_hp : Hyper-parameters for all individuals in training set. Optional, computed if NULL.
  ## mu : Database containing parameters of the mean GP at all prediction timestamps. Optional, computed if NULL.
  ## ini_hp : Initial values of the HP to start the training. Optional, not needed if 'list_hp' is given.
  ## hp_new_i : Hyper-pameters for the new individual to predict. Optional, computed if NULL. 
  ####
  ## return : predicted GP parameters | posterior mean process | all trained hyperparameters
  
  if(is.null(list_hp)){list_hp = training(db, prior_mean, ini_hp, kern_0, kern_i)[c('theta_0', 'theta_i')]}
  
  t_pred = timestamps %>% union(unique(db$Timestamp)) %>% union(unique(new_db$Timestamp)) %>% sort()
  if(is.null(mu)){mu = posterior_mu(db, t_pred, prior_mean, kern_0, kern_i, list_hp)}
  mean_mu = mu$Mean
  t_mu = mu$Timestamp
  obs_input = new_db$Input
  
  if(is.null(hp_new_i)){hp_new_i = train_new_gp(new_db, mean_mu = mean_mu[obs_input,], 
                                                cov_mu = mu$Cov[obs_input,obs_input], ini_hp$theta_i, kern_i)}
  
  pred = pred_gp(new_db, timestamps, mean_mu = mean_mu, cov_mu = mu$Cov, kern_i, hp_new_i$theta, hp_new_i$sigma)
  
  if(plot){(pred %>% plot_gp(data = rbind(new_db, db)) + geom_point(aes(t_mu, mean_mu))) %>% print()}

  list('Prediction' = pred, 'Mean_process' = mu, 'Hyperparameters' = list('other_i' = list_hp, 'new_i' = hp_new_i)) %>% 
  return()
}


##### Testing codes ############

### Testing the full_algo function
# bla = training(db_train, 0, ini_hp, kernel_mu, kernel)
#  list_hp_test = bla[c('theta_0','theta_i')]
# list_hp_test$theta_0 = c(7,1) ## Because hp of K_0 are crucial and often mistrained
# blab = full_algo(db_train, db_obs[3:7,], seq(10, 20, 0.05), kernel, plot = T, prior_mean = 0, kern_0 = kernel_mu,
#           list_hp = list_hp_test, mu = NULL, ini_hp = ini_hp, hp_new_i = hp_one_gp)
# 
# plot_gp(blab$Prediction, db_obs[3:7,])

# hp_one_gp = list('theta' = c(5,2), 'sigma' = 0.2)
# fu = pred_gp(db_obs[3:7,], seq(10, 20, 0.05), mean_mu = 0 , cov_mu = NULL, kernel,
#              theta = c(5,2), 0.2)
# plot_gp(fu, db_obs[3:7,])
## Testing the training function
# bla = training(db_train, 45, ini_hp, kernel_mu, kernel, common_hp = T)
# fu = bla$param$mean$Output
# names(fu) = paste0('X', bla$param$mean$Timestamp)
# hp_pred = train_new_gp(db_obs, bla$param$mean, bla$param$cov, ini_hp$theta_i, kernel)
# 
# pred_gp(db_obs[1:3,], timestamps = bla$param$mean$Timestamp, mean_mu = fu %>% as.matrix(),
#         cov_mu = bla$param$cov, theta = hp_pred$theta, sigma = hp_pred$sigma) %>%
#   plot_gp(data = rbind(db_obs[1:10,], db_train)) + geom_point(aes(bla$param$mean$Timestamp, bla$param$mean$Output))

# ### Testing update_mean and update_inv
# grid = seq(0,101,0.619)
# pred = c()
# for(i in grid){
#   grid_t = c(i,1:100)
#   list_invi = kern_to_inv(db_train, theta = param_test$theta_i, sigma = 0)
#   list_valuei = list()
# 
#   for(i in unique(db_train$ID)){list_valuei[[i]] =  db_train %>% filter(ID == i) %>% pull(Output)}
#   upinv = update_inv(kern_to_inv(grid_t, theta = c(2, 1), sigma = 0.2), list_invi)
#   upcov = upinv %>% solve()
#   upmean = (upinv %>% solve()) %*% update_mean(40, kern_to_inv(grid_t, theta = c(2, 1), sigma = 0.2), list_invi , list_valuei)
# 
#   pred = rbind(pred, pred_gp(data = db_obs %>% filter(Timestamp %in% c(1:40)), timestamps = grid_t, mean_mu = upmean,
#                              cov_mu = upcov, kern = kernel, theta = list(2, 1), sigma = 0.2)[1,])
# }
# 
# plot_gp(pred, data = rbind(db_obs, db_train))
