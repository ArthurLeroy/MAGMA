library(tidyverse)
library(MASS)
library(Matrix)

setwd(dir = 'C:/Users/user/CloudStation/Maths/These/Processus Gaussiens/Code R/Algo multitask GP')
source('Computing_functions.R')

################ DATA IMPORT AND CLEANING #########

#raw_data = read.csv2("E:/Users/Arthur.Royaldinateur/CloudStation/Articles/Clustering courbes natation 07_2018/Data/Swimming data/100 Nage Libre _50_Homme_id.csv", stringsAsFactors = F)
raw_data = read_csv2("C:/Users/user/CloudStation/Articles/Clustering courbes natation 07_2018/Data/Swimming data/100 Nage Libre _50_Homme_id.csv")

db = raw_data %>% filter(ID <= 50) %>% arrange(ID, Age) %>% distinct(Age, .keep_all= T)

db$ID = as.character(db$ID)
ID1 = db %>% filter(ID == 5) %>% arrange(Age)

## Plot raw curves of several individuals
plot_db = function(db = db_train)
{
  ggplot(db) + geom_smooth(aes(Timestamp, Output, color = ID)) + geom_point(aes(Timestamp, Output, color = ID)) 
}

################ SIMULATED DATA ######################
simu_indiv = function(ID, t, kern, theta, mean, var)
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

set.seed(42)
t = c(11:20, 24)
db_train = rbind(simu_indiv(ID = '1', t, kernel, theta = c(2,1), mean = 45, var = 0.2),
           simu_indiv(ID = '2', t, kernel, theta = c(1,1), mean = 45, var = 0.3),
           simu_indiv(ID = '3', t, kernel, theta = c(1,2), mean = 45, var = 0.4),
           simu_indiv(ID = '4', t, kernel, theta = c(1.5,1), mean = 45, var = 0.5),
           simu_indiv(ID = '5', t, kernel, theta = c(3,2), mean = 45, var = 0.6))

db_obs = simu_indiv(ID = '6', t, kernel, theta = c(2,1), mean = 30, var = 0.2)

################ INITIALISATION ######################
ini_param = list('theta_0' = c(1,1), 'theta_i' = c(1, 1, 0.2))
timestamps = 11:20
m_0 = 50

################ TRAINING FUNCTIONS ################## 
training = function(db, prior_mean, ini_param, kern_0 = kernel_mu, kern_i = kernel)
{ ## db : database with all individuals in training set. Column required : 'ID', Timestamp', 'Input', 'Output'
  ## prior_mean : prior mean parameter of the mean GP (mu_0)
  ## ini_param : initial values of HP for the kernels
  ## kern_0 : kernel associated to covariance functions of the mean GP
  ## kern_i : kernel associated to common covariance functions of all individuals GPs
  ####
  ## return : list of trained HP, boolean to indicate convergence
  
  n_loop_max = 25
  list_ID = unique(db$ID)
  hp = list('theta_0' = ini_param$theta_0, 
            'theta_i' = ini_param$theta_i %>% setNames(nm = list_ID))
  cv = 'FALSE'
  
  for(i in 1:n_loop_max)
  {
    param = e_step(db, prior_mean, kern_0, kern_i, hp)   
    print(i)
    new_hp = m_step(db, hp, mean = param$mean, cov = param$cov, kern_0, kern_i, prior_mean)

    #logL_multi_GP(new_hp, db, kern_i, kern_0, param$mean, param$cov, m_0 = prior_mean) %>% print()
    
    logL_new = logL_multi_GP(new_hp, db, kern_i, kern_0, param$mean, param$cov, m_0 = prior_mean)
    eps = (logL_new - logL_multi_GP(hp, db, kern_i, kern_0, param$mean, param$cov, m_0 = prior_mean)) / 
          abs(logL_new)
    #print(eps)
    print(param$cov %>% det())
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

################ PREDICTION FUNCTIONS ################
posterior_mu = function(data, timestamps, prior_mean, kern, theta_0, sigma_0, list_inv_i, list_value_i)
{ ## data : matrix of data columns required ('Timestamp', 'Input', 'Output')(Input format : paste0('X', Timestamp))
  ## timestamps : timestamps on which we want a prediction
  ## prior_mean : prior mean value of the mean GP (scalar value or vector of same length as 'timestamps')
  ## kern : kernel associated to the covariance function of the mean GP
  ## theta_0 : hyperparameters for the kernel of the mean GP
  ## list_inv_i = list of inverse of covariance matrices for all observed individuals
  ## list_value_i = list of (inverse of covariance matrices) %*% (observed output) for all observed individuals
  ## The last two lists are supposed to be retrieved from the training to avoid recomputation of covariance matrices
  ####
  ## return : pamameters of the mean GP at timestamps
  
  
  if(length(prior_mean) == 1){prior_mean = rep(prior_mean, length(timestamps))}
  old_mean = matrix(prior_mean, ncol = 1, nrow = length(timestamps), dimnames = list(paste0('X', timestamps)))
  
  prior_inv = kern_to_inv(timestamps, kern, theta = theta_0, sigma = sigma_0)
  weighted_values = update_mean(old_mean = old_mean, old_inv = prior_inv, list_value_i = list_value_i)
  
  inv_mu = update_inv(old_inv = prior_inv, list_inv_i = list_inv_i)
  cov_mu = inv_mu %>% solve()
  mean_mu = cov_mu %*% weighted_values
  
  return(list('mean' = mean_mu , 'cov' = cov_mu))
}

pred_gp = function(data = db_obs, timestamps = 10:21, mean_mu = NULL , cov_mu = NULL, 
                   kern = kernel, theta = list(1,0.2), sigma = 0.2) 
{ ## data: tibble of data columns required ('Timestamp', 'Input', 'Output')(Input format : paste0('X', Timestamp))
  ## timestamps : timestamps on which we want a prediction
  ## mean_mu : mean value of mean GP at timestamps (obs + pred) (matrix dim: timestamps x 1, with Input rownames)
  ## cov_mu : covariance of mean GP at timestamps (obs + pred) (square matrix, with Input row/colnames)
  ## kernel : kernel associated to the covariance function of the GP
  ## theta : list of hyperparameters for the kernel of the GP
  ## sigma : variance of the error term of the models
  ####
  ## return : pamameters of the gaussian density predicted at timestamps 
  tn = data %>% pull(Timestamp)
  input = data %>% pull(Input)
  yn = data %>% pull(Output)
  
  input_t = paste0('X', timestamps)
  all_times = c(tn,timestamps)
  
  if(is.null(cov_mu))
  {
    cov_mu = matrix(0, length(all_times), length(all_times),
                    dimnames = list(paste0('X', all_times), paste0('X', all_times)))
  }
  if(is.null(mean_mu)){mean_mu = matrix(0, length(all_times), 1, dimnames = list(paste0('X', all_times)))}
  if(length(mean_mu) == 1){mean_mu = matrix(mean_mu, length(all_times), 1, dimnames = list(paste0('X', all_times)))}
  
  inv_mat = (kern_to_cov(tn, kern, theta, sigma) + cov_mu[input, input]) %>% solve()
  cov_tn_t = sapply(timestamps, function(t) sapply(tn, function(s) kern(t, s, theta))) + cov_mu[input,input_t]
  cov_t_t = sapply(timestamps, function(t) sapply(timestamps, function(s) kern(t, s, theta))) + cov_mu[input_t ,input_t]
  
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
  
  if(!is.null(data)){gg = gg + geom_point(data = data, aes(x = Timestamp, y = Output), shape = 4, color = 'red')}
  return(gg)
}

################ APPLICATION #########################
## Test graph
pred_gp(db_obs %>% filter(Timestamp %in% 11:20), timestamps = seq(10,20, 0.03), mean_mu = 50,
        theta = c(3,0.4), sigma = 0.2) %>% plot_gp(data = db_obs)


full_algo = function(db, timestamps, ini_param, prior_mean, kern_0, kern_i)
{
  
}


##### Testing codes ############
bla = training(db_train, 0, param_test, kernel_mu, kernel)
fu = bla$param$mean$Output
names(fu) = colnames(bla$param$cov)

pred_gp(db_obs %>% filter(Timestamp %in% 11:15), timestamps = c(11:20, 24), mean_mu = fu %>% as.matrix(),
        cov_mu = bla$param$cov, theta = c(1,1), sigma = 0.5) %>% plot_gp(data = rbind(db_obs, db_train))

### Testing update_mean and update_inv
# grid = seq(9,21,0.019)
# pred = c()
# for(i in grid){
#   grid_t = c(i,11:20)
#   list_invi = kern_to_inv(db_train)
#   list_valuei = list()
# 
#   for(i in unique(db_train$ID)){list_valuei[[i]] =  db_train %>% filter(ID == i) %>% pull(Output)}
#   upinv = update_inv(kern_to_inv(grid_t, theta = c(2, 1), sigma = 0.2), list_invi)
#   upcov = upinv %>% solve()
#   upmean = (upinv %>% solve()) %*% update_mean(40, kern_to_inv(grid_t, theta = c(2, 1), sigma = 0.2), list_invi , list_valuei)
# 
#   pred = rbind(pred, pred_gp(data = db_obs %>% filter(Timestamp %in% c(11:15)), timestamps = grid_t, mean_mu = upmean,
#                              cov_mu = upcov, kern = kernel, theta = list(2, 1), sigma = 0.2)[1,])
# }
# 
# plot_gp(pred, data = rbind(db_obs, db_train))
