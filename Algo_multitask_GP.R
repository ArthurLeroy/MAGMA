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
plot_db = function(data = db)
{
  ggplot(db) + geom_smooth(aes(Timestamp, Output, color = ID)) + geom_point(aes(Timestamp, Output, color = ID)) 
}

################ SIMULATED DATA ######################
simu_indiv = function(nb, ID, tmin, tmax, a, b, var)
{ ## Return Age and Performance of a simulated individual
  ## nb : Number of timestamps
  ## Id : identifier of the individual
  ## tmin : Value of the first timestamp
  ## tmax : Value of the last timestamp
  ## a : Multiplying coefficient
  ## b : Intercept
  ## var : Variance of the additive noise
  t = seq(tmin, tmax, length.out = nb)
  indiv = tibble('ID' = rep(as.character(ID), nb), 'Timestamp' = t, 'Input' = paste0('X', t),
                 'Output' = a*t + b + rnorm(nb,0,var),'a' = rep(runif(1,0,5) %>% round(2), nb), 
                 'b' = rep(runif(1,0,0.5) %>% round(2), nb))
  return(indiv)
}

set.seed(42)
db = rbind(simu_indiv(10, 1, 11, 20, -1, 60, 2),
           simu_indiv(10, 2, 11, 20, -0.2, 50, 2),
           simu_indiv(10, 3, 11, 20, -1.6, 70, 2),
           simu_indiv(10, 4, 11, 20, -0.8, 65, 2),
           simu_indiv(10, 5, 11, 20, -0.6, 55, 2))

db_obs = simu_indiv(10, 6, 11, 20, -1.4, 54, 2)

################ INITIALISATION ######################
theta_0_0 = 1
theta_i_0 = 1
sigma_0 = 1
timestamps = 11:20
m_0 = 50

################ TRAINING FUNCTIONS ################## 
training = function(data = db, prior_mean, theta_0_0 = 1, theta_i_0 = 1, sigma_0 = 1, kern = kernel)
{ ## data : database with all individuals in training set. Column required : 'ID', Timestamp', 'Input', 'Output'
  ## theta_i_0 : common initial value of HP for the common prior kernel of all individuals
  ## sigma_0 : initial value of the variance of the error process
  ## kern : kernel associated to common covariance functions of all individuals 
  ####
  ## return : list of trained HP, list of inv covariance matrices and weighted data points(useful in posterior_mu())
  
  list_ID = data %>% pull(ID) %>% unique()
  list_theta_i = list()
  list_weighted_values_i = list()
  list_inv_i = list()
  n_loop_max = 10
  
  for(i in list_ID)
  {
    theta_i = theta_i_0
    theta_0 = theta_0_0
    sigma = sigma_0
    y_i = data %>% filter(ID == i) %>% pull(Output)
    t_i = data %>% filter(ID == i) %>% pull(Timestamp)
    
    for(i in 1:n_loop_max)
    {
      param = e_step_GP(theta_i, theta_0, sigma)
      new_hp = m_step_GP(mean = param$mean, inv = param$inv)
      
      eps = logL(y = y_i, t = t_i, theta_0 = new_param$theta_0, theta_i = new_hp$theta_i, sigma = new_hp$sigma) - 
            logL(y= y_i, t = t_i , theta_0 = param$theta_0, theta_i = hp$theta_i , sigma = hp$sigma )
      if(eps > 0 & eps < 1e-5){break}
      if(eps > 0){hp = new_hp}
    }
  }
  return(list('trained_sigma' = sigma,'trained_theta_0' = theta_0,  'list_trained_theta_i' = list_theta_i, 
              'list_weigthed_values_i' = list_weighted_values_i, 'list_inv_i' = list_inv_i)) 
}

################ PREDICTION FUNCTIONS ################
posterior_mu = function(data, timestamps, prior_mean, kern, theta_0, list_inv_i, list_value_i)
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
  
  prior_inv = kern_to_inv(timestamps, kern = kern, theta = theta_0)
  weighted_values = update_mean(old_mean = old_mean, old_inv = prior_inv, list_value_i = list_value_i)
  
  inv_mu = update_inv(old_inv = prior_inv, list_inv_i = list_inv_i)
  cov_mu = inv_mu %>% solve()
  mean_mu = cov_mu %*% weighted_values
  
  return(list('mean' = mean_mu , 'cov' = cov_mu))
}

pred_gp = function(data = db_obs, timestamps = 10:21, mean_mu = NULL , cov_mu = NULL, 
                   kern = kernel, theta = list(1,0.2, 0))
{ ## data: tibble of data columns required ('Timestamp', 'Input', 'Output')(Input format : paste0('X', Timestamp))
  ## timestamps : timestamps on which we want a prediction
  ## mean_mu : mean value of mean GP at timestamps (obs + pred) (matrix dim: timestamps x 1, with Input rownames)
  ## cov_mu : covariance of mean GP at timestamps (obs + pred) (square matrix, with Input row/colnames)
  ## kernel : kernel associated to the covariance function of the GP
  ## theta : list of hyperparameters for the kernel of the GP
  ## sigma : variance of the error term of the model
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
  
  inv_mat = (kern_to_cov(tn, kern,  theta) + cov_mu[input, input]) %>% solve()
  cov_tn_t = sapply(timestamps, function(x) kern(x,tn, theta)) + cov_mu[input,input_t]
  cov_t_t = sapply(timestamps, function(x) kern(x,timestamps, theta)) + cov_mu[input_t ,input_t]
           
  pred = tibble('Timestamp' = timestamps,
                'Mean' = mean_mu[input_t,] + t(cov_tn_t) %*% inv_mat %*% (yn - mean_mu[input,]) %>% as.vector,
                'Var' = (cov_t_t - t(cov_tn_t) %*% inv_mat %*% cov_tn_t) %>% diag())
  return(pred)
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
        theta = c(3,0.4,0.2)) %>%  plot_gp(data = db_obs)


full_algo = function()
{
  
}


##### Testing codes ############
grid = seq(9,21,0.019)
pred = c()
for(i in grid){
  grid_t = c(i,11:20)
  list_invi = kern_to_inv(db)
  list_valuei = list()
  for(i in unique(db$ID)){list_valuei[[i]] =  db %>% filter(ID == i) %>% pull(Output)}
  upinv = update_inv(kern_to_inv(grid_t, theta = c(14, 1, 0.2)), list_invi) 
  upcov = upinv %>% solve()
  upmean = (upinv %>% solve()) %*% update_mean(40, kern_to_inv(grid_t, theta = c(14, 1, 0.2)), list_invi , list_valuei)
  
  pred = rbind(pred, pred_gp(data = db_obs %>% filter(Timestamp %in% c(11:16)), timestamps = grid_t, mean_mu = upmean,
                             cov_mu = upcov, kern = kernel, theta = list(5, 1, 0))[1,])
}

plot_gp(pred, data = db_obs)
