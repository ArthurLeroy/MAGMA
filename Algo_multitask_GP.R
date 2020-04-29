library(tidyverse)
library(MASS)
library(Matrix)
library(optimr)
library(mvtnorm)
library(plotly)
library(gganimate)
library(transformr)
library(gifski)
library(png)


setwd(dir = 'C:/Users/user/CloudStation/Maths/These/Processus Gaussiens/Code R/Algo multitask GP')
source('Computing_functions.R')

################ TRAINING FUNCTIONS ################## 
training = function(db, prior_mean, ini_hp, kern_0, kern_i, common_hp = T)
{ ## db : database with all individuals in training set. Column required : 'ID', Timestamp', 'Output'
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
  list_plot = list()
  t1 = Sys.time()
  for(i in 1:n_loop_max)
  { 
    print(i)
    ## E-Step
    param = e_step(db, prior_mean, kern_0, kern_i, hp)   
    
    ##list_plot[[i]] = param$pred_GP %>% plot_gp(data_train = db)
    ## Return list_plot if you want monitoring graphs of the mean process' learning
    
    ## M-Step
    new_hp = m_step(db, hp, mean = param$mean, cov = param$cov, kern_0, kern_i, prior_mean, common_hp)
  
    ## Testing the stoping condition
    ## Monitoring of the LL
    new_logLL_monitoring = logL_monitoring(hp, db, kern_i, kern_0, param$mean, param$cov, prior_mean) 
                            # + 0.5 * log(det(param$cov)) for an exact likelihood but constant respectively to HPs
    c('logLL = ', new_logLL_monitoring) %>% print()
    diff_moni = new_logLL_monitoring - logLL_monitoring
    if(diff_moni < - 0.1){warning('Likelihood descreased')}

    logL_new_hp = logL_monitoring(new_hp, db, kern_i, kern_0, param$mean, param$cov, prior_mean)
                  # + 0.5 * log(det(param$cov)) for an exact likelihood but constant respectively to HPs
    eps = (logL_new_hp - new_logLL_monitoring) / abs(logL_new_hp)
    c('eps =', eps) %>% print()
    
    if(eps > 0 & eps < 1e-3)
    {
      cv = 'TRUE'
      break
    }
    hp = new_hp
    logLL_monitoring = new_logLL_monitoring
  }
  t2 = Sys.time()
  list('hp' = new_hp, 'convergence' = cv, 'param' = param, 'Time_train' = as.numeric(t2-t1)) %>% 
  return()
}

train_new_gp = function(db, mean_mu, cov_mu, ini_hp, kern_i)
{
  if(is.vector(mean_mu)){mean = mean_mu} 
  else {mean = mean_mu %>% filter(Timestamp %in% db$Timestamp) %>% pull(Output) %>% as.vector}
  
  if(is.matrix(cov_mu)){new_cov = cov_mu[paste0('X', db$Timestamp), paste0('X', db$Timestamp)]}
  else {new_cov = 0}
  
  LL_GP<- function(hp, db, kern) 
  {
    return(-dmvnorm(db$Output, mean, solve(kern_to_cov(db$Timestamp, kern, theta = hp[1:2], sigma = hp[3]) + new_cov),
                    log = T))
  }
  new_hp = opm(ini_hp, fn = LL_GP, gr = NULL, kern = kern_i, db = db, method = "L-BFGS-B", control = list(kkt = FALSE)) 
  list('theta' = new_hp[1:2], 'sigma' = new_hp[3]) %>% return()
}

################ PREDICTION FUNCTIONS ################
posterior_mu = function(db, new_db, timestamps, m_0, kern_0, kern_i, hp)
{ ## db : matrix of data columns required ('Timestamp', 'Output')
  ## timestamps : timestamps on which we want a prediction
  ## prior_mean : prior mean value of the mean GP (scalar value or vector of same length as 'timestamps')
  ## kern_0 : kernel associated to the covariance function of the mean GP
  ####
  ## return : pamameters of the mean GP at timestamps
  t_pred = timestamps %>% union(unique(db$Timestamp)) %>% union(unique(new_db$Timestamp)) %>% sort()

  inv_0 = kern_to_inv(t_pred, kern_0, hp$theta_0, sigma = 0.01)
  inv_i = kern_to_inv(db, kern_i, hp$theta_i, sigma = 0)
  value_i = base::split(db$Output, list(db$ID))

  new_inv = update_inv(prior_inv = inv_0, list_inv_i = inv_i)
  new_cov = tryCatch(solve(new_inv), error = function(e){MASS::ginv(new_inv)}) ## fast or slow matrix inversion if singular
  rownames(new_cov) = rownames(new_inv)
  colnames(new_cov) = colnames(new_inv)

  weighted_mean = update_mean(prior_mean = m_0, prior_inv = inv_0, list_inv_i = inv_i, list_value_i = value_i)
  new_mean = (new_cov %*% weighted_mean) %>% as.vector
  
  #names(mean_mu) = paste0('X', t_mu)
  list('mean' = tibble('Timestamp' = t_pred, 'Output' = new_mean) , 'cov' = new_cov, 
       'pred_GP' = tibble('Timestamp' = t_pred, 'Mean' = new_mean, 'Var' = diag(new_cov))) %>% return()
}

pred_gp = function(db, timestamps = NULL, mean_mu = 0, cov_mu = NULL, 
                   kern = kernel, theta = list(1,0.2), sigma = 0.2)
{ ## db: tibble of data columns required ('Timestamp', 'Output')
  ## timestamps : timestamps on which we want a prediction
  ## mean_mu : mean value of mean GP at timestamps (obs + pred) (matrix dim: timestamps x 1, with Input rownames)
  ## cov_mu : covariance of mean GP at timestamps (obs + pred) (square matrix, with Input row/colnames)
  ## kernel : kernel associated to the covariance function of the GP
  ## theta : list of hyperparameters for the kernel of the GP
  ## sigma : variance of the error term of the models
  ####
  ## return : pamameters of the gaussian density predicted at timestamps 
  tn = db %>% pull(Timestamp)
  #input = db %>% pull(Input)
  input = paste0('X', db$Timestamp)
  yn = db %>% pull(Output)
  
  ## Define a default prediction grid
  if(is.null(timestamps)){timestamps = seq(min(tn), max(tn), length.out = 500)}
  input_t = paste0('X', timestamps)
  all_times = union(tn,timestamps)
  
  if(is.null(cov_mu))
  { ## Case of standard GP regression. Without trained posterior mean process. 
    cov_mu = matrix(0, length(all_times), length(all_times),
                    dimnames = list(paste0('X', all_times), paste0('X', all_times)))
  }
  if(length(mean_mu) == 1)
  { ## If the provided mean is a constant function
    mean_mu_obs = rep(mean_mu, length(tn))
    mean_mu_pred = rep(mean_mu, length(timestamps))
  }
  else
  { ## If the provided mean has defined values at timestamps, typically from training. Format : Timestamp, Output
    mean_mu_obs = mean_mu %>%  filter(Timestamp %in% tn) %>% pull(Output)
    mean_mu_pred = mean_mu %>% filter(Timestamp %in% timestamps) %>% pull(Output)
  }
  
  inv_mat = (kern_to_cov(tn, kern, theta, sigma) + cov_mu[input, input]) %>% solve()
  cov_tn_t = kern(mat_dist(tn, timestamps), theta) + cov_mu[input,input_t]

  cov_t_t = kern_to_cov(timestamps, kern, theta, sigma) + cov_mu[input_t ,input_t]

  tibble('Timestamp' = timestamps,
         'Mean' = (mean_mu_pred + t(cov_tn_t) %*% inv_mat %*% (yn - mean_mu_obs)) %>% as.vector(),
         'Var' = (cov_t_t - t(cov_tn_t) %*% inv_mat %*% cov_tn_t) %>% diag()) %>% return()
}

pred_gp_animate = function(db, timestamps = NULL, mean_mu = 0, cov_mu = NULL, 
                           kern = kernel, theta = list(1,0.2), sigma = 0.2)
{
  ## Inputs : same as for a classic GP prediction
  ####
  ## return : tibble of classic GP predictions but with an inscreasing number of data points considered as 'observed'
  db %>% arrange(Timestamp)
  all_pred = tibble()

  if(is.null(timestamps)){timestamps = seq(min(db$Timestamp), max(db$Timestamp), length.out = 500)}
  
  for(j in 1:nrow(db))
  {
    pred_j = pred_gp(db[1:j,], timestamps, mean_mu, cov_mu, kern = kernel, theta, sigma) %>% mutate(Nb_data = j)
    all_pred = all_pred %>% rbind(pred_j)
  }
  return(all_pred)
}

################ PLOT FUNCTIONS ######################
plot_gp = function(pred_gp, data = NULL, data_train = NULL, mean = NULL, mean_CI = F)
{ ## pred_gp : tibble coming out of the pred_gp() function, columns required : 'Timestamp', 'Mean', 'Var'
  ## data : tibble of observational data, columns required : 'Timestamp', 'Output' (Optional)
  ####
  ## return : plot the predicted curve of the GP with the 0.95 confidence interval (optional display raw data)
  gg = ggplot() +
         geom_line(data = pred_gp, aes(x = Timestamp, y = Mean), color = 'blue') +
         geom_ribbon(data = pred_gp, aes(x = Timestamp, ymin = Mean - 1.96 * sqrt(Var), 
                                         ymax = Mean +  1.96 * sqrt(Var)), alpha = 0.2) + ylab('Output')
      
  ## Display the raw data and/or mean (with or without its CI) if provided
  if(!is.null(data_train)){gg = gg + geom_point(data = data_train, aes(x = Timestamp, y = Output, col = ID), shape = 4)}
  if(!is.null(data)){gg = gg + geom_point(data = data, aes(x = Timestamp, y = Output), size = 2, shape = 18)}
  if(!is.null(mean)){gg = gg + geom_line(data = mean, aes(x = Timestamp, y = Mean), linetype = 'dashed')}
  if(mean_CI){gg = gg + geom_ribbon(data = mean, aes(x = Timestamp, ymin = Mean - 1.96 * sqrt(Var), 
                                           ymax = Mean +  1.96 * sqrt(Var)), alpha = 0.4)}
  
  return(gg)
}

plot_db = function(db = db_train)
{ ## Visualize smoothed raw data. Format : ID, Timestamp, Output
  ggplot(db) + geom_smooth(aes(Timestamp, Output, color = ID)) + geom_point(aes(Timestamp, Output, color = ID))
}

plot_heat =  function(pred_gp, data = NULL, data_train = NULL, mean = NULL, ygrid = NULL, interactive = F, CI = T)
{ ## pred_gp : tibble coming out of the pred_gp() function, columns required : 'Timestamp', 'Mean', 'Var'
  ## data : tibble of observational data for the new individual, columns required : 'Timestamp', 'Output' (optional)
  ## data_train : tibble of training data for the model, , columns required : 'Timestamp', 'Output' (optional)
  ## mean : (optional)
  ## ygrid : grid on y-axis to plot the heatmap on (optional, default = range of data +- 2 sd)
  ## interactive : boolean indicating whether the output plot should be interactive (plotly)
  ## CI : boolean. If True, display the heatmap of Credible Intervals. If False, display the heatmap of likelihood
  ####
  ## return : An heatmap illustrating the GP predition results (Optionnal diplay of raw data and mean process)

  if(is.null(ygrid)){ygrid = seq(min(pred_gp$Mean) - 2 * sqrt(max(pred_gp$Var)) ,
                                 max(pred_gp$Mean) + 2 * sqrt(max(pred_gp$Var)), 0.1)}
  if(CI)
  {
    db_heat = tidyr::expand(pred_gp, nesting(Timestamp, Mean, Var), Ygrid = ygrid) %>% 
              mutate(Proba = 2 * pnorm(abs((Ygrid - Mean)/ sqrt(Var))) - 1) 
    gg = ggplot(db_heat) + geom_tile(aes(Timestamp, Ygrid, fill = Proba)) + scale_fill_distiller(palette = "RdPu") +
      theme_minimal() + ylab("Output") + labs(fill = "Proba CI")
  }
  else
  {
    db_heat = tidyr::expand(pred_gp, nesting(Timestamp, Mean, Var), Ygrid = ygrid) %>% 
              mutate(Proba = dnorm(Ygrid, mean = Mean, sd = sqrt(Var)) ) 
    gg = ggplot(db_heat) + geom_tile(aes(Timestamp, Ygrid, fill = Proba)) + 
         scale_fill_distiller(palette = "RdPu", trans = "reverse") + 
         theme_minimal() + ylab("Output") + labs(fill = "Likelihood")
  }
  ## Display data/training data/mean process if provided
  if(!is.null(data_train)){gg = gg + geom_point(data = data_train, aes(x = Timestamp, y = Output, col = ID), shape = 4)}
  if(!is.null(data)){gg = gg + geom_point(data = data, aes(x = Timestamp, y = Output), size = 2, shape = 18)}
  if(!is.null(mean)){gg = gg + geom_line(data = mean, aes(x = Timestamp, y = Mean), linetype = 'dashed')}

  ## Turn into an interactive plotly (html plot)
  if(interactive){gg = ggplotly(gg)}
  
  return(gg)
}

plot_animate = function(pred_gp, data = NULL, data_train = NULL, mean = NULL, mean_CI = F, file = "gganim.gif")
{ ## pred_gp : tibble coming out of the pred_gp_animate() function, columns required : 'Timestamp', 'Mean', 'Var'
  ## data : tibble of observational data, columns required : 'Timestamp', 'Output' (Optional)
  ####
  ## return : plot the animated curves of the GP with the 0.95 confidence interval (optional display raw data)

  gg = plot_gp(pred_gp, data, data_train, mean, mean_CI)  +
       transition_states(Nb_data, transition_length = 2, state_length = 1)
  animate(gg, renderer = gifski_renderer(file)) %>% return()
}

################ APPLICATION #########################

full_algo = function(db, new_db, timestamps, kern_i, common_hp = T, plot = T, prior_mean,
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

  ## If hp are not provided, train the model
  if(is.null(list_hp)){list_hp = training(db, prior_mean, ini_hp, kern_0, kern_i, common_hp)$hp}
  
  ## If mean GP (mu_0) paramaters at prediction timestamps are not provided , compute them
  if(is.null(mu)){mu = posterior_mu(db, new_db, timestamps, prior_mean, kern_0, kern_i, list_hp)}
  
  ## If hyperparameters of the GP for the new individuals are not provided, learn them
  ## If hyperparameters are common across individuals by hypothesis, simply pick them up from the trained model
  if(is.null(hp_new_i) & common_hp){hp_new_i = list('theta' = list_hp$theta_i[[1]][1:2], 
                                                    'sigma' = list_hp$theta_i[[1]][[3]])}
  else if(is.null(hp_new_i)){hp_new_i = train_new_gp(new_db, mu$mean, mu$cov, ini_hp$theta_i, kern_i)}
  
  pred = pred_gp(new_db, timestamps, mean_mu = mu$mean, cov_mu = mu$cov, kern_i, hp_new_i$theta, hp_new_i$sigma)

  ## If True, display a plot of the resulting GP prediction 
  if(plot)
  {
    (plot_gp(pred, data = new_db, data_train = db) + 
     geom_line(data = mu$pred_GP, aes(x = Timestamp, y = Mean), color = 'black', linetype = 'dashed')) %>% print()
  }  

  list('Prediction' = pred, 'Mean_process' = mu, 'Hyperparameters' = list('other_i' = list_hp, 'new_i' = hp_new_i)) %>% 
  return()
}

############### SIMULATED DATA ######################
simu_indiv_test = function(ID, t, kern = kernel_mu, theta, mean, var)
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
              'Output' = rmvnorm(1, rep(mean,length(t)), kern_to_cov(t, kern, theta, sigma = var)) %>% as.vector())
  return(db)
}

#set.seed(42)
M = 20
N = 20
t_i = sample(seq(10, 20, 0.05),N, replace = F) %>% sort()
t = matrix(0, ncol = N, nrow = M)
for(i in 1:M){t[i,] = t_i}

db_train = simu_indiv_test(ID = '1', t[1,], kernel_mu, theta = c(2,1), mean = 45, var = 0.2)
for(i in 2:M)
{
  k = i %% 5
  if(k == 0){db_train = rbind(db_train, simu_indiv_test(ID = as.character(i), t[i,], kernel_mu, theta = c(2,2), mean = 45, var = 0.6))}
  if(k == 1){db_train = rbind(db_train, simu_indiv_test(ID = as.character(i), t[i,], kernel_mu, theta = c(2,1), mean = 45, var = 0.2))}
  if(k == 2){db_train = rbind(db_train, simu_indiv_test(ID = as.character(i), t[i,], kernel_mu, theta = c(3,2), mean = 45, var = 0.3))}
  if(k == 3){db_train = rbind(db_train, simu_indiv_test(ID = as.character(i), t[i,], kernel_mu, theta = c(1,2), mean = 45, var = 0.4))}
  if(k == 4){db_train = rbind(db_train, simu_indiv_test(ID = as.character(i), t[i,], kernel_mu, theta = c(1,1), mean = 45, var = 0.5))}
}
db_obs = simu_indiv_test(ID = (M+1) %>% as.character(), t_i,
                    kernel_mu, theta = c(2,1), mean = 45, var = 0.2)

# ################ INITIALISATION ######################
ini_hp = list('theta_0' = c(1,1), 'theta_i' = c(1, 1, 0.2))
m_0 = 0

##### Testing codes ############

## Testing the full_algo function
# bla = training(db_train, 0, ini_hp, kernel_mu, kernel)
# list_hp_test = bla[c('theta_0','theta_i')]
# blab = full_algo(db_train, db_obs[1:5,], seq(9, 20, 0.02), kernel, common_hp = T, plot = T, prior_mean = m_0,
#                  kern_0 = kernel_mu, list_hp = list_hp_test, mu = NULL, ini_hp = ini_hp, hp_new_i = NULL)
# plot_gp(blab$Prediction, db_obs[3:7,])

# hp_one_gp = list('theta' = c(5,2), 'sigma' = 0.2)
# fu = pred_gp(db_obs[3:7,], seq(10, 20, 0.05), mean_mu = 0 , cov_mu = NULL, kernel,
#              theta = c(5,2), 0.2)
# plot_gp(fu, db_obs[3:7,])

# Testing the training function
# common_hp = F
# bla = training(db_train, 5, ini_hp, kernel_mu, kernel, common_hp)
# fu = bla$param$mean$Output
# names(fu) = paste0('X', bla$param$mean$Timestamp)
# timestamps = seq(10,20, 0.01)
# post_mu = posterior_mu(db_train, timestamps, m_0 = 40, 
#                        kernel_mu, kernel, list('theta_0' = bla$theta_0, 'theta_i' = bla$theta_i))
# 
# hp_pred = list('theta' = bla$theta_i[['1']][1:2], 'sigma' =  bla$theta_i[['1']][[3]])
# if(!common_hp) hp_pred = train_new_gp(db_obs, post_mu$mean, post_mu$cov, ini_hp$theta_i, kernel)
# 
# pred_gp(db_obs[1:3,], timestamps = timestamps, mean_mu = post_mu$mean,
#         cov_mu = post_mu$cov, theta = hp_pred$theta, sigma = hp_pred$sigma) %>%
#   plot_gp(data = rbind(db_obs[1:10,], db_train)) + geom_point(aes(post_mu$mean$Timestamp, post_mu$mean$Output))

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
#
### Creating Heatmaps and GIF
# bla = training(db_train, 0, ini_hp, kernel_mu, kernel)
# post_mu = posterior_mu(db_train, db_obs, seq(10, 20, length.out = 500), 0, kernel_mu, kernel, bla$hp)
# pred_fu = pred_gp(db_obs[1:5,], seq(10,20, length.out = 500), mean_mu = post_mu$mean, cov_mu = post_mu$cov,
#                      kern = kernel, theta = bla$hp$theta_i$`1`[1:2], sigma = bla$hp$theta_i$`1`[[3]])
# plot_heat(pred_fu, data = db_obs, data_train = db_train, mean = post_mu$pred_GP, ygrid = NULL, interactive = F, CI = T)
# fu = pred_gp_animate(db_obs, seq(10,20, length.out = 500), mean_mu = post_mu$mean, cov_mu = post_mu$cov,
#                      kern = kernel, theta = bla$hp$theta_i$`1`[1:2], sigma = bla$hp$theta_i$`1`[[3]])
# # plot_fu = plot_animate(fu, db_obs, db_train, post_mu$pred_GP, F)
# anim_save('GIF_GP.gif', animation = plot_fu, path = NULL)
