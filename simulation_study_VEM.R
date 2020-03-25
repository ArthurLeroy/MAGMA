setwd(dir = 'C:/Users/user/CloudStation/Maths/These/Processus Gaussiens/Code R/Algo multitask GP')
source('Algo_multitask_GP.R')

##### SIMULATION FUNCTIONS #####
simu_indiv = function(ID, t, mean, kern, a, b, sigma, clust)
{ # ID : identification of the individual
  # t : timestamps on which we observe the GP
  # kern : kernel associated to the covariance function of the GP
  # theta : list of hp of the kernel
  # mean : mean parameter of the GP
  # var : variance parameter of the error
  ##
  # return : a simulated individual
  db = tibble('ID' = ID,
              'Timestamp' = t, 
              'Input' = paste0('X', t),
              'Output' = rmvnorm(1, mean, kern_to_cov(t, kern, theta = c(a,b), sigma)) %>% as.vector(),
              'a' = a,
              'b' = b,
              'sigma' = sigma, 
              'Cluster' = clust)
  return(db)
}

draw = function(int)
{
  return(runif(1,int[1], int[2]) %>% round(2))
}

prior_mean = function(t)
{
  if(rbernoulli(1, 0.5)){b = draw(c(7,10))} else {b = draw(c(0,3))}
  a = draw(c(-1, 1))
  return(a * t + b)
}

simu_scheme = function(M = 10, N = 10, K = 2, G = seq(0, 10, 0.05), kern_0 = kernel_mu, kern_i = kernel,
                       int_mu_a = c(0,5),
                       int_mu_b = c(0,2),
                       int_i_a = c(0,1),
                       int_i_b = c(0,1),
                       int_i_sigma = c(0,0.5))
{
  Z_i = sample(K, M, replace = T)
  
  db_k = tibble()
  for(k in 1:K)
  { 
    m_k = prior_mean(G)
    mu_a = draw(int_mu_a)
    mu_b = draw(int_mu_b)
    
    db_k = rbind(db_k , simu_indiv(paste0('K',k), G, m_k, kern_0, mu_a, mu_b, 0, paste0('K', k)))
  }
  
  floop = function(i)
  { 
    t_i = sample(G, N, replace = F) %>% sort()
    i_a = draw(int_i_a)
    i_b = draw(int_i_b)
    i_sigma = draw(int_i_sigma)
    
    mean_i = db_k %>% filter(ID == paste0('K', Z_i[[i]])) %>% filter(Timestamp %in% t_i) %>% pull(Output)
    
    simu_indiv(as.character(i), t_i, mean_i, kern_i, i_a, i_b, i_sigma, paste0('K', Z_i[[i]])) %>% return()
  }
  db_i = lapply(seq_len(M), floop)
  db = do.call("rbind", db_i) %>% rbind(db_k)
  return(db)
}

plot_db = function(db = db_train)
{
  ggplot(db) + geom_smooth(aes(Timestamp, Output, color = ID, shape = Cluster)) +
    geom_point(aes(Timestamp, Output, color = ID, shape = Cluster))
}

split_train = function(db, ratio_train)
{ ## db : the database of all observed data
  ## ratio : number between 0 and 1 indicating the proportion of individuals in the training set
  ####
  ## return : orginal db with the repartition of individuals between the training and testing set
  n_indiv = db$ID%>% unique()hp
  n_index = (ratio_train * length(n_indiv)) %>% round()
  index_train = sample(n_indiv, n_index, replace = F) 
  
  db %>% mutate(Training = ifelse(ID %in% index_train, 1,0)) %>% return()
}

split_times = function(db, int_test)
{
  
  db = db %>% rowid_to_column("Row")
  row_obs = db %>% group_by(ID) %>% sample_n(sample(int_test[1]:int_test[2],1)) %>% pull(Row)
  
  db %>% mutate(Observed = ifelse(db$Row %in% row_obs, 1, 0)) %>% dplyr::select(- Row) %>% return()
}

##### EVALUATION FUNCTIONS #####

RMSE = function(error)
{ ## return : the RMSE given the vector of errors
  sqrt(mean(error^2)) %>% return()
}

loss = function(x, y)
{ ## return : loss function between x and y
  abs(x - y) %>% return()
}

ratio_IC = function(obs, IC_inf, IC_sup)
{ ## obs : the true observed values
  ## IC_inf : inferior boundary of the predicted IC_0.95
  ## IC_sup : superior boundary of the predicted IC_0.95
  ####
  ## return : the ratio of observed values lying within the predicted IC_0.95
  nb_between = ((IC_inf < obs) & (obs < IC_sup)) %>% sum()
  nb_between/ length(obs) * 100 %>% return()
}

eval_methods = function(db_results, db_test)
{ ## db_results : list of results of the different methods. Format : list('algo','one_gp','gpfda')
  ## db_test : vector of observed values on which we test the predictions
  ####
  ## return : Table of results of the evaluation of the different methods through RMSE and ratio_IC
  pred_algo = db_results$algo$Mean 
  sd_algo = db_results$algo$Var %>% sqrt()
  
  pred_one_gp = db_results$one_gp$Mean
  sd_one_gp = db_results$one_gp$Var %>% sqrt()
  
  pred_gpfda = db_results$gpfda$Mean
  sd_gpfda = db_results$gpfda$Var %>% sqrt()
  
  eval_algo = tibble('RMSE' = loss(pred_algo, db_test) %>% RMSE(),
                     'Ratio_IC' = ratio_IC(db_test, pred_algo - 1.96 * sd_algo, pred_algo + 1.96 * sd_algo))
  eval_one_gp = tibble('RMSE' = loss(pred_one_gp, db_test) %>% RMSE(),
                       'Ratio_IC' = ratio_IC(db_test, pred_one_gp - 1.96 * sd_one_gp, pred_one_gp + 1.96 * sd_one_gp))
  eval_gpfda = tibble('RMSE' = loss(pred_gpfda, db_test) %>% RMSE(),
                      'Ratio_IC' = ratio_IC(db_test, pred_gpfda - 1.96 * sd_gpfda, pred_gpfda + 1.96 * sd_gpfda))
  
  rbind(eval_algo, eval_one_gp, eval_gpfda) %>% mutate(Method = c('Algo', 'One GP', 'GPFDA')) %>% return()
}

##### FULL SIMULATION FUNCTION #####

simu_study = function(M, N, G, prior_mean, kern_0, kern_i, ini_hp, ratio_train, int_mu_a, int_mu_b, int_i_a,
                      int_i_b, int_i_sigma, int_test)
{
  db_full = simu_scheme(M, N, G, kern_0, kern_i, int_mu_a, int_mu_b, int_i_a, int_i_b, int_i_sigma)
  
  mean_process = db_full %>% filter(ID == 0)
  db = db_full %>% filter(ID != 0) %>% split_train(ratio_train)
  
  db_train = db %>% filter(Training == 1)
  db_test = db %>% filter(Training == 0) %>% split_times(int_test)
  
  list_hp = training(db_train, prior_mean, ini_hp, kern_0, kern_i)
  
  floop = function(i)
  {
    db_obs_i = db_test %>% filter(ID == i) %>% filter(Observed == 1)
    db_pred_i = db_test %>% filter(ID == i) %>% filter(Observed == 0) 
    t_i_pred = db_pred_i %>% pull(Timestamp)
    res_algo = full_algo(db_train, db_obs_i, t_i_pred, kern_i, plot = F, prior_mean, kern_0,
                         list_hp, mu = NULL, ini_hp, hp_new_i = NULL)$Prediction
    
    hp_one_gp_i = train_new_gp(db_obs_i, rep(prior_mean, nrow(db_obs_i)), cov_mu = 0, ini_hp$theta_i, kern_i)
    res_one_gp = pred_gp(db_obs_i, t_i_pred, prior_mean, cov_mu = NULL, kern_i, hp_one_gp_i$theta, hp_one_gp_i$sigma)
    
    res_gpfda = full_algo(db_train, db_obs_i, t_i_pred, kern_i, plot = F, prior_mean, kern_0,
                          list_hp, mu = NULL, ini_hp, hp_new_i = NULL)$Prediction
    
    list('algo' = res_algo, 'one_gp' = res_one_gp, 'gpfda' = res_gpfda) %>% 
      eval_methods(db_pred_i %>% pull(Output)) %>% return() 
  }
  list_eval = db_test$ID %>% unique() %>% lapply(floop)
  table_eval = do.call('rbind', list_eval) %>% group_by(Method) %>% summarise_all(list('Mean' = mean, 'SD' = sd))
  
  return(table_eval)
}

##### SIMULATION STUDY ##### 
db_train = simu_scheme(M = 10, N = 10, G = seq(0, 10, 0.05))
plot_db(db_train)

eval = simu_study(M = 10, N = 10, G = seq(0, 10, 0.1), prior_mean = 0, kern_0 = kernel_mu, kern_i = kernel,
                  ini_hp,
                  ratio_train = 0.6,
                  int_mu_a = c(0,5),
                  int_mu_b = c(0,2),
                  int_i_a = c(0,5),
                  int_i_b = c(0,2),
                  int_i_sigma = c(0,1), 
                  int_test = c(2,9))


##### INITIALISATION #####

# M = 10
# N = 10
# G = seq(0,10, 0.05)
# int_mu_a = c(0,10)
# int_mu_b = c(0,10)
# int_i_a = c(0,10)
# int_i_b = c(0,10)
# int_i_sigma = c(0,10)
# int_test = c(2, (N-1))
#ratio_train_set = 0.3

