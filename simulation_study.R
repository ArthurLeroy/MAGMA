setwd(dir = 'C:/Users/user/CloudStation/Maths/These/Processus Gaussiens/Code R/Algo multitask GP')
source('Algo_multitask_GP.R')

library(GPFDA)

##### COMPETING ALGO IN SIMU ##
train_gpfda = function(db)
{
  t1 = Sys.time()
  db = db %>% dplyr::select(ID, Timestamp, Output)
  M = db$ID %>% unique %>% length
  t_obs = db$Timestamp %>% unique
  N = t_obs %>% length
  
  lx = rep(1, M) %>% as.matrix
  
  fy1 = db %>% spread(key = ID, value = Output) %>% dplyr::select(- Timestamp ) %>% as.matrix %>% t()
  fx1 = db$Timestamp %>% matrix(nrow = N, ncol = M) %>% t()
  
  model = gpfr(response=(fy1), lReg=lx, fReg=NULL, gpReg=list(fx1), fyList=list(nbasis=23,lambda=0.1), fbetaList_l=NULL,
           hyper=NULL, Cov=c('pow.ex', 'linear'), fitting=T, time = t_obs, rPreIdx=T, concurrent=T) 
  t2 = Sys.time()
  model[['Time_train']] = as.numeric(t2 - t1)
  
  return(model)
}

pred_gpfda = function(model, new_db, timestamps)
{
  a1 = model
  tfx = timestamps %>% as.matrix
  time_new = timestamps
  
  pfy = new_db$Output
  pfx = new_db$Timestamp %>% as.matrix
  ptime = new_db$Timestamp
  
  b1<-gpfrpred(a1, TestData= tfx, NewTime=time_new, lReg = 1, fReg=NULL,
               gpReg=list('response'=(pfy),'input'=(pfx),'time'= ptime))
  
  tibble('Timestamp' = b1$predtime, 'Mean' = as.vector(b1$ypred.mean), 'Var' = as.vector(b1$ypred.sd)^2) %>% return()
}

##### SIMULATION FUNCTIONS #####
simu_indiv = function(ID, t, mean, kern, a, b, sigma)
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
              'Output' = rmvnorm(1, mean, kern_to_cov(t, kern, theta = c(a,b), sigma)) %>% as.vector(),
              'a' = a,
              'b' = b,
              'sigma' = sigma)
  return(db)
}

draw = function(int)
{
  return(runif(1,int[1], int[2]) %>% round(2))
}

prior_mean = function(t)
{
  a = draw(c(-2, 2))
  b = draw(c(0,10))
  return(a * t + b)
}

simu_scheme = function(M = 10, N = 10, G = seq(0, 10, 0.05), common_times = T, common_hp = T, kern_0 = kernel_mu, kern_i = kernel,
                       int_mu_a = c(0,5),
                       int_mu_b = c(0,2),
                       int_i_a = c(0,5),
                       int_i_b = c(0,2),
                       int_i_sigma = c(0,1))
{
  m_0 = prior_mean(G)
  mu_a = draw(int_mu_a)
  mu_b = draw(int_mu_b)
  
  db_0 = simu_indiv('0', G, m_0, kern_0, mu_a, mu_b, sigma = 0)
  if(common_times){t_i = sample(G, N, replace = F) %>% sort()}
  if(common_hp)
    {
      i_a = draw(int_i_a)
      i_b = draw(int_i_b)
      i_sigma = draw(int_i_sigma)
    }
  
  floop = function(i)
  { 
    if(!common_times){t_i = sample(G, N, replace = F) %>% sort()}
    if(!common_hp)
    {
      i_a = draw(int_i_a)
      i_b = draw(int_i_b)
      i_sigma = draw(int_i_sigma)
    }
    mean_i = db_0 %>% filter(Timestamp %in% t_i) %>% pull(Output)
    
    simu_indiv(as.character(i), t_i, mean_i, kern_i, i_a, i_b, i_sigma) %>% return()
  }
  db_i = lapply(seq_len(M), floop)
  db = do.call("rbind", db_i) %>% rbind(db_0)
  return(db)
}


##### LOOP DATA FUNCTIONS ####

datasets_multi_M = function(rep, vec_M, N, G, common_times, common_hp, kern_0, kern_i,
                            int_mu_a, int_mu_b, int_i_a, int_i_b, int_i_sigma)
{ ## rep : number of dataset to draw for each value of M
  ## vec_M : vector of values of M for which we want data sets
  ## other inputs : same as simu_scheme function
  ####  
  ## return : tibble of rep x length(M) binded datasets. Columns 'nb_M' and 'ID_dataset' are added to distinguish them
  multi_db = tibble()
  for(i in vec_M)
  {
    for(j in seq_len(rep))
    {
      multi_db = rbind(multi_db, simu_scheme(i, N, G, common_times, common_hp, kern_0, kern_i, int_mu_a, int_mu_b,
                                             int_i_a, int_i_b, int_i_sigma) %>% mutate('nb_M' = i, 'ID_dataset' = j))
    }
  }
  return(multi_db)
}


datasets_multi_N = function(rep, M, N, G, common_times, common_hp, kern_0, kern_i,
                            int_mu_a, int_mu_b, int_i_a, int_i_b, int_i_sigma)
{ ## rep : number of dataset to draw
  ## other inputs : same as simu_scheme function
  ####  
  ## return : tibble of rep x length(M) binded datasets. Columns 'nb_M' and 'ID_dataset' are added to distinguish them
  multi_db = tibble()

    for(j in seq_len(rep))
    {
      multi_db = rbind(multi_db, simu_scheme(M, N, G, common_times, common_hp, kern_0, kern_i, int_mu_a, int_mu_b,
                                             int_i_a, int_i_b, int_i_sigma) %>% mutate('ID_dataset' = j))
    }
  return(multi_db)
}

##### EVALUATION FUNCTIONS #####

MSE = function(error)
{ ## return : the RMSE given the vector of errors
  mean(error^2) %>% return()
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
  
  eval_algo = tibble('MSE' = loss(pred_algo, db_test) %>% MSE(),
                       'Ratio_IC' = ratio_IC(db_test, pred_algo - 1.96 * sd_algo, pred_algo + 1.96 * sd_algo),
                     'Time_train' = db_results$Time_train_algo, 'Time_pred' = db_results$Time_pred_algo)
  eval_one_gp = tibble('MSE' = loss(pred_one_gp, db_test) %>% MSE(),
                     'Ratio_IC' = ratio_IC(db_test, pred_one_gp - 1.96 * sd_one_gp, pred_one_gp + 1.96 * sd_one_gp),
                     'Time_train' = 0, 'Time_pred' = db_results$Time_pred_one_gp)
  eval_gpfda = tibble('MSE' = loss(pred_gpfda, db_test) %>% MSE(),
                    'Ratio_IC' = ratio_IC(db_test, pred_gpfda - 1.96 * sd_gpfda, pred_gpfda + 1.96 * sd_gpfda),
                    'Time_train' = db_results$Time_train_gpfda, 'Time_pred' = db_results$Time_pred_gpfda)

  rbind(eval_algo, eval_one_gp, eval_gpfda) %>% mutate(Method = c('Algo', 'One GP', 'GPFDA')) %>% return()
}

##### FULL SIMULATION FUNCTION #####

loop_training = function(db_loop, prior_mean, ini_hp, kern_0, kern_i, diff_M, common_times, common_hp)
{  
  ## Loop over the different datasets
  floop = function(i)
  {
    print(paste0('Dataset n°', i))
    ## Select the i-th dataset and remove mean process and testing individual (ID = 0 and 1)
    db_train = db_M %>% filter(ID_dataset == i) %>% filter(!(ID %in% c(0,1))) %>%
                           dplyr::select('ID', 'Timestamp', 'Output')
    
    if(common_times){model_gpfda = train_gpfda(db_train)}
    list_hp = training(db_train, prior_mean, ini_hp, kern_0, kern_i, common_hp)[c('hp', 'Time_train')]
    list('gpfda' = model_gpfda, 'algo' = list_hp) %>% return()
  }
  
  if(diff_M)
  { 
    list_train = list()
    ## Removing data with only 1 individual (= the testing individual) or 2 indiv (gpfda doesn't run)
    list_value_M = unique(db_loop$nb_M) %>% subset(. %notin% c(1,2))
    ## Loop over the different values of M (optional)
    for(j in list_value_M)
    {
      db_M = db_loop %>% filter(nb_M == j)
      list_train[[paste0('M=', j)]] = unique(db_M$ID_dataset) %>% as.character() %>%
                                      sapply(floop, simplify = FALSE, USE.NAMES = TRUE)
    }
  }
  else
  {
    db_M = db_loop
    list_train = unique(db_loop$ID_dataset) %>% as.character() %>% sapply(floop, simplify = FALSE, USE.NAMES = TRUE)
  }
  
  list_train %>% c(list('prior_mean' = prior_mean, 'ini_hp' = ini_hp, 'kern_0' = kern_0,'kern_i' =  kern_i,
                        'diff_M' = diff_M, 'common_times' = common_times, 'common_hp' = common_hp)) %>% 
  return()
}

loop_pred = function(db_loop, train_loop, nb_obs, nb_test)
{ ## Get the settings used for training
  prior_mean = train_loop$prior_mean
  ini_hp = train_loop$ini_hp
  kern_0 = train_loop$kern_0
  kern_i = train_loop$kern_i
  diff_M = train_loop$diff_M
  common_times = train_loop$common_times
  common_hp = train_loop$common_hp
  
  floop = function(i)
  { 
    ## Get the trained model for GPFDA and our algo
    model_algo = train_M[[i]]$algo
    list_hp = model_algo$hp
    model_gpfda = train_M[[i]]$gpfda
 
    ## Get the corresponding database
    db_i = db_M %>% filter(ID_dataset == i)
    db_train_i = db_i %>% filter(ID %notin% c(0,1))
    ## Select the 'nb_obs' first observations of the testing individual to predict with
    db_obs_i = db_i %>% filter(ID == 1) %>% top_n(- nb_obs, Timestamp)
    ## Select the 'nb_test last observations of the testing individual to evaluate predictions on 
    db_pred_i = db_i %>% filter(ID == 1) %>% top_n(nb_test, Timestamp)
    ## Get timestamps to predict on
    t_i_pred = db_pred_i %>% pull(Timestamp)
    
    t1 = Sys.time()
    ## Prediction for our algo (train new indiv + pred)
    res_algo = full_algo(db_train_i, db_obs_i, t_i_pred, kern_i, common_hp, plot = F, prior_mean, kern_0,
                         list_hp, mu = NULL, ini_hp, hp_new_i = NULL)$Prediction
    t2 = Sys.time()
    ## Train new indiv for one GP model
    hp_one_gp_hp = train_new_gp(db_obs_i, rep(prior_mean, nrow(db_obs_i)), cov_mu = 0, ini_hp$theta_i, kern_i)
    ## Prediction for one GP model
    res_one_gp = pred_gp(db_obs_i, t_i_pred, prior_mean, cov_mu = NULL, kern_i, hp_one_gp_hp$theta, hp_one_gp_hp$sigma)
    t3 = Sys.time()
    ## Prediction for GPFDA (unable for uncommon timestamps)
    if(common_times){res_gpfda = pred_gpfda(model_gpfda, db_obs_i, t_i_pred)}
    else{res_gpfda =  tibble('Timestamp' = t_i_pred, 'Mean' = NA ,   'Var' = NA)}
    t4 = Sys.time()
    
    ### Get MSE, RATIO IC95 and computing times on testing points for all methods 
    list('algo' = res_algo, 'Time_train_algo' = model_algo$Time_train, 'Time_pred_algo' = as.numeric(t2 - t1) ,
         'one_gp' = res_one_gp, 'Time_pred_one_gp' = as.numeric(t3 - t2)  ,
         'gpfda' = res_gpfda, 'Time_train_gpfda' = model_gpfda$Time_train, 'Time_pred_gpfda' = as.numeric(t4-t3)) %>%
    eval_methods(db_pred_i %>% pull(Output)) %>%
    return()
  }
  
  if(diff_M)
  {
    list_eval = tibble()
    ## Removing data with only 1 individual (= the testing individual) or 2 indiv (gpfda doesn't run)
    list_value_M = unique(db_loop$nb_M) %>% subset(. %notin% c(1,2))
    ## Loop over the different values of M (optional)
    for(j in list_value_M)
    { 
      db_M = db_loop %>% filter(nb_M == j)
      train_M = train_loop[[paste0('M=', j)]]
      eval_M = unique(db_M$ID_dataset) %>% as.character() %>% sapply(floop, simplify = FALSE, USE.NAMES = TRUE)
      list_eval = list_eval %>% rbind(do.call('rbind', eval_M) %>% mutate('M' = j))
    }
    return(list_eval)
  }
  else
  {
    db_M = db_loop
    train_M = train_loop
    list_eval = db_M$ID_dataset %>% unique() %>% lapply(floop)
    table_eval = do.call('rbind', list_eval) 
    # %>% group_by(Method) %>%
    #              summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE) %>% return()
    # 
  }
}

simu_var_N = function(db, nb_obs_max, nb_test, prior_mean, ini_hp, kern_0, kern_i, diff_M, common_times, common_hp,
                      plot = T)
{
  loop_train = loop_training(db, prior_mean, ini_hp, kern_0, kern_i, diff_M = F, common_times, common_hp)
  
  floop = function(i)
  {
    loop_pred = loop_pred(db, loop_train, nb_obs = i, nb_test, prior_mean, ini_hp, 
                          kern_0, kern_i, diff_M = F, common_times, common_hp) %>% mutate('N' = i) %>% 
      return()
  }
  res_n = lapply(1:nb_obs_max, floop) %>% do.call('rbind', .)
  
  if(plot)
  {
    ggplot(res_n) + geom_boxplot(aes(x = as.factor(N), y = MSE_Mean, fill = Method))
  }
  return(res_n)
}

simu_var_M = function(db,  prior_mean, ini_hp, kern_0, kern_i, diff_M = T, common_times, common_hp, graph = T)
{
  
}


##### INITIALISATION #####

# M = 10
# N = 10
# G = seq(0,10, 0.05)
# int_mu_a = c(0,10)
# int_mu_b = c(0,10)
# int_i_a = c(0,10)
# int_i_b = c(0,10)
# int_i_sigma = c(0,10)


##### SIMULATION DATASETS ####

# set.seed(42)
# for(i in c(T, F))
# {
#   for(j in c(T, F))
#   {
#     datasets_multi_N(rep = 100, M = 21, N = 30, G = seq(0, 10, 0.05), common_times = i,
#                      common_hp = j, kern_0 = kernel_mu, kern_i = kernel, int_mu_a = c(0,5), int_mu_b = c(0,2),
#                      int_i_a = c(0,5), int_i_b = c(0,2), int_i_sigma = c(0,2)) %>%
#     write_csv2(paste0('Simulations/Data/db_rep_100_M_20_N_30_time_', i, '_hp_', j,'.csv'))
#   }
# }

set.seed(42)
for(i in c(T, F))
{
  for(j in c(T, F))
  {
    datasets_multi_M(rep = 100, vec_M = c(21, 51, 101, 151, 201), N = 30, G = seq(0, 10, 0.05), common_times = i,
                     common_hp = j, kern_0 = kernel_mu, kern_i = kernel, int_mu_a = c(0,5), int_mu_b = c(0,2),
                     int_i_a = c(0,5), int_i_b = c(0,2), int_i_sigma = c(0,2)) %>%
      write_csv2(paste0('Simulations/Data/db_rep_100_M_20to200_N_30_time_', i, '_hp_', j,'.csv'))
  }
}

##### TABLE OF RESULTS ####

test_res_n = simu_var_N(bla2, nb_obs_max = 20, nb_test = 10, prior_mean = 0,
                        ini_hp = list('theta_0' = c(1,1), 'theta_i' = c(1, 1, 0.2)), 
                        kern_0 = kernel_mu, kern_i = kernel, diff_M = F, common_times = T, common_hp = T)


##### PLOT OF RESULTS #### 
# ggplot(test_res_n) + geom_boxplot(aes(x = as.factor(N), y = MSE, fill = Method), outlier.shape = NA)
# ggplot(test_loop_pred) + geom_boxplot(aes(x = as.factor(M), y = MSE, fill = Method))


##### TESTS SIMU ####
# bla_db = datasets_multi_N(rep = 10, M = 21, N = 30, G = seq(0, 10, 0.05), common_times = T,
#                           common_hp = F, kern_0 = kernel_mu, kern_i = kernel, int_mu_a = c(0,5), int_mu_b = c(0,2),
#                           int_i_a = c(0,5), int_i_b = c(0,2), int_i_sigma = c(0,1))
#
# bla_db2 = datasets_multi_M(rep = 10, vec_M = c(3, 11, 21), N = 30, G = seq(0, 10, 0.05), common_times = T,
#                           common_hp = F, kern_0 = kernel_mu, kern_i = kernel, int_mu_a = c(0,5), int_mu_b = c(0,2),
#                           int_i_a = c(0,5), int_i_b = c(0,2), int_i_sigma = c(0,1))
# 
# test_loop_train = loop_training(bla_db, prior_mean = 0, ini_hp = list('theta_0' = c(1,1), 'theta_i' = c(1, 1, 0.2)), 
#                                 kern_0 = kernel_mu, kern_i = kernel, diff_M = T, common_times = T, common_hp = T)
# 
# test_loop_pred = loop_pred(bla_db, test_loop_train, nb_obs =20, nb_test = 10)
# 

##### OLD SIMULATION STUDY ##### 
# db_train = simu_scheme(M = 10, N = 10, G = seq(0, 10, 0.05))
# plot_db(db_train)
# 
# eval = simu_study(M = 10, N = 10, G = seq(0, 10, 0.05), prior_mean = 0, kern_0 = kernel_mu, kern_i = kernel,
#                   ini_hp, common_hp = T, common_times = T,
#                   ratio_train = 0.6,
#                   int_mu_a = c(0,5),
#                   int_mu_b = c(0,2),
#                   int_i_a = c(0,5),
#                   int_i_b = c(0,2),
#                   int_i_sigma = c(0,1),
#                   int_test = c(2,9))


###### TEST GPFDA ####
# M = 20
# N = 11
# t = matrix(0, ncol = N, nrow = M)
# for(i in 1:M){t[i,] = seq(0,10, length.out = N)}
# 
# db_train = simu_indiv(ID = '1', t[1,], kernel_mu, theta = c(2,1), mean = 45, var = 0.2)
# for(i in 2:M)
# {
#   k = i %% 5
#   if(k == 0){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(2,2), mean = 45, var = 0.6))}
#   if(k == 1){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(2,1), mean = 45, var = 0.2))}
#   if(k == 2){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(3,2), mean = 45, var = 0.3))}
#   if(k == 3){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(1,2), mean = 45, var = 0.4))}
#   if(k == 4){db_train = rbind(db_train, simu_indiv(ID = as.character(i), t[i,], kernel_mu, theta = c(1,1), mean = 45, var = 0.5))}
# }
# db_obs = simu_indiv(ID = (M+1) %>% as.character(), seq(0,10, length.out = N),
#                     kernel_mu, theta = c(2,1), mean = 40, var = 0.2)
# 
# 
# plot(b1,type='prediction')
# 
# plot(-1000,col=0,xlim=range(b1$time),ylim=range(b1$ypred),xlab='time',ylab='prediction',
#      main='Prediction by GPFR: type I')
# 
# lines(b1$predtime,b1$ypred[,1])
# lines(b1$predtime,b1$ypred[,2],lty=2,col=2)
# lines(b1$predtime,b1$ypred[,3],lty=2,col=2)
# points(xt,yt)
# 

##### OLD FUNCTIONS ####
# split_train = function(db, ratio_train)
# { ## db : the database of all observed data
#   ## ratio : number between 0 and 1 indicating the proportion of individuals in the training set
#   #
#   ## return : orginal db with the repartition of individuals between the training and testing set
#   n_indiv = db$ID%>% unique()
#   n_index = (ratio_train * length(n_indiv)) %>% round()
#   index_train = sample(n_indiv, n_index, replace = F) 
#   
#   db %>% mutate(Training = ifelse(ID %in% index_train, 1,0)) %>% return()
# }
# 
# split_times = function(db, int_test)
# {
#   
#   db = db %>% rowid_to_column("Row")
#   row_obs = db %>% group_by(ID) %>% sample_n(sample(int_test[1]:int_test[2],1)) %>% pull(Row)
#   
#   db %>% mutate(Observed = ifelse(db$Row %in% row_obs, 1, 0)) %>% dplyr::select(- Row) %>% return()
# }

# simu_study = function(M, N, G, prior_mean, kern_0, kern_i, ini_hp, common_hp, common_times, ratio_train,
#                       int_mu_a, int_mu_b, int_i_a, int_i_b, int_i_sigma, int_test)
# {
#   db_full = simu_scheme(M, N, G, common_times,common_hp, kern_0, kern_i, int_mu_a, int_mu_b, int_i_a, int_i_b, int_i_sigma)
#   
#   mean_process = db_full %>% filter(ID == 0)
#   db = db_full %>% filter(ID != 0) %>% split_train(ratio_train)
#   
#   db_train = db %>% filter(Training == 1)
#   db_test = db %>% filter(Training == 0) %>% split_times(int_test)
#   
#   if(common_times){model_gpfda = train_gpfda(db_train)}
#   list_hp = training(db_train, prior_mean, ini_hp, kern_0, kern_i, common_hp)$hp
#   
#   floop = function(i)
#   {
#     db_obs_i = db_test %>% filter(ID == i) %>% filter(Observed == 1)
#     db_pred_i = db_test %>% filter(ID == i) %>% filter(Observed == 0)
#     t_i_pred = db_pred_i %>% pull(Timestamp)
#     
#     res_algo = full_algo(db_train, db_obs_i, t_i_pred, kern_i, common_hp, plot = F, prior_mean, kern_0,
#                          list_hp, mu = NULL, ini_hp, hp_new_i = NULL)$Prediction
#     
#     hp_one_gp_hp = train_new_gp(db_obs_i, rep(prior_mean, nrow(db_obs_i)), cov_mu = 0, ini_hp$theta_i, kern_i)
#     
#     res_one_gp = pred_gp(db_obs_i, t_i_pred, prior_mean, cov_mu = NULL, kern_i, hp_one_gp_hp$theta, hp_one_gp_hp$sigma)
#     
#     if(common_times){res_gpfda = pred_gpfda(model_gpfda, db_obs_i, t_i_pred)}
#     else{res_gpfda =  tibble('Timestamp' = t_i_pred, 'Mean' = NA , 'Var' = NA)}
#     
#     list('algo' = res_algo, 'one_gp' = res_one_gp, 'gpfda' = res_gpfda) %>%
#       eval_methods(db_pred_i %>% pull(Output)) %>% return()
#   }
#   list_eval = db_test$ID %>% unique() %>% lapply(floop)
#   print(list_eval)
#   table_eval = do.call('rbind', list_eval) %>% group_by(Method) %>%
#     summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE)
#   
#   return(table_eval)
# }
