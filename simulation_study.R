#setwd(dir = 'C:/Users/user/Google Drive/Travail/GitHub/gpclust/')
source('Algo_multitask_GP.R')

library(GPFDA)
library(parallel)

##### COMPETING ALGO IN SIMU ####
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
  model[['Time_train']] =  difftime(t2, t1, units = "secs")
  
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

mean_gpfda =  function(model, timestamps)
{
  a1 = model
  tfx = timestamps %>% as.matrix
  time_new = timestamps
  
  b1<-gpfrpred(a1, TestData= tfx, NewTime=time_new, lReg = 1, fReg=NULL,
               gpReg= NULL)
  
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
      multi_db = rbind(multi_db, simu_scheme(i, N, G, common_times, common_hp, kern_0, kern_i, 
                                             int_mu_a, int_mu_b,int_i_a, int_i_b, int_i_sigma) %>%
                                 mutate('nb_M' = i, 'ID_dataset' = as.character(j)))
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
                                             int_i_a, int_i_b, int_i_sigma) %>% mutate('ID_dataset' = as.character(j)))
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

##### TRAINING AND PRED FUNCTIONS #####

loop_training = function(db_loop, prior_mean, ini_hp, kern_0, kern_i, diff_M, common_times, common_hp)
{  
  ## Loop over the different datasets
  floop = function(i, M = NULL)
  {
    
    if(diff_M){db_M = db_loop %>% filter(nb_M == M)}else{db_M = db_loop}
    if(diff_M){print(paste0('Nb M = ', db_M$nb_M[[1]]))}
    print(paste0('Dataset n°', i))
    
    ## Select the i-th dataset and remove mean process and testing individual (ID = 0 and 1)
    db_train = db_M %>% filter(ID_dataset == i) %>% filter(!(ID %in% c(0,1))) %>%
                           dplyr::select('ID', 'Timestamp', 'Output')
    
    if(common_times){model_gpfda = train_gpfda(db_train)}else{model_gpfda = list('Time_train' = 0)}
    list_hp = training(db_train, prior_mean, ini_hp, kern_0, kern_i, common_hp)[c('hp', 'Time_train')]
    list('gpfda' = model_gpfda, 'algo' = list_hp) %>% return()
  }
  
  if(diff_M)
  { 
    ## Removing data with only 1 individual (= the testing individual) or 2 indiv (gpfda doesn't run)
    list_value_M = unique(db_loop$nb_M) %>% subset(. %notin% c(1,2))
    ## Parallel computing through mclapply not available in windows
    nb_core = ifelse(.Platform$OS.type == 'windows', 1, 1)
    
    ## Parallel computing for different values of M
    list_train <- mclapply(list_value_M, function(j) {
      unique(db_loop$ID_dataset) %>% as.character() %>% 
      sapply(floop, M = j,  simplify = FALSE, USE.NAMES = TRUE) %>% 
      return()
    }, mc.cores= nb_core)
    names(list_train) = paste0('M=', list_value_M)
  }
  else
  {
    list_train = unique(db_loop$ID_dataset) %>% as.character() %>% sapply(floop, simplify = FALSE, USE.NAMES = TRUE)
  }
  
  list_train %>% c(list('prior_mean' = prior_mean, 'ini_hp' = ini_hp, 'kern_0' = kern_0,'kern_i' =  kern_i,
                        'diff_M' = diff_M, 'common_times' = common_times, 'common_hp' = common_hp)) %>% 
  return()
}

loop_pred = function(db_loop, train_loop, nb_obs, nb_test, diff_M)
{
  db_loop$ID = db_loop$ID %>% as.character
  db_loop$ID = db_loop$ID %>% as.character
  ## Get the settings used for training
  prior_mean = train_loop$prior_mean
  ini_hp = train_loop$ini_hp
  kern_0 = train_loop$kern_0
  kern_i = train_loop$kern_i
  diff_M = train_loop$diff_M
  common_times = train_loop$common_times
  common_hp = train_loop$common_hp

  floop = function(i, M = NULL)
  {
    if(diff_M)
    {
      print(paste('M =',M, ', i =',i))
      db_M = db_loop %>% filter(nb_M == M)
      train_M = train_loop[[paste0('M=', M)]]
    }
    else
    {
      print(paste('i =' ,i))
      db_M = db_loop
      train_M = train_loop
    }
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
    hp_one_gp = train_new_gp(db_obs_i, rep(prior_mean, nrow(db_obs_i)), cov_mu = 0, ini_hp$theta_i, kern_i)
    ## Prediction for one GP model
    res_one_gp = pred_gp(db_obs_i, t_i_pred, prior_mean, cov_mu = NULL, kern_i, hp_one_gp$theta, hp_one_gp$sigma)
    t3 = Sys.time()
    ## Prediction for GPFDA (unable for uncommon timestamps)
    if(common_times){res_gpfda = pred_gpfda(model_gpfda, db_obs_i, t_i_pred)}
    else{res_gpfda =  tibble('Timestamp' = t_i_pred, 'Mean' = NA ,   'Var' = NA)}
    t4 = Sys.time()

    ### Get MSE, RATIO IC95 and computing times on testing points for all methods 
    list('algo' = res_algo, 'Time_train_algo' = model_algo$Time_train, 'Time_pred_algo' = difftime(t2, t1, units = "secs"),
         'one_gp' = res_one_gp, 'Time_pred_one_gp' =  difftime(t3, t2, units = "secs"),
         'gpfda' = res_gpfda, 'Time_train_gpfda' = model_gpfda$Time_train, 'Time_pred_gpfda' =  difftime(t4, t3, units = "secs")) %>%
    eval_methods(db_pred_i %>% pull(Output)) %>%
    return()
  }
  
  if(diff_M)
  {
    ## Removing data with only 1 individual (= the testing individual) or 2 indiv (gpfda doesn't run)
    list_value_M = unique(db_loop$nb_M) %>% subset(. %notin% c(1,2))
    ## Parallel computing through mclapply not available in windows
    nb_core = ifelse(.Platform$OS.type == 'windows', 1, 1)
    
    ## Parallel computing for different values of M
    list_M <- mclapply(list_value_M, function(j) {
        list_eval = unique(db_loop$ID_dataset) %>% as.character() %>% 
                    sapply(floop, M = j,  simplify = FALSE, USE.NAMES = TRUE)
      
        rbind(do.call('rbind', list_eval) %>% 
        mutate(Time_train = as.numeric(Time_train), Time_pred = as.numeric(Time_pred), 'M' = j)) %>% 
        return()
    }, mc.cores= nb_core) 
    
    do.call('rbind', list_M) %>%
    return()
  }
  else
  {
    list_eval = db_loop$ID_dataset %>% unique() %>% lapply(floop)
    do.call('rbind', list_eval) %>% 
    mutate(Time_train = as.numeric(Time_train), Time_pred = as.numeric(Time_pred)) %>% 
    return()
  }
}

simu_var_N = function(db, train_loop, nb_obs_max, nb_test)
{
  ## Set IDs as characters in the database
  db$ID = db$ID %>% as.character 
  db$ID_dataset = db$ID_dataset %>% as.character
  
  floop = function(i)
  {
    print(paste('N =', i))
    loop_pred = loop_pred(db, train_loop, nb_obs = i, nb_test) %>% mutate('N' = i) %>% 
    return()
  }
  res_n = lapply(1:nb_obs_max, floop) %>% do.call('rbind', .)
  return(res_n)
}

eval_mu = function(db, train_loop, M = NULL)
{ ## Set IDs as characters in the database
  db$ID = db$ID %>% as.character 
  db$ID_dataset = db$ID_dataset %>% as.character
  ## Set the training parameters
  prior_mean = train_loop$prior_mean
  ini_hp = train_loop$ini_hp
  kern_0 = train_loop$kern_0
  kern_i = train_loop$kern_i
  diff_M = train_loop$diff_M
  common_times = train_loop$common_times
  common_hp = train_loop$common_hp
  
  if(is.null(M))
  {
    train_M = train_loop
    db_M = db
  }
  else
  {
    train_M = train_loop[[paste0('M=', M)]]
    db_M = db %>% filter(nb_M == M)
  }

  floop = function(i)
  { #browser()
    print(paste('M =', M, ', i =', i))
    ## Get the trained model for GPFDA and our algo
    model_algo = train_M[[i]]$algo
    list_hp = model_algo$hp
    model_gpfda = train_M[[i]]$gpfda
    
    ## Get the corresponding database
    db_i = db_M %>% filter(ID_dataset == i)
    db_train_i = db_i %>% filter(ID %notin% c(0,1))
    db_pred_i = db_i %>% filter(ID == 1)
    ## Get timestamps to predict on
    t_pred = db_train_i$Timestamp %>% unique()
    ## Select the true value of the mean process to evaluate predictions on 
    db_pred_mu = db_i %>% filter(ID == 0) %>% filter(Timestamp %in% t_pred)
    
    t1 = Sys.time()
    ## Estimation of the posterior mean process p(mu|data) for our algo
    res_algo = posterior_mu(db_train_i, db_pred_i, t_pred, prior_mean, kern_0, kern_i, list_hp)
    t2 = Sys.time()
    ## No estimation of the mean process with classic GP
    res_one_gp =  tibble('Timestamp' = t_pred, 'Mean' = NA ,   'Var' = NA)
    t3 = Sys.time()
    ## Estimation of the mean process for GPFDA (deterministic function)
    if(common_times){res_gpfda = mean_gpfda(model_gpfda, t_pred)}
    else{res_gpfda =  tibble('Timestamp' = t_pred, 'Mean' = NA ,   'Var' = NA)}
    t4 = Sys.time()

    ### Get MSE, RATIO IC95 and computing times on testing points for all methods 
    list('algo' = res_algo$pred_GP, 'Time_train_algo' = model_algo$Time_train, 'Time_pred_algo' = difftime(t2, t1, units = "secs"),
         'one_gp' = res_one_gp, 'Time_pred_one_gp' =  difftime(t3, t2, units = "secs"),
         'gpfda' = res_gpfda, 'Time_train_gpfda' = model_gpfda$Time_train, 'Time_pred_gpfda' =  difftime(t4, t3, units = "secs")) %>%
      eval_methods(db_pred_mu %>% pull(Output)) %>%
      return()
  }
  
  list_eval = db_M$ID_dataset %>% unique() %>% lapply(floop)
  do.call('rbind', list_eval) %>% 
  mutate(Time_train = as.numeric(Time_train), Time_pred = as.numeric(Time_pred)) %>% 
  return()
}

eval_mu_M = function(db, train_loop)
{
  floop = function(i)
  {
    eval_mu(db, train_loop, M = i) %>% mutate(M = i) %>% return()
  }
  ## Removing data with only 1 individual (= the testing individual) or 2 indiv (gpfda doesn't run)
  M_values = db$nb_M %>% unique %>% subset(. %notin% c(1,2))
  eval_M = M_values %>% lapply(floop)
  do.call('rbind', eval_M) %>% return()
}

##### INITIALISATION #####

# M = 21
# N = 30
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

# set.seed(42)
# for(i in c(T, F))
# {
#   for(j in c(T, F))
#   {
#     datasets_multi_M(rep = 100, vec_M = c(21, 51, 101, 151, 201), N = 30, G = seq(0, 10, 0.05), common_times = i,
#                      common_hp = j, kern_0 = kernel_mu, kern_i = kernel, int_mu_a = c(0,5), int_mu_b = c(0,2),
#                      int_i_a = c(0,5), int_i_b = c(0,2), int_i_sigma = c(0,2)) %>%
#       write_csv2(paste0('Simulations/Data/db_rep_100_M_20to200_N_30_time_', i, '_hp_', j,'.csv'))
#   }
# }


##### TABLES OF DATA ####
tableGPFDA = read_csv2("Simulations/Data/db_GPFDA_M_20_N_30_time_TRUE_hp_TRUE.csv")
tableGPFDA$ID = as.character(tableGPFDA$ID)
tableGPFDA$ID_dataset = as.character(tableGPFDA$ID_dataset)

## Load the data and ensure IDs are filled as characters
tableTT = read_csv2("Simulations/Data/db_rep_100_M_20_N_30_time_TRUE_hp_TRUE.csv")
tableTT$ID = as.character(tableTT$ID)
tableTT$ID_dataset = as.character(tableTT$ID_dataset)
tableTF = read_csv2("Simulations/Data/db_rep_100_M_20_N_30_time_TRUE_hp_FALSE.csv")
tableTF$ID = as.character(tableTF$ID)
tableTF$ID_dataset = as.character(tableTF$ID_dataset)
tableFT = read_csv2("Simulations/Data/db_rep_100_M_20_N_30_time_FALSE_hp_TRUE.csv")
tableFT$ID = as.character(tableFT$ID)
tableFT$ID_dataset = as.character(tableFT$ID_dataset)
tableFF = read_csv2("Simulations/Data/db_rep_100_M_20_N_30_time_FALSE_hp_FALSE.csv")
tableFF$ID = as.character(tableFF$ID)
tableFF$ID_dataset = as.character(tableFF$ID_dataset)

tableM_0to20_TT = read_csv2("Simulations/Data/db_rep_100_M_0to20_N_30_time_TRUE_hp_TRUE.csv")
tableM_0to20_TT$ID = as.character(tableM_0to20_TT$ID)
tableM_0to20_TT$ID_dataset = as.character(tableM_0to20_TT$ID_dataset)
tableM_0to20_TF = read_csv2("Simulations/Data/db_rep_100_M_0to20_N_30_time_TRUE_hp_FALSE.csv")
tableM_0to20_TF$ID = as.character(tableM_0to20_TF$ID)
tableM_0to20_TF$ID_dataset = as.character(tableM_0to20_TF$ID_dataset)
tableM_0to20_FT = read_csv2("Simulations/Data/db_rep_100_M_0to20_N_30_time_FALSE_hp_TRUE.csv")
tableM_0to20_FT$ID = as.character(tableM_0to20_FT$ID)
tableM_0to20_FT$ID_dataset = as.character(tableM_0to20_FT$ID_dataset)
tableM_0to20_FF = read_csv2("Simulations/Data/db_rep_100_M_0to20_N_30_time_FALSE_hp_FALSE.csv")
tableM_0to20_FF$ID = as.character(tableM_0to20_FF$ID)
tableM_0to20_FF$ID_dataset = as.character(tableM_0to20_FF$ID_dataset)

tableM_20to200_TT = read_csv2("Simulations/Data/db_rep_100_M_20to200_N_30_time_TRUE_hp_TRUE.csv")
tableM_20to200_TT$ID = as.character(tableM_20to200_TT$ID)
tableM_20to200_TT$ID_dataset = as.character(tableM_20to200_TT$ID_dataset)
tableM_20to200_TF = read_csv2("Simulations/Data/db_rep_100_M_20to200_N_30_time_TRUE_hp_FALSE.csv")
tableM_20to200_TF$ID = as.character(tableM_20to200_TF$ID)
tableM_20to200_TF$ID_dataset = as.character(tableM_20to200_TF$ID_dataset)
tableM_20to200_FT = read_csv2("Simulations/Data/db_rep_100_M_20to200_N_30_time_FALSE_hp_TRUE.csv")
tableM_20to200_FT$ID = as.character(tableM_20to200_FT$ID)
tableM_20to200_FT$ID_dataset = as.character(tableM_20to200_FT$ID_dataset)
tableM_20to200_FF = read_csv2("Simulations/Data/db_rep_100_M_20to200_N_30_time_FALSE_hp_FALSE.csv")
tableM_20to200_FF$ID = as.character(tableM_20to200_FF$ID)
tableM_20to200_FF$ID_dataset = as.character(tableM_20to200_FF$ID_dataset)

##### TRAIN ALL MODEL ####
# db_to_train = tableGPFDA
# t1 = Sys.time()
# train_loop = loop_training(db_to_train, prior_mean = 0, ini_hp = list('theta_0' = c(1,1), 'theta_i' = c(1, 1, 0.2)),
#                            kern_0 = kernel_mu, kern_i = kernel, diff_M = F, common_times = T, common_hp = T)
# t2 = Sys.time()
# train_loop[['Time_train_tot']] = t2 - t1
# 
# saveRDS(train_loop, 'Simulations/Training/train_GPFDA_TT.rds')

##### RESULTS : evaluation of pred  ####
# train_loop = readRDS('Simulations/Training/train_FT.rds')
# tab_pred = tableGPFDA
# 
# res_pred = loop_pred(tab_pred, train_loop, 20, 10, F)
# 
# write.csv2(res_pred, "Simulations/Results/res_pred_GPFDA_TT.csv", row.names=FALSE)
# 
# res_pred %>% group_by(Method) %>% summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE)
# ggplot(res_pred) + geom_boxplot(aes(x = Method, y = MSE, fill = Method)) + scale_y_continuous(limits = c(0,200))


##### RESULTS : evaluation of mu_0  ####
# train_loop = readRDS('Simulations/Training/train_TT.rds')
# tab_mu = tableGPFDA
# 
# res_mu = eval_mu(tab_mu, train_loop)
# 
# write.csv2(res_mu, "Simulations/Results/res_GPFDA_TT.csv", row.names=FALSE)
# 
# res_mu %>% group_by(Method) %>% summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE)
# ggplot(res_mu) + geom_boxplot(aes(x = Method, y = MSE, fill = Method)) #+ scale_y_continuous(limits = c(0,100))
# 

##### RESULTS : pred with varying values of N* #####
# train_loop = readRDS('Simulations/Training/train_TT.rds')
# tab = tableTT
# 
# res = simu_var_N(tab, train_loop, nb_obs_max = 20, nb_test = 10)
# 
# write.csv2(res, "Simulations/Results/res_pred_N20-10_M20_TT.csv", row.names=FALSE)
# 
# res %>% group_by(Method) %>% summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE)
# ggplot(res) + geom_boxplot(aes(x = as.factor(N), y = MSE, fill = Method)) + scale_y_continuous(limits = c(0,100))

##### RESULTS : pred with varying values of M ####

#train_loop = readRDS('Simulations/Training/train_M_0to20_TT.rds')
# tab_M = tableM_0to20_TT #%>% filter(nb_M != 21)
# 
# res_M = loop_pred(tab_M, train_loop, nb_obs = 20, nb_test = 10, diff_M = T)
# 
# write.csv2(res_M, "Simulations/Results/res_pred_N20-10_M0to20_TT.csv", row.names=FALSE)
# 
# res_M %>% group_by(Method) %>% summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE)
# ggplot(res_M) + geom_boxplot(aes(x = as.factor(M), y = MSE, fill = Method)) + scale_y_continuous(limits = c(0,70))


##### RESULTS : mu_0 with varying M ####

#train = readRDS('Simulations/Training/train_M_0to20_TT.rds')
# tab_mu_M = tableM_0to20_TT #%>% filter(nb_M != 21)
# 
# res_mu_M = eval_mu_M(tab_mu_M, train_loop)
# write.csv2(res_mu_M, "Simulations/Results/res_mu_N20-10_M0to20_TT.csv", row.names=FALSE)
# 
# res_mu_M %>% group_by(Method, M) %>% summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE)
# ggplot(res_mu_M) + geom_boxplot(aes(x = as.factor(M), y = MSE, fill = Method)) + scale_y_continuous(limits = c(0,10))


##### PLOT OF RESULTS #### 

# table_res = res5  %>% dplyr::select(Time_train, Time_pred, Setting) %>% group_by(Setting) %>%
#    summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE) %>% mutate_if(is.numeric, round, 1) %>% 
#    mutate(Train = paste0(Time_train_Mean,' (', Time_train_SD, ')'), Pred = paste0(Time_pred_Mean,' (', Time_pred_SD, ')')) %>%
#   dplyr::select(Setting, Train, Pred) 
#
# res_plot = read_csv2("Simulations/Results/res_mu_N20-10_M20_FF.csv")
# table_res = res_plot %>% group_by(Method) %>%
#   summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE) %>% mutate_if(is.numeric, round, 1)
#
# write_csv2(table_res, 'Simulations/Table/table_pred_FF.csv')
# 
# ggplot(res_plot) + geom_boxplot(aes(x = as.factor(N), y = MSE, fill = Method)) 
# ggplot(res_plot) + geom_boxplot(aes(x = as.factor(M), y = MSE, fill = Method)) +
#   facet_wrap( ~ as.factor(M), scales="free") +
#   scale_y_continuous(limits = c(0,50))

##### TESTS SIMU ####
# bla_db = datasets_multi_N(rep = 10, M = 21, N = 30, G = seq(0, 10, 0.05), common_times = T,
#                           common_hp = F, kern_0 = kernel_mu, kern_i = kernel, int_mu_a = c(0,5), int_mu_b = c(0,2),
#                           int_i_a = c(0,5), int_i_b = c(0,2), int_i_sigma = c(0,1))

# bla_db2 = datasets_multi_M(rep = 10, vec_M = c(3, 11, 21), N = 30, G = seq(0, 10, 0.05), common_times = T,
#                           common_hp = F, kern_0 = kernel_mu, kern_i = kernel, int_mu_a = c(0,5), int_mu_b = c(0,2),
#                           int_i_a = c(0,5), int_i_b = c(0,2), int_i_sigma = c(0,1))
# 
# test_loop_train = loop_training(bla_db, prior_mean = 0, ini_hp = list('theta_0' = c(1,1), 'theta_i' = c(1, 1, 0.2)),
#                                 kern_0 = kernel_mu, kern_i = kernel, diff_M = F, common_times = T, common_hp = T)

# test_loop_pred = loop_pred(bla_db, test_loop_train, nb_obs =20, nb_test = 10)
# 
# test_res = simu_var_N(bla_db, nb_obs_max = 20, nb_test = 10, prior_mean = 0,
#                       ini_hp = list('theta_0' = c(1,1), 'theta_i' = c(1, 1, 0.2)), 
#                       kern_0 = kernel_mu, kern_i = kernel, common_times = T, common_hp = T)
# test_res %>% group_by(Method) %>% summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE)
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


##### TEST GPFDA ####
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
