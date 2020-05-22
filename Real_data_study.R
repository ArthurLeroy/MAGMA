source('Algo_multitask_GP.R')

#### SPLITING TRAINING AND TESTING FUNCTIONS ####
split_train = function(db, ratio_train)
{ ## db : the database of all observed data
  ## ratio_train : number between 0 and 1 indicating the proportion of individuals in the training set
  #
  ## return : orginal db with the repartition of individuals between the training and testing set
  n_indiv = db$ID%>% unique()
  n_index = (ratio_train * length(n_indiv)) %>% round()
  index_train = sample(n_indiv, n_index, replace = F)
  
  db %>% mutate(Training = ifelse(ID %in% index_train, 1,0)) %>% return()
}

split_times = function(db, prop_test)
{ ## db : the database of all observed data
  ## prop_test : number between 0 and 1 indicating the proportion of points we predict on
  #
  ## return : orginal db with the repartition of individuals between the training and testing set
  db = db %>% rowid_to_column("Row")
  row_obs = db %>% group_by(ID) %>% top_n(Timestamp, n = - n() * (1 - prop_test)) %>% pull(Row)
  
  db %>% ungroup() %>% mutate(Observed = ifelse(db$Row %in% row_obs, 1, 0)) %>% dplyr::select(- Row) %>% return()
}

#### EVALUATION FUNCTIONS ####
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

#### RAW DATA TO CLEAN DATABASE ####

# raw_db = read_delim('C:/Users/user/Google Drive/Travail/Data thèse/001 BaseParEpreuve/100 Nage Libre.csv', delim = ';',
#                     col_types = cols_only(Identifiant = col_character(), Age = col_double(),
#                                           INDIVIDU_GENRE = col_integer(), COMPETITION_BASSIN =  col_integer())) %>%
#          bind_cols(read_csv2('C:/Users/user/Google Drive/Travail/Data thèse/001 BaseParEpreuve/100 Nage Libre.csv',
#                                col_types = cols_only(TEMPS  = col_double())))

raw_db = read_csv2('Simulations/Data/db_100m_freestyle.csv')

## Set the width of the observational grid to reduce the size of matrices
size_grid = 1000
age_min = 10
age_max = 20
## Set the approximation rounding step to result in a observational grid of the correct size
round_step = (age_max - age_min) / size_grid

db = raw_db %>% filter(COMPETITION_BASSIN == 50) %>%
     transmute(ID = Identifiant, Timestamp = Age, Output = TEMPS, Gender = INDIVIDU_GENRE) %>% 
     mutate(Timestamp = Timestamp %>% plyr::round_any(round_step)) %>% 
     group_by(ID, Timestamp) %>% summarise_all(mean) %>% 
     filter(Timestamp > age_min, Timestamp < age_max) %>% group_by(ID) %>% filter(n() > 5)
## If you need a subset of db
#db = db %>% ungroup() %>% filter(ID %in% unique(.$ID)[1:200])
 
## Split by gender and draw the training/testing sets
db_m = db %>% filter(Gender == 1) %>% dplyr::select(- Gender) %>% split_train(ratio_train = 0.6)
db_f = db %>% filter(Gender == 2) %>% dplyr::select(- Gender) %>% split_train(ratio_train = 0.6)

### Few info on the databases

## Number of individual in the db 
db_m$ID %>% n_distinct  
db_f$ID %>% n_distinct

## Number of different timestamps in the db 
db_m$Timestamp %>% n_distinct  
db_f$Timestamp %>% n_distinct

## Mean number of Timestamps per individual
db_m %>% count(ID) %>% pull(n) %>% mean
db_f %>% count(ID) %>% pull(n) %>% mean

## Max number of Timestamps per individual
db_m %>% count(ID) %>% pull(n) %>% max
db_f %>% count(ID) %>% pull(n) %>% max

#### Spliting training and testing sets and select observed times
db_m_train = db_m %>% filter(Training == 1)
db_m_test = db_m %>% filter(Training == 0) %>% split_times(prop_test = 0.2)
db_f_train = db_f %>% filter(Training == 1)
db_f_test = db_f %>% filter(Training == 0) %>% split_times(prop_test = 0.2)

db_train = db_f_train
db_test = db_f_test

# model_train = training(db_train, 0, ini_hp, kernel_mu, kernel, common_hp = T)
# saveRDS(model_train 'Simulations/Training/train_real_data_female_TT.rds')
# list_ID = model_train$hp$theta_i %>% names
# db_train = db_f %>% filter(ID %in% list_ID)
# db_test = db_f %>% filter(ID %notin% list_ID) %>% split_times(prop_test = 0.2)
post_mu = posterior_mu(db_train, db_train, db_f$Timestamp, 0, kernel_mu, kernel, model_train$hp)
floop = function(i)
{
  print(i)
  db_obs_i = db_test %>% filter(ID == i) %>% filter(Observed == 1)
  db_pred_i = db_test %>% filter(ID == i) %>% filter(Observed == 0)
  t_i_pred = db_pred_i %>% pull(Timestamp)

  res_algo = full_algo(db_train, db_obs_i, t_i_pred, kern_i = kernel, common_hp = T, plot = F, prior_mean = 0, kern_0 = kernel_mu,
                       list_hp = model_train$hp, mu = post_mu, ini_hp = ini_hp, hp_new_i = NULL)$Prediction

  hp_one_gp = train_new_gp(db_obs_i, rep(0, nrow(db_obs_i)), cov_mu = 0, ini_hp$theta_i, kern_i = kernel)
  ## Prediction for one GP model
  res_one_gp = pred_gp(db_obs_i, t_i_pred, 0, cov_mu = NULL, kernel, hp_one_gp$theta, hp_one_gp$sigma)
  
  pred_algo = res_algo$Mean 
  sd_algo = res_algo$Var %>% sqrt()
  
  pred_one_gp = res_one_gp$Mean 
  sd_one_gp = res_one_gp$Var %>% sqrt()
  
  eval_one_gp = tibble('ID' = i,
                       'MSE' = loss(pred_one_gp, db_pred_i$Output) %>% MSE(),
                       'Ratio_IC' = ratio_IC(db_pred_i$Output, pred_one_gp - 1.96 * sd_one_gp, pred_one_gp + 1.96 * sd_one_gp))
  
  eval_algo = tibble('ID' = i, 
                     'MSE' = loss(pred_algo, db_pred_i$Output) %>% MSE(),
                     'Ratio_IC' = ratio_IC(db_pred_i$Output, pred_algo - 1.96 * sd_algo, pred_algo + 1.96 * sd_algo))
  
  rbind(eval_algo, eval_one_gp) %>% mutate(Method = c('Algo', 'One GP')) %>% return()
}
res_test = db_test$ID %>% unique() %>% lapply(floop)
db_res = do.call('rbind', res_test)
db_res %>% select(-ID) %>% group_by(Method) %>% summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE) %>% 
write_csv2( 'Simulations/Table/res_pred_realdata_women.csv')

### Test on an individual 
<<<<<<< HEAD
indiv = 'LAUMOND Émilie 22/02/2000'

pred_example = full_algo(db_train,(db_test %>% filter(ID == indiv))[1:4,] , seq(10, 20, 0.01), kernel,
                         common_hp = T, plot = F, prior_mean = 0, kernel, list_hp = model_train$hp, mu = NULL,
=======
indiv = unique(db_test$ID)[99]
db_obs = (db_test %>% filter(ID == indiv))

post_mu = posterior_mu(db_train, db_obs[1:6,], seq(10, 20, 0.01), 0, kernel_mu, kernel, model_train$hp)

pred_example = full_algo(db_train, db_obs[1:6,], seq(10, 20, 0.01), kernel,
                         common_hp = T, plot = F, prior_mean = 0, kernel, list_hp = model_train$hp, mu = post_mu,
>>>>>>> ef1ba1098f4a93449ae3451edeb28268e62567cc
                         ini_hp = ini_hp, hp_new_i = NULL)

plot_gp(pred_example$Prediction, data_train = db_train %>% filter(ID %in% unique(db_train$ID)[1:50]), data = db_test %>% filter(ID == indiv),
        mean = pred_example$Mean_process$pred_GP) + guides(color = FALSE) + 
  geom_point(data = db_obs[6:nrow(db_obs),], aes(Timestamp, Output), color ='blue') 

hp_one_gp = train_new_gp((db_test %>% filter(ID == indiv))[1:4,], 0, 0, ini_hp$theta_i, kernel)
pred_gp((db_test %>% filter(ID == indiv))[1:4,], timestamps =  seq(10, 20, 0.01), mean_mu = 0, cov_mu = NULL, 
            kern = kernel, theta = hp_one_gp$theta, sigma = hp_one_gp$sigma) %>% 
       plot_gp(data_train = db_train, data = db_test %>% filter(ID == indiv),
          mean = pred_example$Mean_process$pred_GP) + guides(color = FALSE)
  
