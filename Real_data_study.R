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
     filter(Timestamp > age_min, Timestamp < age_max) %>% group_by(ID) %>% filter(n() > 4)
## If you need a subset of db
db = db %>% ungroup() %>% filter(ID %in% unique(.$ID)[1:200])
 
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

db_train = db_m_train
db_test = db_m_test

model_train = training(db_train, 0, ini_hp, kernel_mu, kernel, common_hp = T)

floop = function(i)
{ print(i)
  db_obs_i = db_test %>% filter(ID == i) %>% filter(Observed == 1)
  db_pred_i = db_test %>% filter(ID == i) %>% filter(Observed == 0)
  t_i_pred = db_pred_i %>% pull(Timestamp)

  res_algo = full_algo(db_train, db_obs_i, t_i_pred, kern_i = kernel, common_hp = T, plot = F, prior_mean = 0, kern_0 = kernel,
                       list_hp = model_train$hp, mu = NULL, ini_hp = ini_hp, hp_new_i = NULL)$Prediction

  pred_algo = res_algo$Mean 
  sd_algo = res_algo$Var %>% sqrt()
  tibble('ID' = i, 
         'MSE' = loss(pred_algo, db_pred_i$Output) %>% MSE(),
         'Ratio_IC' = ratio_IC(db_pred_i$Output, pred_algo - 1.96 * sd_algo, pred_algo + 1.96 * sd_algo)) %>% 
  return()
}
res_test = db_test$ID %>% unique() %>% lapply(floop)
db_res = do.call('rbind', res_test)
db_res %>% select(-ID) %>% summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE)
