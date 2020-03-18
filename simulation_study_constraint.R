simu_scheme2 = function(M = 10, N = 10, G = seq(0, 10, 0.05), kern_0 = kernel_mu, kern_i = kernel,
                       int_mu_a = 1-5,
                       int_mu_b = 1-5,
                       int_i_a = 0.10,
                       int_i_b = 0.10,
                       int_i_sigma = 1e-1)
{
  m_0 = prior_mean(G)
  mu_a = int_mu_a
  mu_b = int_mu_b
  
  db_0 = simu_indiv('0', G, m_0, kern_0, mu_a, mu_b, sigma = 0)
  
  floop = function(i)
  { 
    #t_i = sample(G, N, replace = F) %>% sort()
    t_i = G
    #i_a = draw(int_i_a)
    i_a = int_i_a
    #i_b = draw(int_i_b)
    i_b = int_i_b
    #i_sigma = draw(int_i_sigma)
    i_sigma = int_i_sigma
    mean_i = db_0 %>% filter(Timestamp %in% t_i) %>% pull(Output)
    
    simu_indiv(as.character(i), t_i, mean_i, kern_i, i_a, i_b, i_sigma) %>% return()
  }
  db_i = lapply(seq_len(M), floop)
  db = do.call("rbind", db_i) %>% rbind(db_0)
  return(db)
}


##### SIMULATION STUDY ##### 
db_train = simu_scheme2(M = 10, N = 10, G = seq(0, 10, 0.05))
plot_db(db_train)