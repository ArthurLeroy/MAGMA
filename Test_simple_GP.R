setwd(dir = 'C:/Users/user/CloudStation/Maths/These/Processus Gaussiens/Code R/Algo multitask GP')
source('Algo_multitask_GP.R')

RMSE = function(m, o){sqrt(mean((m - o)^2))}

library(mvtnorm)

########## Package gptk####
# library(gptk)
# 
# ## Inputs
# x = db_test$Timestamp %>% as.matrix()
# xtest = matrix(seq(9,21, 0.03),ncol=1)
# ## Generate y
# 
# options = gpOptions()
# options$kern$comp = list("rbf","white")
# 
# trueKern = kernCreate(x, list(type="cmpnd",comp=list("rbf", "white")))
# K = kernCompute(trueKern, x)
# 
# y = db_test$Output %>% as.matrix()
# y = scale(y,scale=FALSE)
# 
# t3 = Sys.time()
# model = gpCreate(dim(x)[2], dim(y)[2], x, y, options)
# model = gpOptimise(model, display=TRUE, iters = 400)
# t4 = Sys.time()
# print(t4-t3)
# 
# ll_opt = gpLogLikelihood(model)
# 
# meanVar = gpPosteriorMeanVar(model, xtest, varsigma.return=TRUE) ## GP mean and variance.
# 
# plot(x, y, xlim=c(9,21))
# lines(xtest, meanVar$mu)
# lines(xtest, meanVar$mu+sqrt(meanVar$varsigma), col="blue")
# lines(xtest, meanVar$mu-sqrt(meanVar$varsigma), col="blue")
# 
# graphics.off();

########### Gradients of the logL GP #####
deriv_hp1 = function(mat, theta)
{
  exp(theta[[1]] - exp(-theta[[2]])/2 * mat) %>% return()
}

deriv_hp2 = function(mat, theta)
{
  grad = 0.5 * exp(- theta[[2]]) * mat
  
  return( exp(theta[[1]]) * grad * exp(- grad) )
}

gr_fn = function(x, db, mean = 50, kern)
{ 
  y = db$Output
  t = db$Timestamp
  inv = kern_to_inv(t, kern, theta = x[1:2], sigma = x[[3]])
  prod_inv = inv %*% (y - mean) 
  cste_term = prod_inv %*% t(prod_inv) - inv
  
  g_1 = 1/2 * (cste_term %*% kern_to_cov(t, deriv_hp1, theta = x[1:2], sigma = 0)) %>%  diag() %>% sum()
  g_2 = 1/2 * (cste_term %*% as.matrix(deriv_hp2(dist(t)^2, theta = x[1:2]) ))  %>%  diag() %>% sum()
  g_3 = x[[3]] * (cste_term %>% diag() %>% sum() )
  grad = c(g_1, g_2, g_3)
  
  return(- grad)
}


############ logL of the GP ##############
fn<- function (x, db = db_test, mean = 50, kern) 
{
  return(- dmvnorm(db$Output, rep(mean, nrow(db)), 
                   kern_to_inv(db$Timestamp, kern, theta = x[1:2], sigma = x[[3]]), log = T))
}


t = seq(10, 20, 0.05)
db_test = tibble('ID' = '1',
                 'Timestamp' = t, 
                 'Input' = paste0('X', t),
                 'Output' = rmvnorm(1, rep(0,length(t)), kern_to_cov(t, kernel, theta = list(1, 0.5), sigma = 0.2)) %>% as.vector())

meth0 <- c("Nelder-Mead", 'L-BFGS-B')
st0 <- c(3,3,1)  # the standard start
result0 <- opm(st0, fn, db = db_test, gr = gr_fn, method = meth0, control = list(kkt = FALSE))
result0 = summary(result0, order=value)
result0

pred_gp(db_test, timestamps = seq(9,21, 0.03), mean = 0,
        theta = c(result0$p1[1], result0$p2[1]), sigma = result0$p3[1]) %>% plot_gp(data = db_test)

