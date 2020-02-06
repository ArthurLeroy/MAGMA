setwd(dir = 'C:/Users/user/CloudStation/Maths/These/Processus Gaussiens/Code R/Algo multitask GP')
source('Algo_multitask_GP.R')

RMSE = function(m, o){sqrt(mean((m - o)^2))}

library(mvtnorm)

########### Package gptk####
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
# model = gpCreate(dim(x)[2], dim(y)[2], x, y, options)	
# model = gpOptimise(model, display=TRUE, iters=400)
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
deriv1 = function(t1, t2, theta = list(1, 0.5), sigma = 0.2)
{
  a = theta[[1]]
  b = theta[[2]]
  
  return( exp(a) * exp((-1/(2 * exp(b))) * t((t1 - t2)) %*% (t1 - t2)) )
}

deriv2 = function(t1, t2, theta = list(1, 0.5), sigma = 0.2)
{
  a = theta[[1]]
  b = theta[[2]]
  resi = t(t1 - t2) %*% (t1 - t2)
  norm = 0.5 * exp(-b)
  
  return( exp(a) * norm * resi * exp(- norm * resi) )
}

deriv3 = function(t1, t2, theta = list(1, 0.5), sigma = 0.2)
{
  if(all(t1 == t2))
  {
    return(sigma^2)
  }
  else{return(0)}
}

gr_fn = function(x, db = db_test)
{
  y = db$Output
  t = db$Timestamp
  inv = kern_to_inv(t, kernel, theta = x[1:2], sigma = x[3])
  
  cste_term = inv %*% y %*% t(inv %*% y) - inv
  
  g_1 = (-1/2 * cste_term %*% kern_to_cov(t, deriv1, theta = x[1:2], sigma = x[3])) %>%  diag() %>% sum()
  g_2 = (-1/2 * cste_term %*% kern_to_cov(t, deriv2, theta = x[1:2], sigma = x[3])) %>%  diag() %>% sum()
  g_3 = (-1/2 * cste_term %*% kern_to_cov(t, deriv3, theta = x[1:2], sigma = x[3])) %>%  diag() %>% sum()
  grad = c(g_1, g_2, g_3)
  
  return(grad)
}


############ logL of the GP ##############
fn<- function (x, db = db_test, mean = 0, kern = kernel) 
{
  #y = db$Output
  #t = db$Timestamp
  #n = nrow(db)
  #inv = kern_to_inv(t, kernel, theta = x)
  #return(-(t(y) %*% inv %*% y)/2 + log(det(inv))/2 - log(2*pi)*n/2) 
  return(- dmvnorm(db$Output, rep(0, nrow(db)), 
                   kern_to_inv(db$Timestamp, kern, theta = x[1:2], sigma = x[3]), log = T))
}


t = 10:20
db_test = tibble('Timestamp' = t, 
                 'Input' = paste0('X', t),
                 'Output' = rmvnorm(1, rep(0,length(t)), kern_to_cov(t, kernel, theta = list(1, 0.5), sigma = 0.2)) %>% as.vector())

meth0 <- c("Nelder-Mead", "L-BFGS-B")
st0 <- c(1, 0.5, 0.2)  # the standard start
result0 <- opm(st0, fn, method = meth0, control = list(kkt = FALSE))
result0 = summary(result0, order=value)
result0

pred_gp(db_test, timestamps = seq(9,21, 0.03), mean = 0,
        theta = c(result0$p1[1], result0$p2[1]), sigma = result0$p3[1]) %>% plot_gp(data = db_test)

