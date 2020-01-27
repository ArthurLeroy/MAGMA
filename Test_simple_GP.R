setwd(dir = 'C:/Users/user/CloudStation/Maths/These/Processus Gaussiens/Code R/Algo multitask GP')
source('Computing_functions.R')
source('Algo_multitask_GP.R')

RMSE = function(m, o){sqrt(mean((m - o)^2))}

########### Package gptk####
library(gptk)




## Inputs
x = db_test$Timestamp %>% as.matrix()
xtest = matrix(seq(9,21, 0.03),ncol=1)
## Generate y

options = gpOptions()
options$kern$comp = list("rbf","white")

trueKern = kernCreate(x, list(type="cmpnd",comp=list("rbf", "white")))
K = kernCompute(trueKern, x)

y = db_test$Output %>% as.matrix()
y = scale(y,scale=FALSE)

model = gpCreate(dim(x)[2], dim(y)[2], x, y, options)	
model = gpOptimise(model, display=TRUE, iters=400)

ll_opt = gpLogLikelihood(model)

meanVar = gpPosteriorMeanVar(model, xtest, varsigma.return=TRUE) ## GP mean and variance.

plot(x, y, xlim=c(9,21))
lines(xtest, meanVar$mu)
lines(xtest, meanVar$mu+sqrt(meanVar$varsigma), col="blue")
lines(xtest, meanVar$mu-sqrt(meanVar$varsigma), col="blue")

graphics.off();

########### Gradients of the logL GP #####
deriv1 = function(t1, t2, theta = list(1, 0.5, 0.5))
{
  a = theta[[1]]
  b = theta[[2]]
  #c = theta[[3]]
  
  return( 2 * a * exp((-1/(2 * b^2)) * t((t1 - t2)) %*% (t1 - t2)) )
}

deriv2 = function(t1, t2, theta = list(1, 0.5, 0.5))
{
  a = theta[[1]]
  b = theta[[2]]
  # c = theta[[3]]
  
  return( a^2 *  t((t1 - t2)) %*% (t1 - t2) * exp((-1/(2 * b^2)) * t((t1 - t2)) %*% (t1 - t2)) / (b^3) )
}

deriv3 = function(t1, t2, theta = list(1, 0.5, 0.5))
{
  if(all(t1 == t2))
  {
    return(2 * theta[[3]])
  }
  else{return(0)}
}

gr_fn = function(x, db = db_obs)
{
  y = db$Output
  t = db$Timestamp
  inv = kern_to_inv(t, kernel, theta = x)
  
  cste_term = inv %*% y %*% t(inv %*% y) - inv
  
  g_1 = (-1/2 * cste_term %*% kern_to_inv(t, deriv1, as.list(x))) %>%  diag() %>% sum()
  g_2 = (-1/2 * cste_term %*% kern_to_inv(t, deriv2, as.list(x))) %>%  diag() %>% sum()
  #g_3 = (-1/2 * cste_term %*% kern_to_inv(t, deriv3, as.list(x))) %>%  diag() %>% sum()
  
  grad = c(g_1, g_2)
  return(grad)
}


############ logL of the GP ##############
fn<- function (x, db = db_test, kern = kernel) 
{
  #y = db$Output
  #t = db$Timestamp
  #n = nrow(db)
  #inv = kern_to_inv(t, kernel, theta = x)
  #return(-(t(y) %*% inv %*% y)/2 + log(det(inv))/2 - log(2*pi)*n/2) 
  return(- dmvnorm(db$Output, rep(0, nrow(db)), 
                   kern_to_cov(db$Timestamp, kern, theta = x), log = T))
}


t = 10:20
db_test = tibble('Timestamp' = t, 
                 'Input' = paste0('X', t),
                 'Output' = rmvnorm(1, rep(0,length(t)), kern_to_cov(t, kernel, theta = c(2,1, 0.5))) %>% as.vector())

meth0 <- c("Nelder-Mead", "BFGS", "L-BFGS-B")
st0 <- c(5, 5, 0.5) %>% log() # the standard start
result0 <- opm(st0, fn, method= meth0)
result0 = summary(result0, order=value)
result0

pred_gp(db_test, timestamps = seq(9,21, 0.03), theta = c(result0$p1[1], result0$p2[1], result0$p3[1])) %>% 
  plot_gp(data = db_test)


