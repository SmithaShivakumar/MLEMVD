source('C:\\Users\\smitha\\Documents\\R\\RModels\\logdensity2loglik.R')
source('C:\\Users\\smitha\\Documents\\R\\RModels\\mymle.R')
source('C:\\Users\\smitha\\Documents\\R\\RModels\\Univariate\\ModelU1.R')
require(nloptr)
require(numDeriv)

args<-list(maxiter=50, eps=1e-8, print_level=3)
eps <-1e-8
#k1,rho, kappa, theta, sigma
args$l <- c(eps, -1.0 + eps,eps,eps, eps)
args$u <- c(1.0-eps,1.0-eps,4.0-eps,1.0-eps,2.0-eps)
#k1,rho, kappa, theta, sigma
eval_g_ineq <- function (x) {
  grad <- c(0,0, -2.0*x[4],-2.0*x[3],2.0*x[5])
  return(list("constraints"=c(x[5]*x[5] - 2.0*x[3]*x[4]), "jacobian"=grad))  
}
# Model = B11;
# dx1 = (k1 + k2*x2)*dt + sqrt(x2)*(sqrt(1 - rho^2)*dW1 + rho*dW2)
# dx2 = kappa*(theta - x2)*dt + sigma*x2*dW2
# 5 parameters to be estimated: (k1,k2,rho,kappa,theta,sigma)

# Heston model:
# dln(S_t) = \mu dt + \sqrt{V_t}dW_t^{1}
# dV_t = \kappa(\theta - V_t)dt  + \sigma \sqrt{V_t}dW_t^{2}
# k2 = 0 for Heston model


# starting values for MLE algorithm and simulated series
# Step 1: Simulating Heston data

S_0      <-50 
v_0      <-0.1  # the initial variance is assumed to be unknown
rate     <-0.1  # the risk free rate is assumed to be known
q        <-0

a_0      <- rate -q
b_0      <- -0.5
rho_0    <- -0.8
kappa_0  <-3
theta_0  <-0.2
sigma_0  <-0.25 

param_0<-c(v_0,rho_0,kappa_0,theta_0,sigma_0)

args$mode = 'state' # calibration to ATM option prices 
# daily data: del = 1/52
delta <- 1/252
n     <- 500
# s_t=ln(S_t)
x1_0 <- log(S_0) #s_0 <- ln S_0;
x2_0 <- v_0
factor  <- 10
nsimul  <- factor*n
delsimul<- delta/factor

rndn1 <- rnorm(nsimul)
rndn2 <- rnorm(nsimul)

x1simul <- rep(0, nsimul)   
x2simul <- rep(0, nsimul)  

x1simul[1] <- x1_0 
x2simul[1] <- x2_0
# calibrate to option prices
if(args$mode=='state'){
  T_0  <- 0.1
  callput <-'C'
}



for (i in 2:nsimul){
  # dx2 = kappa*(theta - x2)*dt + sigma*sqrt(x2)*dW2
  x2simul[i] <- x2simul[i-1]+kappa_0*(theta_0 - x2simul[i-1])*delsimul + sigma_0*sqrt(x2simul[i-1])*sqrt(delsimul)*rndn2[i]
  # dx1 = (a_0 + b_0*x2)*dt + sqrt(x2)*(sqrt(1 - rho^2)*dW1 + rho*dW2)
  x1simul[i] <- x1simul[i-1]+(a_0 + b_0*x2simul[i-1])*delsimul + sqrt(x2simul[i-1])*(sqrt(1 - rho_0^2)*sqrt(delsimul)*rndn1[i] + rho_0*sqrt(delsimul)*rndn2[i])    
}

x   <- cbind(rep(0,n),rep(0,n))
x_0  <- c(x1_0,x2_0)
x[1,]<- x_0
for (i in 2:n){
  x[i,1] <- x1simul[1+(i-1)*factor]
  x[i,2] <- x2simul[1+(i-1)*factor]
}

output <- mymle(ModelU1, x, delta, param_0)
print(output)




# print('\nComputing Variance...\n')
# # Compute the var-cov matrix
# # The derivatives are computed numerically. The Hessian may not be
# # positive definite. We report the inverse[Infomation], as well as the
# # robust sandwich matrix.
# H <- hessian(objfun,param_0)
# InfoMatrix <- logdensity2info(logdensity,x,del,param_0)
# Variance <- solve(InfoMatrix)
# invH <- solve(H)
# Variance_Robust <- invH * InfoMatrix * invH
# print('Estimation finished.\n\n')
# 
# #Output
# # output.param = param_mle
# # output.variance = Variance
# # output.variance_robust = Variance_Robust
# # output.se = sqrt(diag(Variance))
# # output.se_robust = sqrt(diag(Variance_Robust))
# # output.objfun = objfun
# # output.loglik = -fval
# # output.exitflag = exitflag
# # Printing the MLE result
# kronecker(matrix(1,1,50),'*')
# print('%10s%10s%10s\n','COEF','SE','Robust SE')
# for (k in 1:param_mle)
# {
#   print('%10.4f%10.4f%10.4f\n',output.param(k),output.se(k),output.se_robust(k))
# }
# print(kronecker(matrix(1,1,50),'*'))
# print(output.exitmsg)
# # remove the search path
# # rmpath(genpath(cd))
# genpath <- getwd()
# setwd(new.dir)
#  ################################################################################
# # Helper functions
# logdensity2info <- function(logdensity,x,del,param){
#    # Compute the information matrix = Sum_{i=1}^n Score_i' Score_i
#   n <- length(x) - 1
#   output <- 0
#   for (i in 1:n)
#   {
#     tmpfun <- (theta)(logdensity(x(i+1),x(i),del,theta))
#     score_i <- gradest(tmpfun,param) # This is a row vector
#     output <- output + score_i * score_i
#   }
#   
#}