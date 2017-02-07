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

S_0    <-50 
k1_0   <-0.05 
k2_0   <-0
rho_0  <- -0.8
kappa_0<-3
theta_0<-0.2
sigma_0<-0.25 
  
param_0 <- c(k1_0,rho_0,kappa_0,theta_0,sigma_0)
  
# daily data: del = 1/52
delta <- 1/252
n     <- 500
# s_t=ln(S_t)
x1_0 <- log(S_0) #s_0 <- ln S_0;
x2_0 <- theta_0; #v_0 <- theta_0;
  
factor  <- 10
nsimul  <- factor*n
delsimul<- del/factor
  
rndn1 <- rnorm(nsimul)
rndn2 <- rnorm(nsimul)
  
x1simul <- rep(0,nsimul)   
x2simul <- rep(0, nsimul)  
  
x1simul[1] <- x1_0 
x2simul[1] <- x2_0
  
for (i in 2:nsimul){
  # dx2 = kappa*(theta - x2)*dt + sigma*x2*dW2
  x2simul[i] <- x2simul[i-1]+kappa_0*(theta_0 - x2simul[i-1])*delsimul + sigma_0*x2simul[i-1]*sqrt(delsimul)*rndn2[i]
  
  # dx1 = (k1 + k2*x2)*dt + sqrt(x2)*(sqrt(1 - rho^2)*dW1 + rho*dW2)
  x1simul[i] <- x1simul[i-1]+(k1_0 + k2_0*x2simul[i-1])*delsimul + sqrt(x2simul[i-1])*(sqrt(1 - rho_0^2)*sqrt(delsimul)*rndn1[i] + rho_0*sqrt(delsimul)*rndn2[i])    
}
  
x    <- cbind(rep(0,n),rep(0,n))
x_0  <- c(x1_0,x2_0)
x[1,]<- x_0

for (i in 2:n){
  x[i,1] <- x1simul[1+(i-1)*factor]
  x[i,2] <- x2simul[1+(i-1)*factor]    
}
output <- mymle(ModelB11, x, del, param_0)