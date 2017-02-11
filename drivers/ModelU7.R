
ModelHeston<- function(x,x0,del,param,args){
  
  output <- 0
  if (args$mode=='option')
  {
    e<- 1e-5
    S<-exp(x[1])
    rho   <- param[1]
    kappa <- param[2]
    theta <- param[3]
    sigmaMH <- param[4]
    
    
    objfun<-function(vega){
      objfun<- abs(HestonCOS(S,S,T_0,rate,q,sigmaMH,kappa,theta,vega,rho,args$callput)-x[3]) 
    }
    # Infer v_0 from observed option prices using a root finding method.
    res<- nloptr( x0=v_0, 
                  eval_f=objfun, 
                  lb = c(0.01), 
                  ub = c(5), 
                  opts = list("algorithm"="NLOPT_LN_COBYLA", "maxeval" = 50, "xtol_rel" = args$eps, "print_level"=0))
    
    x[2] < - res$solution # the implied volatility
    
    # calculate the vega to obtain the Jacobian
    dVdv0 <- HestonCOS_vega(S,S,T_0,rate,q,sigmaMH,kappa,theta,x[2],rho,args$callput)
    
    J <- dVdv0 
    if (is.nan(log(J))){
      print('Warning: NAN occured')
      print(J)
      print(x[2])
    }
    else
      output <- -log(J)
  }
  
  
  b        <- param[1]  # kappa  
  a        <- param[2]  # alpha
  g        <- param[3]  # eta
  
  
  param_prime <- c(b,a,g)
  output <- output + ModelU7(x,x0,del,param_prime)
  
  return(output)
}


ModelU7 <- function(x,x0,del,param)
{
  #   Vasicek (Ornstein-Uhlenbeck) Model  
  # \(Mu)(x)=kappa (alpha-x); 
  # \(Sigma)(x)=eta;  
  # m = 1;
  
  kappa = param[1]
  alpha = param[2]
  eta = param[3]
  m = 1
  
  output = (-m/2)*log(2*pi*del) - log(eta)  
  -((x - x0)^2/(2*eta^2))/del  
  + ((-(x^2/2) + x0^2/2 + x*alpha - x0*alpha)*kappa)/eta^2  
  - ((1/(6*eta^2))*(kappa*(-3*eta^2 + (x^2 + x0^2 + x*(x0 - 3*alpha) - 3*x0*alpha + 3*alpha^2)*kappa)))*del 
  - (1/2)*(kappa^2/6)*del^2  
  + (1/6)*((4*x^2+7*x*x0+4*x0^2-15*x*alpha-15*x0*alpha+15*alpha^2)*kappa^4)/(60*eta^2)*del^3
  
  print(output)
  return(output)
  
}

# ModelU7(4,5,1/23,c(2,3,4))
