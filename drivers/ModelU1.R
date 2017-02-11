
source('C:\\Users\\smitha\\Documents\\R\\MLEMVD-master\\MLEMVD-master\\HestonFourierCosine.R')
require(nloptr)
require(pracma)

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
  
  
  
  a        <- param[1]  # theta
  b        <- param[2]  # kappa
  g        <- param[3]  # sigma
  
  
  param_prime <- c(a,b,g)
  output <- output + ModelU1(x,x0,del,param_prime)
  
  return(output)
}

ModelU1 <- function(x,x0,del,param,args){
  
  a <- param[1] 
  b <- param[2] 
  sigma <- param[3] 
  
  y <- 2/(sqrt(x)*sigma) 
  y0 <- 2/(sqrt(x0)*sigma) 
  sx <- sigma*x^(3/2) 
  
  output <- -log(2*pi*del)/2 - log(sx)   
  + (-(1/2))*(y - y0)^2 / del   
  + (a*(-y^2 + y0^2)* sigma^2 + (-8*b + 6* sigma^2)*log(y/y0))/(4* sigma^2)  
  + ( -((1/(24*y*y0* sigma^4))* (48*b^2 + 24*b*(-2 + a*y*y0)* sigma^2 + (9 - 24*a*y*y0 + a^2*y*y0*(y^2 + y*y0 + y0^2))* sigma^4)) )*del   
  + ( -((48*b^2 - 48*b* sigma^2 + (9 + a^2*y^2*y0^2)* sigma^4)/ (24*y^2*y0^2* sigma^4)) )*(del^2/2) 
  
  print(output)
  return(output) 
  
}

# ModelHeston(3,4,1/52,c(0.3,0.4,0.5))

