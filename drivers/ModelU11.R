

ModelHeston<- function(x,x0,del,param,args){
  
  output <- 0
  if (args$mode=='option')
  {
    e<- 1e-5
    S<-exp(x[1])
    rho   <- param[1]
    kappa <- param[2]
    theta <- param[3]
    sigma <- param[4]
    
    
    objfun<-function(vega){
      objfun<- abs(HestonCOS(S,S,T_0,rate,q,sigma,kappa,theta,vega,rho,args$callput)-x[3]) 
    }
    # Infer v_0 from observed option prices using a root finding method.
    res<- nloptr( x0=v_0, 
                  eval_f=objfun, 
                  lb = c(0.01), 
                  ub = c(5), 
                  opts = list("algorithm"="NLOPT_LN_COBYLA", "maxeval" = 50, "xtol_rel" = args$eps, "print_level"=0))
    
    x[2] < - res$solution # the implied volatility
    
    # calculate the vega to obtain the Jacobian
    dVdv0 <- HestonCOS_vega(S,S,T_0,rate,q,sigma,kappa,theta,x[2],rho,args$callput)
    
    J <- dVdv0 
    if (is.nan(log(J))){
      print('Warning: NAN occured')
      print(J)
      print(x[2])
    }
    else
      output <- -log(J)
  }
  
  
  
  a        <- param[1]  # a
  b        <- param[2]  # b
  g        <- param[3]  # f
  d        <- param[4]  # d
  
  param_prime <- c(a,b,g)
  output <- output + ModelU11(x,x0,del,param_prime)
  
  return(output)
}


ModelU11 <- function(x,x0,del,param)
{

  a <- param[1]
  b <- param[2]
  f <- param[3]
  d <- param[4]
  
  sx <- f + d*x 
  y <- log(1 + (d*x)/f)/d 
  y0 <- log(1 + (d*x0)/f)/d 
  
  E <- exp(1) 
  
  cYm1 <- (-(1/2))*(y - y0)^2  
  
  cY0 <- (E^((-d)*y) - E^((-d)*y0))*((b*f - a*d)/(d^2*f)) + (y - y0)*((2*b - d^2)/(2*d))  
  
  if (y != y0)
  {
    cY1 <- (1/(2*d))*(b^2/(2*d^2) + a^2/(2*f^2) - (a*b)/(d*f))*((E^(-2*d*y) - E^(-2*d*y0))/(y - y0)) + 
      ((a*b)/(d^2*f) - b^2/d^3 + b/d - a/f)*((E^((-d)*y) - E^((-d)*y0))/(y - y0)) - (2*b - d^2)^2/(8*d^2) 
  }
   
  else
  {
    cY1 <- -((a*d - b*f)^2/(E^(2*d*y)*(2*f^2*d^2))) + 
     ((b - d^2)*((-a)*d + b*f))/(E^(d*y)*(f*d^2)) - (2*b - d^2)^2/(8*d^2)  
  }
   
  
  output <- (-(1/2))*log(2*pi*del) - log(sx) + cYm1/del + cY0 + cY1*del  
  
  print(output)
  return(output)
}

# ModelU11(3,4,1/52,c(0.3,0.4,0.5,0.6))
