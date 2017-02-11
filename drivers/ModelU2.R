
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
  
  
  
  a        <- param[1]  # theta
  b        <- param[2]  # kappa
  g        <- param[3]  # sigma
  
  
  param_prime <- c(a,b,g)
  output <- output + ModelU2(x,x0,del,param_prime)
  
  return(output)
}


ModelU2 <- function(x,x0,del,param)
  
{
  a <- param[1]
  b <- param[2]
  d <- param[3]
  
  y <- log(x)/d
  y0 <- log(x0)/d
  
  
  E <- exp(1) 
  sx <- d*x
  
  cYm1 <- (-(1/2))*(y - y0)^2
  cY0 <- (E^((-d)*y) - E^((-d)*y0))*(-(a/d^2)) + (y - y0)*(b/d - d/2) 
  
  if (y != y0)
  {
    cY1 <- (a^2/(4*d^3))*((E^(-2*d*y) - E^(-2*d*y0))/(y - y0)) + ((a*b)/d^3 - a/d)*((E^((-d)*y) - E^((-d)*y0))/(y - y0)) -  (2*b - d^2)^2/(8*d^2)  
    cY2 <- (-(a^2/(2*d^3)))*((E^(-2*d*y) - E^(-2*d*y0))/(y - y0)^3) + 
      ((2*a)/d - (2*a*b)/d^3)*((E^((-d)*y) - E^((-d)*y0))/(y - y0)^3) + (-(a^2/(2*d^2)))*((E^(-2*d*y) + E^(-2*d*y0))/(y - y0)^2) + 
      (a - (a*b)/d^2)*((E^((-d)*y) + E^((-d)*y0))/(y - y0)^2)  
  }
  else
  {
    cY1 <- (-4*a^2 - 8*a*(b - d^2)*E^(d*y) - (-2*b + d^2)^2*E^(2*d*y))/(E^(2*d*y)*(8*d^2))
    cY2 <- ((1/6)*a*(-2*a + (-b + d^2)*E^(d*y)))/E^(2*d*y) 
  }
  
  output <- (-(1/2))*log(2*pi*del) - log(sx) + cYm1/del + cY0 +  cY1*del + cY2*(del^2/2)

  print(output)
  return(output)
}

# ModelU2(3,4,1/25,c(0.1,0.2,0.3))