


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
  
  
  
  
  a        <- param[1]  # am1
  b        <- param[2]  # a0
  g        <- param[3]  # a1
  h        <- param[4]  # a2
  t        <- param[5]  # a3
  sig      <- param[6]  # sigma
  r        <- param[7]  # rho
  
  param_prime <- c(a,b,g,h,t,sig,r)
  output <- output + ModelU13(x,x0,del,param_prime)
  
  return(output)
}


ModelU13 <- function(x,x0,del,param)
{
  
  a0  <- param[1]
  am1  <- param[2]
  a1  <- param[3]
  a2  <- param[4]
  a3  <- param[5]
  sigma  <- param[6]
  rho  <- param[7]
  
  s_x  <- sigma*x^rho 
  
  y  <- x^(1 - rho)/((-1 + rho)*sigma)  
  y0  <- x0^(1 - rho)/((-1 + rho)*sigma)  
  
  cm1  <- (-(1/2))*(y - y0)^2  
  
  
  c0  <- (1/2)*(a1*(1 - rho)*y^2 - a1*(1 - rho)*y0^2 - (2*a2*(-1 + rho)^(1 + (-2 + rho)/(-1 + rho))*sigma^(1/(1 - rho))*       (y^(2 + 1/(1 - rho)) - y0^(2 + 1/(1 - rho))))/(-3 + 2*rho) - (a3*(-1 + rho)^(1 + (-3 + rho)/(-1 + rho))*(y^(2 - 2/(-1 + rho)) - y0^(2 - 2/(-1 + rho))))/      (sigma^(2/(-1 + rho))*(-2 + rho)) - (2*a0*(-1 + rho)^(1 + rho/(-1 + rho))*sigma^(1/(-1 + rho))*(y^(2 + 1/(-1 + rho)) - y0^(2 + 1/(-1 + rho))))/      (-1 + 2*rho) - (am1*(-1 + rho)^(1 + (1 + rho)/(-1 + rho))*sigma^(2/(-1 + rho))*(y^((2*rho)/(-1 + rho)) - y0^((2*rho)/(-1 + rho))))/rho +      (2*rho*log(y/y0))/(-2 + 2*rho))  
  
  
  c1  <- (-(1/(y - y0)))*(a1*(1 - rho)*y + (a1*(1 - rho)*rho*y)/(-2 + 2*rho) -      (1/3)*a3*am1*(-1 + rho)^((-3 + rho)/(-1 + rho) + (1 + rho)/(-1 + rho))*y^3 + (rho*(1/y - 1/y0))/(-2 + 2*rho) + (rho^2*(-(1/y) + 1/y0))/(2*(-2 + 2*rho)^2) +      (rho^2*(-(1/y) + 1/y0))/(2*(-1 + rho)*(-2 + 2*rho)) + (a1*(1 - rho)*(y - y0))/(2*(-1 + rho)) - (a1*(1 - rho)*rho*(y - y0))/(2*(-1 + rho)) -      a1*(1 - rho)*y0 - (a1*(1 - rho)*rho*y0)/(-2 + 2*rho) + (1/3)*a3*am1*(-1 + rho)^((-3 + rho)/(-1 + rho) + (1 + rho)/(-1 + rho))*y0^3 +      ((rho*y)/(-2 + 2*rho) - (rho*y0)/(-2 + 2*rho))/(2*y*y0 - 2*rho*y*y0) + (1/6)*a1^2*(1 - rho)^2*(y^3 - y0^3) +      (2/3)*a3*am1*(-1 + rho)^((-3 + rho)/(-1 + rho) + (1 + rho)/(-1 + rho))*(y^3 - y0^3) + (1/3)*a0*a2*(-1 + rho)^((-2 + rho)/(-1 + rho) + rho/(-1 + rho))*      sigma^(1/(1 - rho) + 1/(-1 + rho))*(y^3 - y0^3) + (a2*(-1 + rho)^((-2 + rho)/(-1 + rho))*sigma^(1/(1 - rho))*(y^(1 + 1/(1 - rho)) - y0^(1 + 1/(1 - rho))))/      (-2 + rho) - (a2*(-1 + rho)^((-2 + rho)/(-1 + rho))*rho*sigma^(1/(1 - rho))*(y^(1 + 1/(1 - rho)) - y0^(1 + 1/(1 - rho))))/(2*(-2 + rho)) -      (a2*(-1 + rho)^(1 + (-2 + rho)/(-1 + rho))*rho*sigma^(1/(1 - rho))*(y^(1 + 1/(1 - rho)) - y0^(1 + 1/(1 - rho))))/((-2 + rho)*(-2 + 2*rho)) -      
      (a1*a2*(1 - rho)*(-1 + rho)^(1 + (-2 + rho)/(-1 + rho))*sigma^(1/(1 - rho))*(y^(3 + 1/(1 - rho)) - 
  y0^(3 + 1/(1 - rho))))/(-4 + 3*rho) +      (a0*a3*(-1 + rho)^(1 + (-3 + rho)/(-1 + rho) + rho/(-1 + rho))*(y^(3 + 1/(1 - rho)) - y0^(3 + 1/(1 - rho))))/(sigma^(1/(-1 + rho))*(-4 + 3*rho)) +      (a3^2*(-1 + rho)^(1 + (2*(-3 + rho))/(-1 + rho))*(y^(3 - 4/(-1 + rho)) - y0^(3 - 4/(-1 + rho))))/(sigma^(4/(-1 + rho))*(2*(-7 + 3*rho))) +      (a2*a3*(-1 + rho)^(1 + (-3 + rho)/(-1 + rho) + (-2 + rho)/(-1 + rho))*sigma^(1/(1 - rho) - 2/(-1 + rho))*(y^(3 - 3/(-1 + rho)) - y0^(3 - 3/(-1 + rho))))/      (3*(-2 + rho)) + (a2^2*(-1 + rho)^(1 + (2*(-2 + rho))/(-1 + rho))*sigma^(2/(1 - rho))*(y^(3 - 2/(-1 + rho)) - y0^(3 - 2/(-1 + rho))))/(2*(-5 + 3*rho)) -      (a1*a3*(1 - rho)*(-1 + rho)^(1 + (-3 + rho)/(-1 + rho))*(y^(3 - 2/(-1 + rho)) - y0^(3 - 2/(-1 + rho))))/(sigma^(2/(-1 + rho))*(-5 + 3*rho)) +      (a2*am1*(-1 + rho)^(1 + (-2 + rho)/(-1 + rho) + (1 + rho)/(-1 + rho))*sigma^(1/(1 - rho) + 2/(-1 + rho))*(y^(3 + 1/(-1 + rho)) - y0^(3 + 1/(-1 + rho))))/      (-2 + 3*rho) - (a0*a1*(1 - rho)*(-1 + rho)^(1 + rho/(-1 + rho))*sigma^(1/(-1 + rho))*(y^(3 + 1/(-1 + rho)) - y0^(3 + 1/(-1 + rho))))/(-2 + 3*rho) +      
    (a0^2*(-1 + rho)^(1 + (2*rho)/(-1 + rho))*sigma^(2/(-1 + rho))*(y^(3 + 2/(-1 + rho)) - y0^(3 + 2/(-1 + rho))))/(-1 + 3*rho) -      (a1*am1*(1 - rho)*(-1 + rho)^(1 + (1 + rho)/(-1 + rho))*sigma^(2/(-1 + rho))*(y^(3 + 2/(-1 + rho)) - y0^(3 + 2/(-1 + rho))))/(-1 + 3*rho) -      (a0^2*(-1 + rho)^(1 + (2*rho)/(-1 + rho))*sigma^(2/(-1 + rho))*(y^(3 + 2/(-1 + rho)) - y0^(3 + 2/(-1 + rho))))/(-2 + 6*rho) +      (am1^2*(-1 + rho)^(1 + (2*(1 + rho))/(-1 + rho))*sigma^(4/(-1 + rho))*(y^(3 + 4/(-1 + rho)) - y0^(3 + 4/(-1 + rho))))/(1 + 3*rho) -      (am1^2*(-1 + rho)^(1 + (2*(1 + rho))/(-1 + rho))*sigma^(4/(-1 + rho))*(y^(3 + 4/(-1 + rho)) - y0^(3 + 4/(-1 + rho))))/(2 + 6*rho) -      (3*a3*(-1 + rho)^((-3 + rho)/(-1 + rho))*(y^((-3 + rho)/(-1 + rho)) - y0^((-3 + rho)/(-1 + rho))))/(sigma^(2/(-1 + rho))*(3 - rho)) -      (a3*(-1 + rho)^((-3 + rho)/(-1 + rho))*rho*(y^((-3 + rho)/(-1 + rho)) - y0^((-3 + rho)/(-1 + rho))))/(sigma^(2/(-1 + rho))*(6 - 2*rho)) -      (1/2)*a0*(-1 + rho)^(rho/(-1 + rho))*sigma^(1/(-1 + rho))*(y^(rho/(-1 + rho)) - y0^(rho/(-1 + rho))) -      (a0*(-1 + rho)^(1 + rho/(-1 + rho))*sigma^(1/(-1 + rho))*(y^(rho/(-1 + rho)) - y0^(rho/(-1 + rho))))/(-2 + 2*rho) -      
    (3*a3*(-1 + rho)^((-3 + rho)/(-1 + rho))*(y^(rho/(-1 + rho))*y0^(3/(-1 + rho)) - y^(3/(-1 + rho))*y0^(rho/(-1 + rho))))/(sigma^(2/(-1 + rho))*(y*y0)^(3/(-1 + rho))*(2*(-3 + rho))) - (a3*(-1 + rho)^((-3 + rho)/(-1 + rho))*rho*       (y^(rho/(-1 + rho))*y0^(3/(-1 + rho)) - y^(3/(-1 + rho))*y0^(rho/(-1 + rho))))/(sigma^(2/(-1 + rho))*(y*y0)^(3/(-1 + rho))*(-3 + rho)) -      (a3*(-1 + rho)^(1 + (-3 + rho)/(-1 + rho))*rho*(y^(rho/(-1 + rho))*y0^(3/(-1 + rho)) - y^(3/(-1 + rho))*y0^(rho/(-1 + rho))))/      (sigma^(2/(-1 + rho))*(y*y0)^(3/(-1 + rho))*((-3 + rho)*(-2 + 2*rho))) + (a0*am1*(-1 + rho)^(1 + rho/(-1 + rho) + (1 + rho)/(-1 + rho))*       sigma^(3/(-1 + rho))*(y^((3*rho)/(-1 + rho)) - y0^((3*rho)/(-1 + rho))))/(3*rho) -      (am1*(-1 + rho)^((1 + rho)/(-1 + rho))*sigma^(2/(-1 + rho))*(y^((1 + rho)/(-1 + rho)) - y0^((1 + rho)/(-1 + rho))))/(2*(1 + rho)) -      (am1*(-1 + rho)^((1 + rho)/(-1 + rho))*rho*sigma^(2/(-1 + rho))*(y^((1 + rho)/(-1 + rho)) - y0^((1 + rho)/(-1 + rho))))/(2*(1 + rho)) -      (am1*(-1 + rho)^(1 + (1 + rho)/(-1 + rho))*rho*sigma^(2/(-1 + rho))*(y^((1 + rho)/(-1 + rho)) - y0^((1 + rho)/(-1 + rho))))/((1 + rho)*(-2 + 2*rho)))  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
  output  <- -log(2*pi*del)/2 - log(s_x) + cm1/del + c0 + c1*del 
  
  print(output)
  return(output)
} 

# ModelU13(3,4,1/65,c(0.1,0.2,0.3,0.4,0.4,0.5,0.6))