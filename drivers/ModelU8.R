

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
  sig      <- param[5]  # sigma
  r        <- param[6]  # rho
  
  
  param_prime <- c(a,b,g, h, sig, r)
  output <- output + ModelU8(x,x0,del,param_prime)
  
  return(output)
}


ModelU8 <- function(x,x0,del,param)
  
  # (* DO NOT LET THE PARAMETER rho GO BELOW 1! THE EXPANSION IS DIFFERENT THERE. 
  # 
  # If you use this code, please acknowledge the source:Ait-Sahalia,Yacine,1999, Transition Densities for Interest Rate and Other Nonlinear Densities, Journal of Finance 54,1361-1395.
  # Ait-Sahalia,Yacine,2002, Maximum-Likelihood Estimation of Discretely Sampled Diffusions: A Closed-Form Approximation Approach, Econometrica 70,223-262. *)
  
{
  am1 <- param[1]
  a0 <- param[2] 
  a1 <- param[3] 
  a2 <- param[4] 
  sigma <- param[5] 
  rho <- param[6] 
  
  if (rho<1)
    # If rho<1, set the log density to -Infinity for estimation purpose.
  {
    output <- -Inf 
  }
  
  
  y <- x^(1 - rho)/(sigma*(rho - 1)) 
  y0 <- x0^(1 - rho)/(sigma*(rho - 1)) 
  sx <- sigma*x^rho 
  
  cYm1 <- -(1/2)*(y - y0)^2 
  cY0 <- -(1/2) * a1 * (-1 + rho) * (y - y0) * (y + y0) - ( a2 * (-1 + rho)^((-3 + 2 * rho)/(-1 + rho)) * 
                                                           sigma^(1/(  1 - rho)) * (y^(2 + 1/(1 - rho)) - y0^(2 + 1/(1 - rho))))/(-3 + 2 * rho) - ( a0 * (-1 + rho)^(1 + rho/(-1 + rho)) * 
                                                                                                                                                    sigma^(  1/(-1 + rho)) * (y^(2 + 1/(-1 + rho)) - y0^(2 + 1/(-1 + rho))))/(-1 +   2 * rho) 
  - ( am1 * (-1 + rho)^((2 * rho)/(-1 + rho))* sigma^(  2/(-1 + rho)) * (y^((2 * rho)/(-1 + rho)) - y0^((2 * rho)/(-1 + rho))))/( 2 * rho) 
  + (rho * log(y/y0))/(-2 + 2 * rho)  
  cY1 <- -(1/(y - y0)) * (( am1^2 * (-1 + rho)^((1 + 3 * rho)/(-1 + rho)) * sigma^( 4/(-1 + rho)) * 
                           (y^(3 + 4/(-1 + rho)) - y0^(3 + 4/(-1 + rho))))/( 2 + 6 * rho) + 1/6 * (a1^2 * (-1 + rho)^2 * (y^3 - y0^3) +  
                                                                                                   ( 3 * a0^2 * (-1 + rho)^(1 + (2 * rho)/(-1 + rho)) * sigma^( 2/(-1 + rho)) * (y^(3 + 2/(-1 + rho)) - y0^(3 + 2/(-1 + rho))))/(-1 + 3 * rho)  
                                                                                                   + ( 3 * a1 * (1 - rho) * (((4 - 11 * rho + 6 * rho^2) * (y - y0))/(-1 + rho) - 2 * a2 * (-1 + rho)^((-3 + 2 * rho)/(-1 + rho)) * 
                                                                                                                             sigma^(1/( 1 - rho)) * (y^(3 + 1/(1 - rho)) - y0^(3 + 1/(1 - rho)))))/(-4 + 3 * rho)  
                                                                                                   + 3 * (-(((-2 + rho) * rho * (y - y0))/(4 * (-1 + rho)^2 * y * y0))  
                                                                                                          + ( a2^2 * (-1 + rho)^((-5 + 3 * rho)/(-1 + rho)) * sigma^(-( 2/(-1 + rho))) * (y^(3 - 2/(-1 + rho)) - y0^( 3 - 2/(-1 + rho))))/(-5 + 3 * rho) 
                                                                                                          - ( 2 * a2 * (-1 + rho)^((-3 + 2 * rho)/(-1 + rho)) * sigma^(1/( 1 - rho)) * (y^((-2 + rho)/(-1 + rho)) 
                                                                                                                                                                                        - y0^((-2 + rho)/(-1 + rho))))/(-2 + rho)) - (1/(rho * (-2 + 3 * rho))) * a0 * (-1 + rho)^(rho/(-1 + rho)) * sigma^( 1/(-1 + rho)) * 
                                                                                                   (-6 * a1 * (-1 + rho)^2 * rho * (y^(3 + 1/(-1 + rho)) - y0^( 3 + 1/(-1 + rho))) + (-2 + 3 * rho) * (-2 * a2 *(-1 + rho)^((-2 + rho)/(-1 + rho)) * 
                                                                                                                                                                                                       rho * sigma^(1/( 1 - rho)) * (y^3 - y0^3) + 6 * rho * (y^(rho/(-1 + rho)) - y0^(rho/(-1 + rho))))))  
                         - (1/( 6 * rho * (2 - 7 * rho + 9 * rho^3))) * am1 * (-1 + rho)^((1 + rho)/(-1 + rho)) * sigma^( 2/(-1 + rho)) * 
                         (-2 * a0 * (-1 + rho)^( rho/(-1 + rho)) * (-2 + 9 * rho - 7 * rho^2 - 9 * rho^3 + 9 * rho^4) * sigma^( 1/(-1 + rho)) * 
                           (y^((3 * rho)/(-1 + rho)) - y0^((3 * rho)/(-1 + rho))) + 3 * rho * (-2 * a1 * (-1 + rho)^2 * (-2 + rho + 3 * rho^2) * (y^(3 + 2/(-1 + rho)) 
                                                                                                                                                  - y0^(3 + 2/(-1 + rho))) + (-1 + 3 * rho) * (-2 * a2 * (-1 + rho)^((-2 + rho)/(-1 + rho)) * (-1 + rho^2)  
                                                                                                                                                                                               * sigma^(1/( 1 - rho)) * (y^(3 + 1/(-1 + rho)) - y0^(3 + 1/(-1 + rho))) + (1 + 2 * rho) * (-2 + 3 * rho) * (y^((1 + rho)/(-1 + rho)) 
                                                                                                                                                                                                                                                                                                           - y0^(( 1 + rho)/(-1 + rho))))))) 
  output <- (-(1/2))*log(2*pi*del) - log(sx) + cYm1/del + cY0 + cY1*del 
  
  print(output)
  return(output)
}

# ModelU8(4,5,1/32,c(0.1,0.2,0.3,0.4,0.5,4))