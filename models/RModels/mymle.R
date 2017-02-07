library('nloptr')
library('numDeriv')


mymle<-function(logdensity,x,del,param0){
  
  objfun <- function(param){
    # Set the objective function
    objfun <- -logdensity2loglik(logdensity,x,del,param)
  }  
  
  # Optimize
  # nloptr.print.options()  
  
  # res <- optim()
  res<- nloptr( x0=param0,
                eval_f=objfun, #f(x)
                eval_g_ineq=eval_g_ineq, #g(x)
                lb = args$l,
                ub = args$u,
                opts = list("algorithm"="NLOPT_LN_COBYLA", "maxeval" = args$maxiter, "xtol_rel" = args$eps, "print_level"=args$print_level))
  # print(res)
  
  
  H <- hessian(objfun,param_0)
  InfoMatrix <- logdensity2info(logdensity,x,del,param_0)
  Variance <- solve(InfoMatrix)
  invH <- solve(H)
  Variance_Robust <- invH * InfoMatrix * invH
  print('Estimation finished.\n\n')
  print(InfoMatrix)
  
  return(res)
  return(InfoMatrix)
  
}  
################################################################################################################  

logdensity2info <- function(logdensity,x,del,param){
  # Compute the information matrix = Sum_{i=1}^n Score_i' Score_i
  n <- length(x) - 1
  output <- 0
  for (i in 1:n)
  {
    tmpfun <- (theta)(logdensity(x(i+1),x(i),del,theta))
    score_i <- gradest(tmpfun,param) # This is a row vector
    output <- output + score_i * score_i
  }
}

function [grad,err,finaldelta] = gradest(fun,x0)
sx = size(x0);

% total number of derivatives we will need to take
nx = numel(x0);

grad = zeros(1,nx);
err = grad;
finaldelta = grad;
for ind = 1:nx
  [grad(ind),err(ind),finaldelta(ind)] = derivest( ...
    @(xi) fun(swapelement(x0,ind,xi)), ...
    x0(ind),'deriv',1,'vectorized','no', ...
    'methodorder',2);
end

end

function vec = swapelement(vec,ind,val)
% swaps val as element ind, into the vector vec
vec(ind) = val;

end % sub-function end



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




# % Set default options for the optimization procedure
# if nargin == 4
# optim_options = optimset('Display','iter','TolX',1e-2,'TolFun',1e-3);
# elseif nargin == 5
# optim_options = varargin{1};
# else
#   error('The number of input is wrong');
# end
# 
# 
# % Set the objective function
# objfun = @(param)(-logdensity2loglik(logdensity,x,del,param));
# 
# % Optimize
# [param_mle,fval,exitflag] = fminsearch(objfun,param0,optim_options);
# 
# % Print
# if exitflag == 1
# output.exitmsg = 'fminsearch converged to a solution.';
# elseif exitflag == 0
# output.exitmsg = 'Maximum number of function evaluations or iterations was reached.';
# elseif exitflag == -1
# output.exitmsg = 'Algorithm was terminated by the output function.';
# end