

diam_hpoly_zono_inter <- function(G,A,b){
  
  d = dim(G)[1]
  k = dim(G)[2]
  
  #G = t(P$G)
  
  ub = rep(1, k)
  lb = rep(-1, k)
  
  DD = t(G)%*%G
  DD = (DD+t(DD))/2
  
  x0 = runif(k)
  
  relvar <- function(x){
    return(-t(x)%*%DD%*%x)
  }
  
  eval_grad_f <- function(x){
    return(-2*(DD%*%x))
  }
  
  eval_g_ineq <- function(x) {
    return(A%*%x-b)
  }
  
  eval_jac_g_ineq <- function(x) {
    return(A)
  }
  
  x1=x0
  #retList=list()
  #retList$solution = 10*x0
  #retList = solnl(x0=x0,objfun=relvar.rev,lb=lb, ub=ub)NLOPT_LD_SLSQP
  #while (sqrt(sum((x0-retList$solution)^2))>0.001){
  #print(sqrt(sum((x0-retList$solution)^2)))
  #x0 = x1
  #retList = nloptr(x0=x0,eval_f=relvar,lb=lb, ub=ub, opts = list("algorithm"="NLOPT_GN_DIRECT"),eval_g_eq=eval_g_eq)
  #x1 = retList$solution
  #retList = fminsearch(x0=x0,fn=relvar.rev,lower=lb, upper=ub,method = "Hooke-Jeeves")NLOPT_GN_CRS2_LM
  #}
  #print(x1)
  retList = nloptr::nloptr(x0=x1, eval_f = relvar, eval_grad_f = eval_grad_f, lb = lb, ub = ub,
                           opts = list("algorithm"="NLOPT_LD_SLSQP"),eval_g_ineq = eval_g_ineq,
                           eval_jac_g_ineq = eval_jac_g_ineq)
  x1 = retList$solution
  #print(sqrt(sum((G%*%x1)^2)))
  return(2*sqrt(sum((G%*%x1)^2)))
}


