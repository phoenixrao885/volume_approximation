library(volesti)
#library(NlcOptim)
#library(pracma)
library(nloptr)

d=100
k=200

P = GenZonotope(d,k)
HP = GenCube(k, 'H')

ub = rep(1, k)
lb = rep(-1, k)
#x0 = sample_points(HP, N=1, WalkType = "RDHR", walk_step = 10*d)
x0 = runif(k,min=-1, max=1)
#for (j in 1:k) {
#  if (runif(1)>0.5){
#    x0[j]=1
#  } else {
#    x0[j]=-1
#  }
#}
print(x0)
#x0 = rep(0, 6)

#A = HP$A
#b = H$b
G = t(P$G)
G = t(G)%*%G

G = (G+t(G))/2

#x <- vector("numeric", length = k)

relvar <- function(x){
  #rv = rep(0,k)
  #for (i in 1:k) {
    #for (j in 1:k) {
      #rv[i] = rv[i]+x[j]*G[j,i]
    #}
  #}
  #rv2 = 0
  #for (i in 1:k) {
    #rv2 = rv2 + rv[i]*x[i]
  #}
  #rv2 = -rv2
  return(-t(x)%*%G%*%x)
}

eval_grad_f <- function(x){
  return(-2*(G%*%x))
}
x1=x0
retList$solution = 10*x0
#retList = solnl(x0=x0,objfun=relvar.rev,lb=lb, ub=ub)NLOPT_LD_SLSQP
while (sqrt(sum((x0-retList$solution)^2))>0.001){
  print(sqrt(sum((x0-retList$solution)^2)))
  x0 = x1
  retList = nloptr(x0=x0,eval_f=relvar,lb=lb, ub=ub, opts = list("algorithm"="NLOPT_GN_DIRECT"))
  x1 = retList$solution
  #retList = fminsearch(x0=x0,fn=relvar.rev,lower=lb, upper=ub,method = "Hooke-Jeeves")NLOPT_GN_CRS2_LM
}
print(x1)
retList = nloptr(x0=x1,eval_f=relvar,eval_grad_f=eval_grad_f,lb=lb, ub=ub, opts = list("algorithm"="NLOPT_LD_SLSQP"))
x1 = retList$solution
print(x1)
print(x1%*%G%*%x1)
