library(volesti)
library(NlcOptim)

d=10
k=20

P = GenZonotope(d,k)
HP = GenCube(k, 'H')

ub = rep(1, k)
lb = rep(-1, k)
x0 = sample_points(HP, N=1, WalkType = "RDHR", walk_step = 4*d)
#x0 = rep(0, 6)

#A = HP$A
#b = H$b
G = t(P$G)
G = t(G)%*%G

#x <- vector("numeric", length = k)

relvar.rev <- function(x){
  rv = rep(0,k)
  for (i in 1:k) {
    for (j in 1:k) {
      rv[i] = rv[i]+x[j]*G[j,i]
    }
  }
  rv2 = 0
  for (i in 1:k) {
    rv2 = rv2 + rv[i]*x[i]
  }
  return(-rv2)
}

retList = solnl(x0,objfun=relvar.rev,lb=lb, ub=ub)
