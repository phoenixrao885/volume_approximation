library(volesti)
#library(NlcOptim)
#library(pracma)
library(nloptr)

d=20
k=40

P = GenZonotope(d,k)

G = t(P$G)
D = G%*%t(G)

D = (D+t(D))/2

zz = eigen(D, symmetric = TRUE)

Q = zz$vectors

a = Q[,1]

aa = a%*%G
solution = c()#vector("numeric", length = k)

for (j in 1:k) {
  if (aa[j]>0){
    solution = c(solution,1)
  } else {
    solution = c(solution,-1)
  }
}

print(solution)
#sol = as.matrix(solution)%*%G
DD = t(G)%*%G

DD = (DD+t(DD))/2
print(t(solution)%*%DD%*%solution)

