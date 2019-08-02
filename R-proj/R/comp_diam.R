comp_diam <- function(P){
  
  d = dim(P$G)[2]
  k = dim(P$G)[1]
  
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
  
  #print(solution)
  
  DD = t(G)%*%G
  
  DD = (DD+t(DD))/2
  return(2*sqrt(abs(t(solution)%*%DD%*%solution)))
  
}