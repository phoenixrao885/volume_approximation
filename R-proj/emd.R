emd <- function(F1, F2, W1, W2) {
  
  f = gdm(F1, F2)
  
  m = dim(F1)[1]
  n = dim(F2)[1]
  
  A1 = matrix(0, m, m * n)
  A2 = matrix(0, n, m * n)
  
  for (i in 1:m) {
    for (j in 1:n) {
      k = j + (i - 1) * n
      A1[i,k] = 1
      A2[j,k] = 1
    }
  }
  
  A = rbind(A1, A2)
  b = rbind(W1, W2)
  
  Aeq = matrix(1, m + n, m * n)
  beq = matrix(1, m + n, 1) * min(sum(W1), sum(W2))
  
  l = c(matrix(0, m * n))
  
  return(solve_lp(A,b,Aeq,beq,l,f))
  
}