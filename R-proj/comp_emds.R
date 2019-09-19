s1 = get_sgn(sgns,1)
s2 = get_sgn(sgns,2)
F1 = s1$F
F2 = s2$F
W1 = s1$W
W2 = s2$W

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


s3 = get_sgn(sgns,3)
F3 = s3$F
W3 = s3$W
b2 = rbind(W1, W3)

D = matrix(rep(0,6616*6616), ncol = 6616, nrow = 6616)

for (i in 1:1) {
  for (j in (i+1):66) {
    s1 = get_sgn(sgns,i)
    s2 = get_sgn(sgns,j)
  
    print(paste0("i = ",i," j = ",j))
  
    b=rbind(s1$W, s2$W)
    D[i,j] = solve_lp22(A,b,Aeq,beq,l,f)
  }
  writeMat(paste0(getwd(),"/D_",i,"_",j,".mat"),D=D)
}


