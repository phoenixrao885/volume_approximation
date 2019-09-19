solve_lp <- function(A,b,Aeq,beq,l,f) {
  
  m = dim(A)[1]
  n = dim(A)[2]
  
  m2 = dim(Aeq)[1]
  
  AA = rbind(A,Aeq)
  bb = rbind(b,beq)
  
  lprec <- lpSolveAPI::make.lp(nrow(AA), ncol(AA))
  lpSolveAPI::set.objfn(lprec, f)
  
  lpSolveAPI::row.add.mode(lprec, "on")
  
  for (i in 1:m) {
    lpSolveAPI::add.constraint(lprec, A[i,], type = "<=", b[i])
  }
  
  for (i in 1:m2) {
    lpSolveAPI::add.constraint(lprec, Aeq[i,], type = "=", beq[i])
  }
  lpSolveAPI::row.add.mode(lprec, "off")
  
  lpSolveAPI::set.bounds(lprec, lower = l, upper = NULL)
  
  lpSolveAPI::solve(lprec)
  #print(sum(c(get.variables(lprec))))
  #x = get.variables(lprec)
  #num=get.solutioncount(lprec)
  
  #print(matrix(c(f),nrow = 1)%*%matrix(c(get.variables(lprec)),ncol = 1)/sum(c(get.variables(lprec))))
  
  val = lpSolveAPI::get.objective(lprec)/sum(c(lpSolveAPI::get.variables(lprec)))
  #delete.lp(lprec)
  
  return(val)
}