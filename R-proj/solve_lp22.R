solve_lp22 <-function(A,b,Aeq,beq,l,f) {
  
  m = dim(A)[1]
  n = dim(A)[2]
  
  m2 = dim(Aeq)[1]
  
  lprec = lpSolve::lp(direction = "min", objective.in = f, const.mat = rbind(A,Aeq), 
                      const.dir = c(rep("<=",m), rep("=",m2)), const.rhs = c(rbind(b,beq)))
  
  return(lprec$objval)
  
}