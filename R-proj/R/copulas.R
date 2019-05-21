#' @export
copulas <- function(MatReturns, Win = NULL, M = NULL, N = NULL, error = NULL, prob = NULL) {
  
  if (is.null(W)) Win =60
  nrows = dim(MatReturns)[1]
  wl = Win-1
  E = cov_shrink(MatReturns[1:Win,])
  
  for (i in 1:(nrows-wl)) {
    W=i:(i+wl)
    E = rbind(E, cov_shrink(MatReturns[W,]))
  }
  
  return( get_copulas(MatReturns, E, Win, M, N, error, prob) )
  
}