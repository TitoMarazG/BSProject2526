log.p_XA <- function(S, n, a, U, A){
  source("log.norm_cost.R")
  nonA <- setdiff(1:ncol(S), A)
  A.c <- length(A)  # cardinalità di A
  nonA.c <- length(nonA)  # cardinalità di nonA
  UAA <- U[A, A, drop = FALSE] #####
  SAA <- S[A, A, drop = FALSE]
  
  return(log(2*pi) * (-n*A.c/2)
         + log.norm_cost(a-nonA.c, UAA)
         - log.norm_cost(a+n-nonA.c, UAA+SAA))
}