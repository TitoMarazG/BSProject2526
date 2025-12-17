log.p_XA <- function(XX, n, a, U, A){
  source("log.norm_cost.R")
  nonA <- setdiff(1:ncol(XX), A)
  A.c <- length(A)  # cardinalità di A
  nonA.c <- length(nonA)  # cardinalità di nonA
  UAA <- U[A, A, drop = FALSE] 
  XXAA <- XX[A, A, drop = FALSE]
  
  return(log(2*pi) * (-n*A.c/2)
         + log.norm_cost(a-nonA.c, UAA)
         - log.norm_cost(a+n-nonA.c, UAA+XXAA))
}

#XX = X(T) * X 
#a,U hyperparemters
#A: set di variabili selezionate (parents/family)
#n sample size
#calcola P di un subset A di variabili