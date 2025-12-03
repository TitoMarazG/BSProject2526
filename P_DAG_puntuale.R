P_DAG_puntuale <- function(X, DAG, a, U, nodo_da_testare) {
  
  #install.packages("bnlearn")
  library("bnlearn")
  source(Function_p(XA))
  
  pa <- parents(DAG,nodo_da_testare)
  fa <- c(pa , nodo_da_testare)
  
  P_DAG = P_XA(X,a,U,fa)/P_XA(X,a,U,pa)
  return(P_DAG)
}
  