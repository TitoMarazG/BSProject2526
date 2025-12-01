P_DAGM_puntuale <- function(S, DAGM, a, U, nodo_da_testare) { 
  #VERSIONE DI P_DAG_puntuale CON MATRICE DI ADIACENZA: DAGM: 
  source(Function_p(XA))
  
  pa <- which(DAGM[, nodo_da_testare] == 1)
  fa <- c(pa,i)
  
  return(P_XA(S,a,U,fa)/P_XA(S,a,U,pa))
}