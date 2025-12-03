P_DAGM_puntuale <- function(S, n, DAGM, a, U, nodo_da_testare) { 
  #VERSIONE DI P_DAG_puntuale CON MATRICE DI ADIACENZA: DAGM: 
  source("log.p_XA.R")
  
  pa <- which(DAGM[, nodo_da_testare] == 1)
  fa <- c(pa,i)
  if(any(pa != 0)){
  return(log.p_XA(S,n,a,U,fa) - log.p_XA(S,n,a,U,pa))
  }
  else {
    return(log.p_XA(S,n,a,U,fa))
  }
}