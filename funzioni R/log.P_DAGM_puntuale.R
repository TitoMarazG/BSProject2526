log.P_DAGM_puntuale <- function(XX, n, a, U, v, S) { 
  #VERSIONE DI P_DAG_puntuale CON MATRICE DI ADIACENZA: DAGM: 
  source("log.p_XA.R")
  return(log.p_XA(XX,n,a,U,c(S, v)) - log.p_XA(XX,n,a,U,S))
}

#S:parents di entrambi
#v:nodo da testare