BF <- function(S, n, a, U, D0, D1, nodo_da_testare){
  source("log.P_DAGM_puntuale.R")
  return (exp(log.P_DAGM_puntuale(S,n, a, U, D0, nodo_da_testare) - log.P_DAGM_puntuale(S,n, a, U, D1, nodo_da_testare)))
}


# chiedere al prof: la lasciamo in log o facciamo l'exp?