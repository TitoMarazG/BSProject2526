BF <- function(S, DAGM0, DAGM1, a, U, nodo_da_testare){
  source("P_DAGM_puntuale.R")
  return(P_DAGM_puntuale(S, DAGM0, a, U, nodo_da_testare)/P_DAGM_puntuale(S, DAGM1, a, U, nodo_da_testare))
}
