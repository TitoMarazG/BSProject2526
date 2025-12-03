BF <- function(S, n, D0, D1, a, U, nodo_da_testare){
  source("P_DAGM_puntuale.R")
  return (exp(P_DAGM_puntuale(S,n, D0, a, U, nodo_da_testare) - P_DAGM_puntuale(S,n, D1, a, U, nodo_da_testare)))
}
