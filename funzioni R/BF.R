BF <- function(S, n, a, U, D1, D0){
  source("log.P_DAGM_puntuale.R")

  for(i in 1:ncol(D1)){
    for (j in ncol(D1)) {
      if(D0[i,j]!=D1[i,j]){
        nodo_da_testare = j
      }
    }
  }
  
  return (exp(log.P_DAGM_puntuale(S,n, a, U, D0, nodo_da_testare) - log.P_DAGM_puntuale(S,n, a, U, D1, nodo_da_testare)))
    
  }
 
}


# chiedere al prof: la lasciamo in log o facciamo l'exp?
## chiedere se la matrice di adiacenza è meglio che sia completa o triangolare, diagonale di zeri o uni?