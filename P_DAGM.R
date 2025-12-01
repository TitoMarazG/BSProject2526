P_DAGM <- function(S, DAGM, a, U) { 
  #VERSIONE DI P_DAG CON MATRICE DI ADIACENZA: DAGM: 
  source("P_(XA).R")
  
  j=1:ncol(DAGM)
  
  sum_pa = 0
  sum_fa = 0
  
  for(i in 1:length(j)){
    
    pa <- which(DAGM[, i] == 1)
    fa <- c(pa,i)  
    
    if(any(pa != 0)){
      
      sum_pa <- sum_pa + (P_XA(S,a,U,pa))
      
    }
    
    sum_fa <- sum_fa + (P_XA(S,a,U,fa))
    
    #SERVE HANDLER PER IL CASO: pa= Vuoto
    #se pa=Vuoto la funzione parents restituisce un vettore char di lunghezza zero
    #check: if(lenght(pa)>0) => esegui prod_pa
    
    
    #ipotesi: se pa non è vuoto eseguo prod_pa, altrimenti skippo
    #equivale a dire P_XA=1 quando A=empty set (è giusto?)
    
    #LABEL PROBLEM: I DAG DI BNLEARN POSSONO AVERE LABELS LETTERE
    #OPZIONI PER RISOLVERE:
    #1) IMPORRE CHE I DAG VADANO DEFINITI CON LABEL NUMERICI INTERI (NO MODIFICHE IN PXA)
    #2) IMPORRE CORRISPONDENZA LETTERE <-> INTERI ("A" == "1") E COSTRUIRE HANDLER IN PXA
    # LA CORRSISPONDENZA LETTERE NUMERI DEVE RISPECCHIARSI NEL DATASET (A=1 = 1a col del dataset)
    
    
    
  }
  
  P_DAG = sum_fa - sum_pa
  
  return(P_DAG)
  
}
