log.P_DAGM <- function(S,n, DAGM, a, U) { 
  #VERSIONE DI P_DAG CON MATRICE DI ADIACENZA: DAGM: 
  source("log.p_XA.R")
  
  j=1:ncol(DAGM)
  
  sum_pa = 0
  sum_fa = 0
  
  for(i in 1:length(j)){
    
    pa <- which(DAGM[, i] == 1)
    fa <- c(pa,i)  
    
    if(any(pa != 0)){
      
      sum_pa <- sum_pa + (log.p_XA(S,n,a,U,pa))
      
    }
    
    sum_fa <- sum_fa + (log.p_XA(S,n,a,U,fa))
   
  }
  
  P_DAG = sum_fa - sum_pa
  
  return(P_DAG)
}
