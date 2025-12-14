gaussCItest -> function (x, y, S, suffStat) 
{
  z <- zStat(x, y, S, C = suffStat$C, n = suffStat$n)
  2 * pnorm(abs(z), lower.tail = FALSE)
  
  # calcola p = 2 * P(Z > |z|) (p-value del test Z bilaterale con Z normale std)
  # sotto H0 (indipendenza)  , la statistica di Fisher è una  N(0,1)
  # quindi quel p-value misura quanto è incompatibile il dato con l’indipendenza condizionata
  # in skeleton: if (p-value >= alpha) remove edge
 
}


zStat -> function (x, y, S, C, n) 
{
  r <- pcorOrder(x, y, S, C)
  
  # Calcola la correlazione parziale  usando la matrice C (correlazione o covarianza).
  # Se S è vuoto, pcorOrder ti restituisce di fatto la correlazione marginale.
  # Se S non è vuoto, usa la precision matrix (l’inversa della sottomatrice su {x,y,S}) per ottenere il parziale.
  
  res <- sqrt(n - length(S) - 3) * 0.5 * logQ1pm(r)
  
  # statistica Z di Fisher
  # sotto l’ipotesi di indipendenza condizionale e gaussianità, Z è una normale standard
  
  if (is.na(res)) 
    0
  else res
}
  
  