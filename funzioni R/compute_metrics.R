# Funzione helper per calcolare True Positives (TP), False Positives (FP), SHD
compute_metrics <- function(estimated_skel, true_dag) {
  # Convertiamo tutto in matrici di adiacenza simmetriche (scheletro)
  true_amat <- as(true_dag, "matrix")
  true_amat[true_amat != 0] <- 1
  true_amat <- pmax(true_amat, t(true_amat)) # Simmetrizza
  
  est_amat <- as(estimated_skel@graph, "matrix")
  est_amat[est_amat != 0] <- 1
  est_amat <- pmax(est_amat, t(est_amat)) # Simmetrizza (anche se dovrebbe già esserlo)
  
  # Calcolo SHD (Structural Hamming Distance) sugli scheletri
  # SHD = FP + FN (sugli scheletri non contiamo errori di direzione)
  diff <- abs(est_amat - true_amat)
  shd_val <- sum(diff) / 2 # Diviso 2 perché la matrice è simmetrica
  
  # Precision e Recall (sugli archi)
  # TP: Arco stimato presente anche nel vero
  TP <- sum((est_amat == 1) & (true_amat == 1)) / 2
  # FP: Arco stimato ma non presente nel vero
  FP <- sum((est_amat == 1) & (true_amat == 0)) / 2
  # FN: Arco vero non trovato
  FN <- sum((est_amat == 0) & (true_amat == 1)) / 2
  
  precision <- if((TP+FP)>0) TP / (TP + FP) else 0  # TRUE / TRUE estimated
  recall    <- if((TP+FN)>0) TP / (TP + FN) else 0  # TRUE / TRUE real
  
  return(c(SHD = shd_val, Precision = precision, Recall = recall))
}

# # ==============================================================================
# # 7. ANALISI DEI BAYES FACTOR (Opzionale)
# # ==============================================================================
# # Esempio: vediamo la storia dei BF per una coppia specifica se esiste
# if (!is.null(res_bayes$BF_history)) {
#   cat("\nNota: I valori di Bayes Factor sono salvati in res_bayes$BF_history.\n")
#   cat("Puoi ispezionarli per capire 'quanto' forte era l'indipendenza.\n")
# }