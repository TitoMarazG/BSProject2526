# ==============================================================================
# FUNZIONI AUSILIARIE PER VALUTAZIONE PERFORMANCE
# ==============================================================================

# Estrae la matrice di adiacenza da un oggetto graph
# 
# g Oggetto di tipo graphNEL o matrice di adiacenza
# return Matrice di adiacenza binaria
extract_adjacency <- function(g) {
  if (is.matrix(g)) {
    adj <- g
  } else if (class(g) == "graphNEL") {
    adj <- as(g, "matrix")
  } else {
    stop("Input deve essere una matrice o un oggetto graphNEL")
  }
  
  # Assicurati che sia binaria
  adj[adj != 0] <- 1
  
  return(adj)
}

# Calcola metriche di confusione per grafi
# 
# true_graph Grafo vero (graphNEL o matrice)
# estimated_graph Grafo stimato (graphNEL o matrice)
# undirected Se TRUE, considera i grafi come non orientati
# return Lista con TP, FP, FN, TN
confusion_metrics <- function(true_graph, estimated_graph, undirected = TRUE) {
  
  # Estrai matrici di adiacenza
  true_adj <- extract_adjacency(true_graph)
  est_adj <- extract_adjacency(estimated_graph)
  
  # Per grafi non orientati, considera solo il triangolo inferiore
  if (undirected) {
    true_adj[upper.tri(true_adj, diag = TRUE)] <- 0
    est_adj[upper.tri(est_adj, diag = TRUE)] <- 0
  }
  
  # Calcola metriche di confusione
  TP <- sum(true_adj == 1 & est_adj == 1)  # Veri positivi
  FP <- sum(true_adj == 0 & est_adj == 1)  # Falsi positivi
  FN <- sum(true_adj == 1 & est_adj == 0)  # Falsi negativi
  TN <- sum(true_adj == 0 & est_adj == 0)  # Veri negativi
  
  return(list(TP = TP, FP = FP, FN = FN, TN = TN))
}

# Calcola precision per un grafo stimato
# 
# true_graph Grafo vero (graphNEL o matrice)
# estimated_graph Grafo stimato (graphNEL o matrice)
# undirected Se TRUE, considera i grafi come non orientati
# return Valore di precision (0-1)
precision <- function(true_graph, estimated_graph, undirected = TRUE) {
  
  metrics <- confusion_metrics(true_graph, estimated_graph, undirected)
  
  # Precision = TP / (TP + FP)
  if (metrics$TP + metrics$FP == 0) {
    return(NA)  # Non ci sono archi stimati
  }
  
  prec <- metrics$TP / (metrics$TP + metrics$FP)
  return(prec)
}

# Calcola recall (sensitivity) per un grafo stimato
# 
# true_graph Grafo vero (graphNEL o matrice)
# estimated_graph Grafo stimato (graphNEL o matrice)
# undirected Se TRUE, considera i grafi come non orientati
# return Valore di recall (0-1)
recall <- function(true_graph, estimated_graph, undirected = TRUE) {
  
  metrics <- confusion_metrics(true_graph, estimated_graph, undirected)
  
  # Recall = TP / (TP + FN)
  if (metrics$TP + metrics$FN == 0) {
    return(NA)  # Non ci sono archi veri
  }
  
  rec <- metrics$TP / (metrics$TP + metrics$FN)
  return(rec)
}

# Calcola F1-score per un grafo stimato
# 
# true_graph Grafo vero (graphNEL o matrice)
# estimated_graph Grafo stimato (graphNEL o matrice)
# undirected Se TRUE, considera i grafi come non orientati
# return Valore di F1-score (0-1)
f1_score <- function(true_graph, estimated_graph, undirected = TRUE) {
  
  prec <- precision(true_graph, estimated_graph, undirected)
  rec <- recall(true_graph, estimated_graph, undirected)
  
  if (is.na(prec) || is.na(rec) || (prec + rec == 0)) {
    return(NA)
  }
  
  f1 <- 2 * (prec * rec) / (prec + rec)
  return(f1)
}

# Calcola tutte le metriche di valutazione
# 
# true_graph Grafo vero (graphNEL o matrice)
# estimated_graph Grafo stimato (graphNEL o matrice)
# undirected Se TRUE, considera i grafi come non orientati
# return Data frame con tutte le metriche
evaluate_graph <- function(true_graph, estimated_graph, undirected = TRUE) {
  
  metrics <- confusion_metrics(true_graph, estimated_graph, undirected)
  prec <- precision(true_graph, estimated_graph, undirected)
  rec <- recall(true_graph, estimated_graph, undirected)
  f1 <- f1_score(true_graph, estimated_graph, undirected)
  shd_val <- shd(true_graph, estimated_graph)
  
  results <- data.frame(
    SHD = shd_val,
    Precision = round(prec, 4),
    Recall = round(rec, 4),
    F1_Score = round(f1, 4),
    TP = metrics$TP,
    FP = metrics$FP,
    FN = metrics$FN,
    TN = metrics$TN
  )
  
  return(results)
}