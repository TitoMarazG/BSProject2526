rm(list=ls()); graphics.off(); cat("\014")

library(hdi)
library(pcalg)
library(Rgraphviz) 
library(methods)

source("/Users/leomarcellopoli/Documents/Bayesian/Project/BSProject2526/funzioni R/personalized_my_skeleton_senza_stable_fast.R")

# DATASET LOADING & PREPARATION ####
data(riboflavin)
X_all <- as.matrix(riboflavin$x) 
y <- as.numeric(riboflavin$y)    

gene_variances <- apply(X_all, 2, var)

# Bühlmann usa q=100.
q_selected <- 1000  # nota: il frequentista potrebbe dare problemi

top_indices <- order(gene_variances, decreasing = TRUE)[1:q_selected]
X_reduced <- X_all[, top_indices]

cat(paste("Selezionati i top", q_selected, "geni con maggiore varianza.\n"))

# Dataset finale
data_final <- cbind(X_reduced, y)
colnames(data_final)[ncol(data_final)] <- "logRiboflavin" # La variabile target

# Standardizzazione
data_final_scaled <- scale(data_final)

n <- nrow(data_final_scaled)       # 71
q_total <- ncol(data_final_scaled) # q_selected + 1

cat("Dimensioni finali dataset: n =", n, "| q =", q_total, "\n")


# METODO BAYESIANO ####
cat("\n--- AVVIO METODO BAYESIANO ---\n")

# Parametri
a_prior <- q_total + 1  
U_prior <- diag(1, q_total)
p_posterior <- 0.9 # Conservative
alpha_bayes <- 1/p_posterior - 1 

time_start_B <- Sys.time()

res_bayes <- tryCatch({
  my_skeleton(X = data_final_scaled, 
              a = a_prior, 
              U = U_prior, 
              alpha = alpha_bayes, 
              BayesTest = BF_Gaussian, 
              saveBF = FALSE, 
              verbose = FALSE) 
}, error = function(e) {
  print(paste("ERRORE BAYES:", e$message))
  return(NULL)
})

time_end_B <- Sys.time()

# ==============================================================================
# 3. METODO FREQUENTISTA (Benchmark: SOLO SKELETON)
# ==============================================================================
cat("\n--- AVVIO METODO FREQUENTISTA (Skeleton) ---\n")

# Bühlmann 2014 usa alpha = 0.01
alpha_freq <- 0.01 

# Calcolo sufficient statistics
suffStat <- list(C = cor(data_final_scaled), n = n)

time_start_F <- Sys.time()

# Recuperiamo i nomi corretti delle variabili
node_labels <- colnames(data_final_scaled)

res_freq <- tryCatch({
  # USARE SKELETON INVECE DI PC
  pcalg::skeleton(suffStat = suffStat,
                  indepTest = gaussCItest,
                  labels = node_labels,    
                  alpha = alpha_freq,
                  verbose = FALSE)
}, error = function(e) {
  print(paste("ERRORE FREQ:", e$message))
  return(NULL)
})

time_end_F <- Sys.time()


# ==============================================================================
# 4. CONFRONTO E ANALISI RISULTATI (SKELETON vs SKELETON)
# ==============================================================================

if (!is.null(res_bayes) && !is.null(res_freq)) {
  
  # --- FASE 1: CALCOLI E PREPARAZIONE DATI ---
  
  # Estrazione grafi
  g_bayes <- res_bayes@graph
  g_freq  <- res_freq@graph
  
  # Metriche Globali
  edges_B <- numEdges(g_bayes)
  edges_F <- numEdges(g_freq)
  
  # --- CALCOLO SHD (Distanza Strutturale tra Skeleton) ---
  shd_val <- pcalg::shd(g_bayes, g_freq)
  
  # Estrazione Matrici di Adiacenza
  # Essendo Skeleton, sono matrici simmetriche (1 se c'è arco, 0 altrimenti)
  adj_B <- as(g_bayes, "matrix")
  adj_F <- as(g_freq,  "matrix")
  
  # Ordiniamo le matrici per essere sicuri al 100% che siano allineate
  common_nodes <- sort(colnames(data_final_scaled))
  adj_B <- adj_B[common_nodes, common_nodes]
  adj_F <- adj_F[common_nodes, common_nodes]
  
  # Calcolo Jaccard
  intersect_edges <- sum((adj_B == 1) & (adj_F == 1)) / 2
  union_edges     <- sum((adj_B == 1) | (adj_F == 1)) / 2
  jaccard_index   <- if(union_edges > 0) intersect_edges / union_edges else 0
  
  # Analisi Hubs (RowSums va benissimo su matrici simmetriche)
  deg_B <- rowSums(adj_B)
  deg_F <- rowSums(adj_F)
  
  max_deg_B <- if(length(deg_B)>0) max(deg_B) else 0
  max_gene_B <- if(length(deg_B)>0) names(which.max(deg_B)) else "N/A"
  
  max_deg_F <- if(length(deg_F)>0) max(deg_F) else 0
  max_gene_F <- if(length(deg_F)>0) names(which.max(deg_F)) else "N/A"
  
  # Vicini di logRiboflavin (Target)
  # Essendo skeleton simmetrici, basta guardare la riga del target
  target_node <- "logRiboflavin"
  
  neigh_B <- if(target_node %in% rownames(adj_B)) {
    names(which(adj_B[target_node, ] == 1))
  } else { character(0) }
  
  neigh_F <- if(target_node %in% rownames(adj_F)) {
    names(which(adj_F[target_node, ] == 1))
  } else { character(0) }
  
  common_genes <- intersect(neigh_B, neigh_F)
  
  # Check Geni Noti (Bühlmann)
  buhlmann_genes <- c("YXLD_at", "LYSC_at", "YOAB_at", "YCIC_at")
  selected_buhlmann <- intersect(buhlmann_genes, colnames(data_final_scaled))
  
  
  # --- FASE 2: SALVATAGGIO FILE E PLOTTING SEPARATO ---
  
  saveRDS(res_bayes, paste0("Riboflavin_Bayes_q", q_selected, ".rds"))
  saveRDS(res_freq,  paste0("Riboflavin_Freq_q",  q_selected, ".rds"))
  
  # Definizione funzione colori
  get_colors <- function(g, target) {
    node_names <- nodes(g)
    fill_col <- rep("white", length(node_names)); names(fill_col) <- node_names
    shapes   <- rep("ellipse", length(node_names)); names(shapes) <- node_names
    
    if(target %in% node_names) {
      fill_col[target] <- "tomato"
      shapes[target]   <- "box"
    }
    return(list(fillcolor = fill_col, shape = shapes))
  }
  
  # 1. PLOT BAYESIANO
  pdf(paste0("Riboflavin_Bayes_Plot_q", q_selected, ".pdf"), width = 12, height = 12)
  if (edges_B > 0 && edges_B < 2000) { 
    attrs_B <- get_colors(g_bayes, "logRiboflavin")
    plot(g_bayes, nodeAttrs = attrs_B, main = paste("Bayes Skeleton (Edges:", edges_B, ")"))
  } else {
    plot(new("graphNEL", nodes=nodes(g_bayes)), main="Bayes (Too dense or Empty)")
  }
  dev.off()
  
  # 2. PLOT FREQUENTISTA
  pdf(paste0("Riboflavin_Freq_Plot_q", q_selected, ".pdf"), width = 12, height = 12)
  if (edges_F > 0 && edges_F < 2000) {
    attrs_F <- get_colors(g_freq, "logRiboflavin")
    plot(g_freq, nodeAttrs = attrs_F, main = paste("Freq Skeleton (Edges:", edges_F, ")"))
  } else {
    plot(new("graphNEL", nodes=nodes(g_freq)), main="Freq (Too dense or Empty)")
  }
  dev.off() 
  
  
  # --- FASE 3: STAMPA FINALE ---
  
  # Calcolo durate
  duration_B <- difftime(time_end_B, time_start_B, units = "mins")
  duration_F <- difftime(time_end_F, time_start_F, units = "mins")
  
  cat("\n\n")
  cat("==========================================================\n")
  cat("               REPORT FINALE DI ANALISI                   \n")
  cat("==========================================================\n")
  
  cat("\n--- 0. TEMPI DI ESECUZIONE ---\n")
  cat(sprintf("Metodo Bayesiano:      %.2f minuti\n", duration_B))
  cat(sprintf("Metodo Frequentista:   %.2f minuti\n", duration_F))
  
  cat("\n--- 1. FILE SALVATI ---\n")
  cat("RDS Bayes:     ", paste0("Riboflavin_Bayes_q", q_selected, ".rds"), "\n")
  cat("RDS Freq:      ", paste0("Riboflavin_Freq_q",  q_selected, ".rds"), "\n")
  cat("PDF Bayes:     ", paste0("Riboflavin_Bayes_Plot_q", q_selected, ".pdf"), "\n")
  cat("PDF Freq:      ", paste0("Riboflavin_Freq_Plot_q",  q_selected, ".pdf"), "\n")
  
  cat("\n--- 2. CONFRONTO SKELETON ---\n")
  cat(sprintf("Archi totali Bayesiano:    %.0f\n", edges_B))
  cat(sprintf("Archi totali Frequentista: %.0f\n", edges_F))
  cat(sprintf("Archi in comune:           %.0f\n", intersect_edges))
  cat(sprintf("Jaccard Similarity:        %.4f\n", jaccard_index))
  cat(sprintf("SHD (Distance):            %.0f\n", shd_val))
  
  cat("\n--- 3. ANALISI HUBS (Max Degree) ---\n")
  cat(sprintf("Bayes Max: %.0f connessioni (Gene: %s)\n", max_deg_B, max_gene_B))
  cat(sprintf("Freq Max:  %.0f connessioni (Gene: %s)\n", max_deg_F, max_gene_F))
  
  cat("\n--- 4. VICINI DI 'logRiboflavin' (Target) ---\n")
  cat("BAYESIANO (", length(neigh_B), "): ", paste(neigh_B, collapse=", "), "\n")
  cat("FREQUENTISTA (", length(neigh_F), "): ", paste(neigh_F, collapse=", "), "\n")
  if(length(common_genes) > 0) {
    cat("GENI COMUNI: ", paste(common_genes, collapse=", "), "\n")
  } else {
    cat("GENI COMUNI: Nessuno\n")
  }
  
  cat("\n--- 5. CHECK GENI NOTI (Bühlmann et al. 2014) ---\n")
  if(length(selected_buhlmann) > 0) {
    cat("Geni noti presenti nello screening:", paste(selected_buhlmann, collapse=", "), "\n")
    for(gene in selected_buhlmann) {
      in_B <- ifelse(gene %in% neigh_B, "SI", "NO")
      in_F <- ifelse(gene %in% neigh_F, "SI", "NO")
      cat(sprintf(" - Gene %s -> Riboflavin? Bayes: %s | Freq: %s\n", gene, in_B, in_F))
    }
  } else {
    cat("Nessuno dei geni noti (YXLD_at, ecc.) è nei top", q_selected, "per varianza.\n")
  }
  
  cat("\n==========================================================\n")
  cat("                     FINE ANALISI                         \n")
  cat("==========================================================\n")
}