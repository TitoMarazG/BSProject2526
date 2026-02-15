library(pcalg)
library(Rgraphviz)
library(igraph)
library(graph)

res_bayes <- readRDS("C:/Users/arian/Documents/GitHub/BSProject2526/test/RDS/Riboflavin_Bayes_q1000.rds")

skltn.bayes <- res_bayes@graph
adj <- as(skltn.bayes, "matrix")

BF_max = res_bayes@pMax
all_node_names <- nodes(res_bayes@graph)
rownames(BF_max) <- all_node_names
colnames(BF_max) <- all_node_names

target_node <- "logRiboflavin"

if (target_node %in% rownames(adj)) {
  neigh <- colnames(adj)[adj[target_node, ] == 1]
} else { character(0) }

if(length(neigh)>0){
  pos_neigh = which(adj[target_node, ] == 1)
  
  # neighbors of neighbors
  neigh2 <- unique(unlist(lapply(neigh, function(nm) {
    colnames(adj)[adj[nm, ] == 1]
  })))
  
  nodes_keep <- unique(c(target_node, neigh, neigh2))
  adj_red = adj[nodes_keep, nodes_keep, drop = FALSE]
  
  # subgrafo
  g_sub <- graph_from_adjacency_matrix(adj_red, mode = "undirected")
  # converti in graphNEL per plot
  g_sub_nel <- as_graphnel(g_sub)    
  
  
  # P = matrice di posterior probability of edge inclusion
  nodes_in_sub = nodes(g_sub_nel)
  BF_sub <- BF_max[nodes_keep, nodes_keep, drop = FALSE]
  # Rendiamo la matrice simmetrica (necessario per grafi non orientati)
  BF_sub[is.na(BF_sub)] <- Inf # Se non c'è evidenza, BF è infinito (P -> 0)
  BF_sub_sym <- pmin(BF_sub, t(BF_sub)) # Usiamo il minimo BF per l'inclusione più probabile
  P <- 1 / (1 + BF_sub_sym) # ipotizzando prior probability of edge inclusion = 0.5
  diag(P) <- NA
  
  # ------ PLOT ------
  
  # cambio la label del nodo logriboflavin in Riboflavin
  lbls <- nodes(g_sub_nel)
  lbls[lbls == "logRiboflavin"] <- "Riboflavin"
  nAttrs <- list(label = setNames(lbls, nodes(g_sub_nel)))
  
  # # cambio colore edges in base posterior probability of edge inclusion
  
  edge_list <- edgeL(g_sub_nel)
  edge_cols <- character()
  
  for (node_from in names(edge_list)) {
    tos <- edge_list[[node_from]]$edges
    if (length(tos) == 0) next
    
    node_tos <- nodes(g_sub_nel)[tos]
    
    for (node_to in node_tos) {
      # Per grafi non orientati, processiamo l'arco solo una volta (ordine alfabetico)
      if (node_from < node_to) {
        # Recuperiamo la probabilità dalla matrice P
        # Se P non è simmetrica, potresti voler usare max(P[from,to], P[to,from])
        p_val <- (P[node_from, node_to] - 0.9)/(1-0.9)
        p_val <- p_val^2 #per aumentare contrasto
        
        # Gestione NA o valori fuori range
        if (is.na(p_val)) p_val <- 0
        
        # Creazione colore: 1 = Nero (p=1), 0.9 = Grigio quasi bianco (p=0)
        # alpha determina la trasparenza, o puoi usare gray()
        col_hex <- gray(1 - p_val) 
        
        # La chiave DEVE essere "nodo1~nodo2"
        edge_key <- paste0(node_from, "~", node_to)
        edge_cols[edge_key] <- col_hex
      }
    }
  }
  
  edgeAttrs <- list(color = edge_cols)
 
  # # # ---- Plot finale "graphNEL style" ---- 
  plot(g_sub_nel, nodeAttrs = nAttrs, edgeAttrs = edgeAttrs, attrs = list(edge = list(lwd = 2)))
  #plot(g_sub_nel, nodeAttrs = nAttrs)
}

