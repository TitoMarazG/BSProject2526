library(pcalg)
library(Rgraphviz)
library(igraph)
library(graph)

res_bayes <- readRDS("C:/Users/Ascolani/Documents/bayesian Proj/BSProject2526/test/RDS/Riboflavin_Bayes_q1000.rds")

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
  
  all_edge_names <- edgeNames(g_sub_nel)
  edge_cols <- character(length(all_edge_names))
  names(edge_cols) <- all_edge_names
  
  # 2. Cicliamo sulla matrice P (che abbiamo reso simmetrica)
  nodes_sub <- nodes(g_sub_nel)
  
  for (i in 1:length(nodes_sub)) {
    for (j in i:length(nodes_sub)) {
      node_i <- nodes_sub[i]
      node_j <- nodes_sub[j]
      
      # Se esiste un arco tra i e j nel grafo ridotto
      edge_key <- paste0(node_i, "~", node_j)
      rev_edge_key <- paste0(node_j, "~", node_i)
      
      # Troviamo quale delle due chiavi è quella corretta nel grafo
      actual_key <- NULL
      if (edge_key %in% all_edge_names) {
        actual_key <- edge_key
      } else if (rev_edge_key %in% all_edge_names) {
        actual_key <- rev_edge_key
      }
      
      if (!is.null(actual_key)) {
        p_val_raw <- P[node_i, node_j]
        p_norm <- (p_val_raw - 0.85) / (1 - 0.85)
        edge_cols[actual_key] <- gray(1 - p_norm)
      }
    }
  }
  
  edgeAttrs <- list(color = edge_cols)
  
  
  # Prepariamo gli attributi
  attrs <- list(
    graph = list(
      layout = "dot",
      rankdir = "BT",       # Alto -> Basso
      nodesep = "0.5",      # Spazio tra nodi dello stesso livello
      ranksep = "1"         # Spazio tra i livelli (altezza del grafo)
    ),
    edge = list(lwd = 3),
    node = list(fontsize = "10", fillcolor = "white", style = "filled")
  )
  
  # Per far risaltare il Riboflavin in basso, coloriamolo diversamente
  node_colors <- rep("white", length(nodes(g_sub_nel)))
  names(node_colors) <- nodes(g_sub_nel)
  node_colors["logRiboflavin"] <- "#E3F2FD" # Così lo vedi subito in fondo
  
  nAttrs$fillcolor <- node_colors
  
  plot(g_sub_nel, 
       nodeAttrs = nAttrs, 
       edgeAttrs = edgeAttrs, 
       attrs = attrs)
 
}
  

