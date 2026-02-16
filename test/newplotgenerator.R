library(pcalg)
library(igraph)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(scales) 

# --- (ASSUMIAMO CHE I DATI 'P' E 'adj' SIANO GIA' CARICATI) ---
# Se necessario, riesegui la parte iniziale di caricamento dati.

# --- 1. DEFINIZIONE LIVELLI (Identificazione Automatica) ---
L0 <- target_node

# Livello 1
if (target_node %in% rownames(adj)) {
  L1 <- colnames(adj)[adj[target_node, ] == 1]
} else { L1 <- character(0) }

# Livello 2
if(length(L1) > 0) {
  L2_raw <- unique(unlist(lapply(L1, function(nm) {
    colnames(adj)[adj[nm, ] == 1]
  })))
  L2 <- setdiff(L2_raw, c(L0, L1))
} else { L2 <- character(0) }

nodes_keep <- c(L0, L1, L2)
adj_red = adj[nodes_keep, nodes_keep, drop = FALSE]
g_sub <- graph_from_adjacency_matrix(adj_red, mode = "undirected")

BF_sub <- BF_max[nodes_keep, nodes_keep, drop = FALSE]
BF_sub[is.na(BF_sub)] <- Inf 
BF_sub_sym <- pmin(BF_sub, t(BF_sub))
P <- 1 / (1 + BF_sub_sym)

t_graph <- as_tbl_graph(g_sub) %>%
  activate(nodes) %>%
  mutate(
    display_name = ifelse(name == "logRiboflavin", "RIBO", gsub("_at", "", name)),
    type = ifelse(name == "logRiboflavin", "Target", "Gene")
  ) %>%
  activate(edges) %>%
  mutate(
    from_name = .N()$name[from],
    to_name = .N()$name[to],
    prob = map2_dbl(from_name, to_name, ~ {
      val <- P[.x, .y]
      if(is.na(val)) 0 else val
    })
  )

# --- 2. POSIZIONAMENTO GEOMETRICO MIRATO ---
node_order <- V(g_sub)$name
n_nodes <- length(node_order)
x_coords <- numeric(n_nodes)
y_coords <- numeric(n_nodes)
get_idx <- function(names_vec) which(node_order %in% names_vec)

# --- A. LIVELLO 0 (Riboflavina) ---
idx_L0 <- get_idx(L0)
x_coords[idx_L0] <- 0
y_coords[idx_L0] <- 0

# --- B. LIVELLO 1 (YCKE e YDDK) ---
# Identifichiamo quale dei due è il "genitore" del Livello 2
# (Ovvero, chi tra i nodi L1 ha connessioni con i nodi L2?)
idx_L1 <- get_idx(L1)

# Controllo automatico: chi è collegato a L2?
is_parent <- sapply(L1, function(n1) any(adj[n1, L2] == 1))
parent_name <- names(is_parent)[is_parent]       # Dovrebbe essere YDDK
isolated_name <- names(is_parent)[!is_parent]    # Dovrebbe essere YCKE

idx_parent <- get_idx(parent_name)
idx_isolated <- get_idx(isolated_name)

# POSIZIONAMENTO L1 (Più vicini orizzontalmente: +/- 0.5)
# Mettiamo il genitore (YDDK) a destra (0.5) e l'altro a sinistra (-0.5)
x_coords[idx_parent] <- 0.5
x_coords[idx_isolated] <- -0.5

y_coords[idx_L1] <- 2 # Altezza Intermedia

# --- C. LIVELLO 2 (YDCL e RAPI) ---
# Devono essere equidistanti dal LORO genitore (YDDK che è a 0.5)
idx_L2 <- get_idx(L2)

# Calcoliamo le posizioni relative al genitore
# Genitore X = 0.5. 
# Figlio 1 va a (0.5 - 0.3) = 0.2
# Figlio 2 va a (0.5 + 0.3) = 0.8
parent_x <- x_coords[idx_parent]
spread <- 0.3 # Quanto si allargano i rami

# Assegniamo le coordinate centrate sul genitore
if(length(idx_L2) > 0) {
  x_coords[idx_L2] <- seq(parent_x - spread, parent_x + spread, length.out = length(idx_L2))
  y_coords[idx_L2] <- 4 # Altezza Massima
}

# --- 3. PLOT FINALE ---
ggraph(t_graph, layout = "manual", x = x_coords, y = y_coords) + 
  
  # Archi
  geom_edge_link(aes(color = prob, width = prob), alpha = 1) +
  
  scale_edge_color_gradientn(
    colors = c("#FFEDA0", "#FEB24C", "#F03B20", "#800026", "#4A004A"),
    name = "Probabilità (P)",
    limits = c(0.85, 1.0),   
    oob = scales::squish    
  ) +
  
  scale_edge_width(range = c(1, 4), guide = "none") +
  
  # Nodi
  geom_node_point(aes(fill = type), 
                  size = 22, 
                  shape = 21, 
                  color = "black", 
                  stroke = 1) + 
  
  scale_fill_manual(
    values = c("Target" = "#FFD700", "Gene" = "white"), 
    guide = "none"
  ) +
  
  # Testo
  geom_node_text(aes(label = display_name), 
                 size = 5, 
                 fontface = "bold", 
                 color = "black") +
  
  # Layout
  coord_cartesian(clip = "off") + 
  theme_void() + 
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey40"),
    plot.margin = margin(1, 2, 1, 2, "cm") 
  ) 
