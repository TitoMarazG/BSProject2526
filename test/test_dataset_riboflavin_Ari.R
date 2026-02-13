library(hdi)

library(pcalg)
library(Rgraphviz)
library(igraph)
source("C:/Users/arian/Desktop/Magistrale/Primo anno/Bayesian Statistics/Progetto/personalized_my_skeleton_senza_stable_fast.R")

data(riboflavin)
 
X <- as.matrix(riboflavin$x)   # matrice delle covariate
y <- as.numeric(riboflavin$y)   # vettore delle risposte

data = cbind(X[,300:350],y)
colnames(data)[ncol(data)] <- "Riboflavin"

data_sc = scale(data)

q = ncol(data)
alpha_bayes = 0.11 # BF associated to posterior probability of edge inclusion of 0.9 [given prior probability of edge inclusion of 0.5]

skltn.bayes <- my_skeleton(X = data_sc, a = q+1, U = diag(1, q), 
                           alpha = alpha_bayes, 
                           BayesTest = BF_Gaussian,
                           saveBF = TRUE, verbose = FALSE)

skltn <- skltn.bayes$object
grafo <- skltn@graph

plot(grafo)
adj <- as(grafo, "matrix")

# Vicini di logRiboflavin (Target)
# Essendo skeleton simmetrici, basta guardare la riga del target

target_node <- "Riboflavin"
neigh <- if(target_node %in% rownames(adj)) {
 names(which(adj[target_node, ] == 1))
} else { character(0) }

if(length(neigh)>0){
  pos_neigh = which(adj[target_node, ] == 1)
  
  # neighbors of neighbors
  neigh2 <- if(neigh %in% rownames(adj)) {
    names(which(adj[neigh, ] == 1))
  } else { character(0) }
  
  adj_red = adj[c(target_node,neigh,neigh2),c(target_node,neigh,neigh2)]
  
  g <- graph_from_adjacency_matrix(adj_red, mode = "undirected")
  g_sub <- induced_subgraph(g, vids = c(target_node,neigh,neigh2))
  # converti in graphNEL per plot
  g_sub_nel <- as_graphnel(g_sub)
  
  # plot "graphNEL style"
  plot(g_sub_nel)  # già molto meglio del plot igraph
  # plot(g_sub)
  
  BF = skltn.bayes$BF_history
  for(i in pos_neigh){
    BF_ribof_parents = BF[q][i]
    post_prob_edge_inclusion[i] = 1/(1 + BF_ribof_parents)
  }
}











