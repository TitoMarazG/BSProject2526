# Questo script è per fare test su subsample a partire da un grosso sample

# SETUP INIZIALE================================================================

setwd("/Users/leomarcellopoli/Documents/Bayesian/Project/")

library(pcalg)
library(BiocGenerics)
library(BCDAG)
library(ggplot2)

rm(list = ls()) 
graphics.off() 

source("/Users/leomarcellopoli/Documents/Bayesian/Project/BSProject2526/funzioni R/personalized_my_skeleton_senza_stable_fast.R")
source("/Users/leomarcellopoli/Documents/Bayesian/Project/BSProject2526/funzioni R/compute_metrics.R")
source("/Users/leomarcellopoli/Documents/Bayesian/Project/BSProject2526/funzioni R/generate_DAG.R")

# ==============================================================================
# PARAMETRI PER SIMULAZIONE SUBSAMPLING
# ==============================================================================
q_list <- c(40)          # Numero di nodi nel DAG
n_max  <- 400           # Dimensione totale del dataset
n_folds <- 4             # Numero di subsample (1600 / 4 = 400 righe ciascuno)
sub_size <- n_max / n_folds

repetitions <- 1:20      # Numero di dataset diversi generati
w <- 0.6
seed = 3
alpha_freq <- 0.01
alpha_bayes <- 10

# Crea lista vuota per salvare i risultati
results_list <- list()

# LOOP ESTERNO: varia q
for (q in q_list) {
  
  cat("\n========================================\n")
  cat("Testing Subsamples with q =", q, "| Total n =", n_max, "| Folds =", n_folds, "\n")
  cat("========================================\n")
  
  # LOOP REPETITIONS (Generiamo 'repetition' scenari diversi)
  for (rep in repetitions) {
    
    cat("  Simulation: q =", q, ", repetition =", rep, "... ")
    
    # 1. Genera il DAG vero e il dataset COMPLETO (n_max)
    seed = seed + 1
    DD = generate_DAG(q = q, seed = seed, w = w, n = n_max)
    D0 = DD[[1]]      # Matrice di adiacenza vera
    X_full = DD[[2]]  # Dataset completo (1600 righe)
    
    # Skeleton vero (unico per tutti i subsample)
    D0_skeleton <- D0 + t(D0)
    D0_skeleton[D0_skeleton > 0] <- 1
    skltn.true <- as(D0_skeleton, "graphNEL")
    
    # 2. LOOP SUI SUBSAMPLE (Tagliamo X_full a fette)
    for (k in 1:n_folds) {
      
      # Calcola indici per il subsample k
      # Es: k=1 -> 1:400; k=2 -> 401:800; ...
      idx_start <- (k - 1) * sub_size + 1
      idx_end   <- k * sub_size
      
      # Estrai il subsample
      X_sub <- X_full[idx_start:idx_end, ]
      
      # --- A. Stima Frequentista ---
      skltn.freq <- skeleton(suffStat = list(C = cor(X_sub), n = nrow(X_sub)), 
                             indepTest = gaussCItest,
                             p = ncol(X_sub),
                             alpha = alpha_freq)
      
      # --- B. Stima Bayesiana ---
      skltn.bayes <- my_skeleton(X = X_sub, a = q+1, U = diag(1, q), 
                                 alpha = alpha_bayes, 
                                 BayesTest = BF_Gaussian,
                                 saveBF = FALSE, verbose = FALSE)
      
      # --- C. Calcola metriche ---
      eval_freq <- evaluate_graph(skltn.true, skltn.freq@graph, undirected = TRUE)
      eval_bayes <- evaluate_graph(skltn.true, skltn.bayes@graph, undirected = TRUE)
      
      # --- D. Salva risultati ---
      results_list[[length(results_list) + 1]] <- data.frame(
        q = q,
        n_total = n_max,
        sub_size = sub_size,
        subsample_id = k,         # Identificativo del subsample (1, 2, 3, 4)
        repetition = rep,
        method = c("Frequentist", "Bayesian"),
        SHD = c(eval_freq$SHD, eval_bayes$SHD),
        Precision = c(eval_freq$Precision, eval_bayes$Precision),
        Recall = c(eval_freq$Recall, eval_bayes$Recall),
        F1_Score = c(eval_freq$F1_Score, eval_bayes$F1_Score)
      )
    } # Fine loop subsamples
    
    cat("Done!\n")
  } # Fine loop repetitions
}

# Risultati in un unico data frame
results_df <- do.call(rbind, results_list)
results_df <- as.data.frame(results_df)
rownames(results_df) <- NULL

# Trasformazioni per il plot
df_plot <- transform(
  results_df,
  subsample_id = factor(subsample_id), # Ora l'asse X sarà il subsample
  q = factor(q),
  method = factor(method)
)

View(results_df)

# ==============================================================================
# FUNZIONE DI PLOTTING (Adattata per Subsamples)
# ==============================================================================
plot_subsample_metric <- function(data, y_var, y_label, title_suffix) {
  
  p <- ggplot(data, aes(x = subsample_id, y = .data[[y_var]], fill = method)) +
    geom_boxplot(
      width = 0.7,
      position = position_dodge(width = 0.8),
      outlier.shape = NA
    ) +
    geom_jitter(
      aes(color = method),
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      alpha = 0.35, size = 1,
      show.legend = FALSE
    ) +
    labs(title = paste("Stability across Subsamples -", title_suffix),
         subtitle = paste("n_max =", n_max, "| Fold size =", sub_size, "| q =", q_list),
         x = "Subsample ID (Disjoint Sets)",
         y = y_label,
         fill = "Method:") +
    theme_bw() +
    theme(
      legend.position = "top",
      legend.justification = "center",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    guides(
      fill = guide_legend(
        title.position = "top", 
        title.hjust = 0.5,      
        direction = "vertical", 
        nrow = 2                
      )
    )
  
  return(p)
}

# ==============================================================================
# GENERAZIONE GRAFICI
# ==============================================================================

metrics_to_plot <- list(
  list(var = "Precision", label = "Precision"),
  list(var = "Recall", label = "Recall"),
  list(var = "F1_Score", label = "F1 Score"),
  list(var = "SHD", label = "SHD to true CPDAG")
)

for (m in metrics_to_plot) {
  cat(paste("\nGenerazione plot per:", m$var, "\n"))
  p <- plot_subsample_metric(df_plot, 
                             y_var = m$var, 
                             y_label = m$label, 
                             title_suffix = m$var)
  print(p)
}
