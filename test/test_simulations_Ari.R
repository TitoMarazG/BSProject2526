# SETUP INIZIALE================================================================

setwd("/Users/leomarcellopoli/Documents/Bayesian/Project/BSProject2526")

library(pcalg)
library(BiocGenerics)
library(BCDAG)
library(ggplot2)

rm(list = ls()); graphics.off(); cat("\014")

source("/Users/leomarcellopoli/Documents/Bayesian/Project/BSProject2526/funzioni R/personalized_my_skeleton_senza_stable_fast.R")
source("/Users/leomarcellopoli/Documents/Bayesian/Project/BSProject2526/funzioni R/compute_metrics.R")
source("/Users/leomarcellopoli/Documents/Bayesian/Project/BSProject2526/funzioni R/generate_DAG.R")

# Liste di parametri da testare
# q_list <- c(1000)                 # Nota per me stesso: con q = 1000, n = 100 ci mette un minuto e mezzo a ripetizione
# n_list <- c(100)
q_list <- c(50)#, 200, 1000)
n_list <- c(100, 200, 500)
repetitions <- 1:40
#w <- 0.3  viene data in input come 2/q
seed = 3
alpha_freq <- 0.01
p_posterior <- 0.9
alpha_bayes <- 1/p_posterior - 1  # come prior mettiamo

# Crea lista vuota per salvare i risultati
results_list <- list()

# LOOP ESTERNO: varia q
for (q in q_list) {
  
  cat("\n========================================\n")
  cat("Testing with q =", q, "nodes\n")
  cat("========================================\n")
  
  # LOOP INTERNO: varia n
  for (n in n_list) {
    
    # LOOP repetitions (stessa (n, q) per #repetitions volte)
    
    for (rep in repetitions) {
      
      cat("  Running simulation: q =", q, ", n =", n, ", repetition =", rep, "... ")
      
      # Genera DAG e dati
      seed = seed + 1
      DD = generate_DAG(q = q, seed = seed, w = 2/q, n = n)
      D0 = DD[[1]]
      X = DD[[2]]
      
      # Skeleton vero
      D0_skeleton <- D0 + t(D0)
      D0_skeleton[D0_skeleton > 0] <- 1
      skltn.true <- as(D0_skeleton, "graphNEL")
      
      # Stima frequentista
      skltn.freq <- skeleton(suffStat = list(C = cor(X), n = nrow(X)), 
                             indepTest = gaussCItest,
                             p = ncol(X),
                             alpha = alpha_freq)
      
      # Stima bayesiana
      #res.bayes <- my_skeleton(X = X, a = q, U = U, 
      #                         alpha = alpha_bayes, 
      #                         BayesTest = BF_Gaussian,
      #                         saveBF = TRUE, verbose = FALSE)
      #skltn.bayes <- res.bayes$object
      skltn.bayes <- my_skeleton(X = X, a = q+1, U = diag(1, q), 
                                 alpha = alpha_bayes, 
                                 BayesTest = BF_Gaussian,
                                 saveBF = FALSE, verbose = FALSE)
      
      # Calcola metriche
      eval_freq <- evaluate_graph(skltn.true, skltn.freq@graph, undirected = TRUE)
      eval_bayes <- evaluate_graph(skltn.true, skltn.bayes@graph, undirected = TRUE)
      
      # Salva risultati
      results_list[[length(results_list) + 1]] <- data.frame(
        q = q,
        n = n,
        method = c("Frequentist", "Bayesian"),
        SHD = c(eval_freq$SHD, eval_bayes$SHD),
        Precision = c(eval_freq$Precision, eval_bayes$Precision),
        Recall = c(eval_freq$Recall, eval_bayes$Recall),
        F1_Score = c(eval_freq$F1_Score, eval_bayes$F1_Score),
        TP = c(eval_freq$TP, eval_bayes$TP),
        FP = c(eval_freq$FP, eval_bayes$FP),
        FN = c(eval_freq$FN, eval_bayes$FN),
        TN = c(eval_freq$TN, eval_bayes$TN),
        test_dataset = repetitions[rep]
      )
      
      cat("Done!\n")
    }
  }
}

# Risultati in un unico data frame
results_df <- do.call(rbind, results_list)
results_df <- as.data.frame(results_df)
rownames(results_df) <- NULL

View(results_df)
#write.csv(results_df, "simulation_results.csv", row.names = FALSE)

df_plot <- transform(
  results_df,
  n = factor(n),
  q = factor(q),
  method = factor(method)
)

##### FUNZIONE DI PLOTTING UNIVERSALE #####
plot_simulation_metric <- function(data, y_var, y_label, title_suffix, q_filter = NULL) {
  
  # Se q_filter è specificato, filtra i dati
  plot_data <- if (!is.null(q_filter)) subset(data, q == q_filter) else data
  
  p <- ggplot(plot_data, aes(x = n, y = .data[[y_var]], fill = method)) +
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
    # Aggiungi titolo e assi
    labs(title = paste("Confronto performance", title_suffix),
         subtitle = paste("Prob. edge generation =", 2/q, "| Bayes threshold =", round(alpha_bayes, digits = 2), "| Significance level =", alpha_freq),
         x = "Sample size (n)",
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
  
  if (is.null(q_filter)) {   # Aggiungi facet_wrap SOLO se stiamo plottando tutti i q insieme
    p <- p + facet_wrap(~ q, nrow = 1)
  }
  
  return(p)
}

##### GENERAZIONE GRAFICI #####

# Crea una cartella "Plots" se non esiste già
if (!dir.exists("Plots varying q & n")) {
  dir.create("Plots varying q & n")
}

# Definiamo le metriche da stampare
metrics_to_plot <- list(
  list(var = "Precision", label = "Precision"),
  list(var = "Recall", label = "Recall"),
  list(var = "F1_Score", label = "F1 Score"),
  list(var = "SHD", label = "SHD to true CPDAG")
)

# CICLO UNICO PER TUTTE LE METRICHE
for (m in metrics_to_plot) {
  
  # Stampa i grafici SINGOLI per ogni q
  for (q_val in q_list) {
    cat(paste("\nGenerazione plot per:", m$var, "\n con q =", q_val))
    p_single <- plot_simulation_metric(df_plot, 
                                       y_var = m$var, 
                                       y_label = m$label, 
                                       title_suffix = paste(m$var, "per q =", q_val),
                                       q_filter = q_val)
    print(p_single)
    
    # # Costruzione nome file: es. "Plots/Precision_q10.png"
    # filename_single <- paste0("Plots/", m$var, "_q", q_val, ".png")
    # 
    # # Salvataggio
    # ggsave(filename = filename_single, plot = p_single, 
    #        width = 7, height = 5, dpi = 300)
  }
  
  # Stampa il grafico GLOBALE (con facet_wrap)
  cat(paste("\nGenerazione plot globale per:", m$var))
  p_global <- plot_simulation_metric(df_plot, 
                                     y_var = m$var, 
                                     y_label = m$label, 
                                     title_suffix = m$var) # Titolo es: "Confronto performance SHD"
  print(p_global)
  
  # # Costruzione nome file: es. "Plots/Precision_Global_Comparison.png"
  # filename_global <- paste0("Plots/", m$var, "_Global_Comparison.png")
  # 
  # # Salvataggio (lo facciamo più largo perché contiene più grafici affiancati)
  # ggsave(filename = filename_global, plot = p_global, 
  #        width = 12, height = 6, dpi = 300)

  cat(paste("\n"))
}
