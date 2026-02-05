#==============================================================================#
# PLOT PERFORMANCES #
#==============================================================================#


##### SHD #####

## PRINT CONFRONTO BAYESIAN VS FREQUENTIST PER TUTTI I q INSIEME #
#
#p <- ggplot(df_plot, aes(x = n, y = SHD, fill = method)) +
#  geom_boxplot(
#    width = 0.7,
#    position = position_dodge(width = 0.8),
#    outlier.shape = NA # per non mostrare gli outlier nei boxplot
#  ) +
#  geom_jitter( # funzione che aggiunge i punti nei boxplot
#    aes(color = method),
#    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), # dodge allinea i punti al boxplot giusto
#    # jitter li sposta lateralmente per non farli sovrapporre alla linea dei boxplot
#    alpha = 0.35, size = 1,
#    show.legend = FALSE
#  ) +
#  facet_wrap(~ q, nrow = 1) +
#  labs(title = "Confronto performance SHD",
#       x = "Sample size (n)",
#       y = "SHD to true CPDAG",
#       fill = "Method:") +
#  theme_bw() +
#  theme(
#    legend.position = "top",
#    legend.justification = "center",
#    legend.title = element_text(face = "bold"),
#    panel.grid.minor = element_blank(),
#    plot.title = element_text(hjust = 0.5, face = "bold") 
#  ) +
#  guides(
#    fill = guide_legend(
#      title.position = "top",  # Mette "Method:" sopra le voci
#      title.hjust = 0.5,       # Centra la scritta "Method:" rispetto ai box
#      direction = "vertical",  # Dispone Bayesian e Frequentist in verticale
#      nrow = 2                 # Forza due righe
#    )
#  )
#
#print(p)
#
## PRINT SEPARATI PER q #
#
#for (q_val in q_list) {
#  p <- ggplot(subset(df_plot, q == q_val), aes(x = n, y = SHD, fill = method)) +
#    geom_boxplot(
#      width = 0.7,
#      position = position_dodge(width = 0.8),
#      outlier.shape = NA # per non mostrare gli outlier nei boxplot
#    ) +
#    geom_jitter( # funzione che aggiunge i punti nei boxplot
#      aes(color = method),
#      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), # dodge allinea i punti al boxplot giusto
#      # jitter li sposta lateralmente per non farli sovrapporre alla linea dei boxplot
#      alpha = 0.35, size = 1,
#      show.legend = FALSE
#    ) +
#    labs(title = paste("Confronto performance SHD per q =", q_val),
#         x = "Sample size (n)",
#         y = "SHD to true CPDAG",
#         fill = "Method:") +
#    theme_bw() +
#    theme(
#      legend.position = "top",
#      legend.justification = "center",
#      legend.title = element_text(face = "bold"),
#      panel.grid.minor = element_blank(),
#      plot.title = element_text(hjust = 0.5, face = "bold") 
#    ) +
#    guides(
#      fill = guide_legend(
#        title.position = "top",  # Mette "Method:" sopra le voci
#        title.hjust = 0.5,       # Centra la scritta "Method:" rispetto ai box
#        direction = "vertical",  # Dispone Bayesian e Frequentist in verticale
#        nrow = 2                 # Forza due righe
#      )
#    )
#  
#  print(p)
#}


##### Precision #####

## PRINT CONFRONTO BAYESIAN VS FREQUENTIST PER TUTTI I q INSIEME #
#
#p_precision <- ggplot(df_plot, aes(x = n, y = Precision, fill = method)) +
#  geom_boxplot(
#    width = 0.7,
#    position = position_dodge(width = 0.8),
#    outlier.shape = NA # per non mostrare gli outlier nei boxplot
#  ) +
#  geom_jitter( # funzione che aggiunge i punti nei boxplot
#    aes(color = method),
#    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), # dodge allinea i punti al boxplot giusto
#    # jitter li sposta lateralmente per non farli sovrapporre alla linea dei boxplot
#    alpha = 0.35, size = 1,
#    show.legend = FALSE
#  ) +
#  facet_wrap(~ q, nrow = 1) +
#  labs(title = "Confronto performance Precision",
#       x = "Sample size (n)",
#       y = "Precision",
#       fill = "Method:") +
#  theme_bw() +
#  theme(
#    legend.position = "top",
#    legend.justification = "center",
#    legend.title = element_text(face = "bold"),
#    panel.grid.minor = element_blank(),
#    plot.title = element_text(hjust = 0.5, face = "bold") 
#  ) +
#  guides(
#    fill = guide_legend(
#      title.position = "top",  # Mette "Method:" sopra le voci
#      title.hjust = 0.5,       # Centra la scritta "Method:" rispetto ai box
#      direction = "vertical",  # Dispone Bayesian e Frequentist in verticale
#      nrow = 2                 # Forza due righe
#    )
#  )
#
#print(p_precision)
#
## PRINT SEPARATI PER q #
#
#for (q_val in q_list) {
#  p_precision <- ggplot(subset(df_plot, q == q_val), aes(x = n, y = Precision, fill = method)) +
#    geom_boxplot(
#      width = 0.7,
#      position = position_dodge(width = 0.8),
#      outlier.shape = NA # per non mostrare gli outlier nei boxplot
#    ) +
#    geom_jitter( # funzione che aggiunge i punti nei boxplot
#      aes(color = method),
#      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), # dodge allinea i punti al boxplot giusto
#      # jitter li sposta lateralmente per non farli sovrapporre alla linea dei boxplot
#      alpha = 0.35, size = 1,
#      show.legend = FALSE
#    ) +
#    labs(title = paste("Confronto performance Precision per q =", q_val),
#         x = "Sample size (n)",
#         y = "Precision",
#         fill = "Method:") +
#    theme_bw() +
#    theme(
#      legend.position = "top",
#      legend.justification = "center",
#      legend.title = element_text(face = "bold"),
#      panel.grid.minor = element_blank(),
#      plot.title = element_text(hjust = 0.5, face = "bold") 
#    ) +
#    guides(
#      fill = guide_legend(
#        title.position = "top",  # Mette "Method:" sopra le voci
#        title.hjust = 0.5,       # Centra la scritta "Method:" rispetto ai box
#        direction = "vertical",  # Dispone Bayesian e Frequentist in verticale
#        nrow = 2                 # Forza due righe
#      )
#    )
#  
#  print(p_precision)
#}

##### Recall #####

## PRINT CONFRONTO BAYESIAN VS FREQUENTIST PER TUTTI I q INSIEME #
#
#p_recall <- ggplot(df_plot, aes(x = n, y = Recall, fill = method)) +
#  geom_boxplot(
#    width = 0.7,
#    position = position_dodge(width = 0.8),
#    outlier.shape = NA # per non mostrare gli outlier nei boxplot
#  ) +
#  geom_jitter( # funzione che aggiunge i punti nei boxplot
#    aes(color = method),
#    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), # dodge allinea i punti al boxplot giusto
#    # jitter li sposta lateralmente per non farli sovrapporre alla linea dei boxplot
#    alpha = 0.35, size = 1,
#    show.legend = FALSE
#  ) +
#  facet_wrap(~ q, nrow = 1) +
#  labs(title = "Confronto performance Recall",
#       x = "Sample size (n)",
#       y = "Recall",
#       fill = "Method:") +
#  theme_bw() +
#  theme(
#    legend.position = "top",
#    legend.justification = "center",
#    legend.title = element_text(face = "bold"),
#    panel.grid.minor = element_blank(),
#    plot.title = element_text(hjust = 0.5, face = "bold") 
#  ) +
#  guides(
#    fill = guide_legend(
#      title.position = "top",  # Mette "Method:" sopra le voci
#      title.hjust = 0.5,       # Centra la scritta "Method:" rispetto ai box
#      direction = "vertical",  # Dispone Bayesian e Frequentist in verticale
#      nrow = 2                 # Forza due righe
#    )
#  )
#
#print(p_recall)
#
## PRINT SEPARATI PER q #
#
#for (q_val in q_list) {
#  p_recall <- ggplot(subset(df_plot, q == q_val), aes(x = n, y = Recall, fill = method)) +
#    geom_boxplot(
#      width = 0.7,
#      position = position_dodge(width = 0.8),
#      outlier.shape = NA # per non mostrare gli outlier nei boxplot
#    ) +
#    geom_jitter( # funzione che aggiunge i punti nei boxplot
#      aes(color = method),
#      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), # dodge allinea i punti al boxplot giusto
#      # jitter li sposta lateralmente per non farli sovrapporre alla linea dei boxplot
#      alpha = 0.35, size = 1,
#      show.legend = FALSE
#    ) +
#    labs(title = paste("Confronto performance Recall per q =", q_val),
#         x = "Sample size (n)",
#         y = "Recall",
#         fill = "Method:") +
#    theme_bw() +
#    theme(
#      legend.position = "top",
#      legend.justification = "center",
#      legend.title = element_text(face = "bold"),
#      panel.grid.minor = element_blank(),
#      plot.title = element_text(hjust = 0.5, face = "bold") 
#    ) +
#    guides(
#      fill = guide_legend(
#        title.position = "top",  # Mette "Method:" sopra le voci
#        title.hjust = 0.5,       # Centra la scritta "Method:" rispetto ai box
#        direction = "vertical",  # Dispone Bayesian e Frequentist in verticale
#        nrow = 2                 # Forza due righe
#      )
#    )
#  
#  print(p_recall)
#}

##### F1_Score #####

## PRINT CONFRONTO BAYESIAN VS FREQUENTIST PER TUTTI I q INSIEME #
#
#p_f1 <- ggplot(df_plot, aes(x = n, y = F1_Score, fill = method)) +
#  geom_boxplot(
#    width = 0.7,
#    position = position_dodge(width = 0.8),
#    outlier.shape = NA # per non mostrare gli outlier nei boxplot
#  ) +
#  geom_jitter( # funzione che aggiunge i punti nei boxplot
#    aes(color = method),
#    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), # dodge allinea i punti al boxplot giusto
#    # jitter li sposta lateralmente per non farli sovrapporre alla linea dei boxplot
#    alpha = 0.35, size = 1,
#    show.legend = FALSE
#  ) +
#  facet_wrap(~ q, nrow = 1) +
#  labs(title = "Confronto performance F1 Score",
#       x = "Sample size (n)",
#       y = "F1 Score",
#       fill = "Method:") +
#  theme_bw() +
#  theme(
#    legend.position = "top",
#    legend.justification = "center",
#    legend.title = element_text(face = "bold"),
#    panel.grid.minor = element_blank(),
#    plot.title = element_text(hjust = 0.5, face = "bold") 
#  ) +
#  guides(
#    fill = guide_legend(
#      title.position = "top",  # Mette "Method:" sopra le voci
#      title.hjust = 0.5,       # Centra la scritta "Method:" rispetto ai box
#      direction = "vertical",  # Dispone Bayesian e Frequentist in verticale
#      nrow = 2                 # Forza due righe
#    )
#  )
#
#print(p_f1)
#
## PRINT SEPARATI PER q #
#
#for (q_val in q_list) {
#  p_f1 <- ggplot(subset(df_plot, q == q_val), aes(x = n, y = F1_Score, fill = method)) +
#    geom_boxplot(
#      width = 0.7,
#      position = position_dodge(width = 0.8),
#      outlier.shape = NA # per non mostrare gli outlier nei boxplot
#    ) +
#    geom_jitter( # funzione che aggiunge i punti nei boxplot
#      aes(color = method),
#      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), # dodge allinea i punti al boxplot giusto
#      # jitter li sposta lateralmente per non farli sovrapporre alla linea dei boxplot
#      alpha = 0.35, size = 1,
#      show.legend = FALSE
#    ) +
#    labs(title = paste("Confronto performance F1 Score per q =", q_val),
#         x = "Sample size (n)",
#         y = "F1 Score",
#         fill = "Method:") +
#    theme_bw() +
#    theme(
#      legend.position = "top",
#      legend.justification = "center",
#      legend.title = element_text(face = "bold"),
#      panel.grid.minor = element_blank(),
#      plot.title = element_text(hjust = 0.5, face = "bold") 
#    ) +
#    guides(
#      fill = guide_legend(
#        title.position = "top",  # Mette "Method:" sopra le voci
#        title.hjust = 0.5,       # Centra la scritta "Method:" rispetto ai box
#        direction = "vertical",  # Dispone Bayesian e Frequentist in verticale
#        nrow = 2                 # Forza due righe
#      )
#    )
#  
#  print(p_f1)
#}


# ... (tutta la parte di simulazione resta uguale) ...
