#install.packages("hdi")
#install.packages("varbvs")

library(hdi)
library(varbvs)

data(riboflavin)

X <- riboflavin$x   # matrice delle covariate
y <- riboflavin$y   # vettore delle risposte


# X: n x p matrix, y: length n
Xsc <- scale(as.matrix(X))
ysc <- as.numeric(scale(y))

fit <- varbvs(X = Xsc, Z = NULL, y = ysc, family = "gaussian")

pip <- fit$pip  # P(gamma_j = 1 | y) posterior prob
names(pip) <- colnames(Xsc)

# Selezione delle covariate con posterior probability più alta di 0.5, oppure, se queste ultime sono poche, teniamo le prime 200
sel_pip <- names(pip)[pip > 0.5]          # soglia "mediana" classica -> include solo 4 covariate 
sel_top <- names(sort(pip, decreasing=TRUE))[1:200]  # tieni top 200

# Scegli una delle due:
selected <- if(length(sel_pip) >= 30) sel_pip else sel_top

X_reduced <- X[, selected, drop = FALSE]
