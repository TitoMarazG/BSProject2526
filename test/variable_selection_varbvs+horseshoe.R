library(hdi)
library(varbvs)
library(rstanarm)

data(riboflavin)

X <- riboflavin$x   # matrice delle covariate
y <- riboflavin$y   # vettore delle risposte

# X: n x p matrix, y: length n
Xsc <- scale(as.matrix(X))
ysc <- as.numeric(scale(y))

fit <- varbvs(X = Xsc, Z = NULL, y = ysc, family = "gaussian")
pip <- fit$pip
names(pip) <- colnames(Xsc)
keep <- names(sort(pip, decreasing = TRUE))[1:300]

X_red <- Xsc[, keep, drop=FALSE]
df <- data.frame(y = ysc, X_red)

# contiene campioni posteriori
fit_hs <- stan_glm(
  y ~ .,        # tutte le altre colonne del data.frame (colonne di X_red)
  data = df,
  prior = hs(global_scale = 1),
  chains = 4, iter = 2000, seed = 1  # 4 catene MCMC indipendenti, 2000 iterazioni per chain
)

post <- as.matrix(fit_hs)  # Righe = campioni MCMC (tutte le chain concatenate); Colonne = coefficienti dei predittori

B <- post[, !colnames(post) %in% c("(Intercept)", "sigma"), drop = FALSE]

# il problema è che mi si sono sballati tutti i nomi delle colonne dopo aver usato STAN che quindi non coincidono più con quelli di Xsc (all'inizio del nome ha aggiunto X_red)
B_names <- colnames(B)
# rimuovi il prefisso X_red
B_names_new <- sub("^X_red", "", B_names)
colnames(B) <- B_names_new

# P(|beta| > eps | y): probabilità a posteriori che coeff di regressione sia “sufficientemente lontano da 0”
eps <- 0.1  # in scala standardizzata: 0.1 è un effetto "piccolo ma reale"
p_big <- colMeans(abs(B) > eps) # media per colonna dei TRUE/FALSE (1/0): circa probabilità a posteriori che beta sia > soglia

# Selezione:
selected <- names(sort(p_big, decreasing = TRUE))[1:50]
X_red <- Xsc[, selected, drop = FALSE]  

