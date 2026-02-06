BiocManager::install("colonCA")
library(colonCA)
?colonCA
data(colonCA)
X <- t(exprs(colonCA))
dim(X)

# VANNO CENTRATI