# Dopo aver letto il paper ho buttato giù uno script che fa vedere come 
#   funziona il pacchetto pcalg.
# Tutto quello che c'è scritto qui dentro è stato preso dal paper Pcalg_JSS.

#### INSTALLAZIONE pcalg package ####
# To install, you first need to install the package that manages the
#   Bioconductor repositories, as this is where the packages "graph", "RBGL"
#   and "Rgraphviz", which pcalg heavily depends on, are located.
#install.packages("BiocManager")
#BiocManager::install(c("graph", "RBGL", "Rgraphviz"))

# Now I install "pcalg"
#install.packages("pcalg")

# Necessary libraries to us pcalg
library("pcalg")
library("graph")
library("RBGL")
library("Rgraphviz")

#### Example of how to use pcalg ####

# Starting objects needed:
#   - a function to compute conditional independence test,
#   - a summary of the data on which the conditional independence test function
#   can work.
# We will use different test as conditional independence function, and the
#   correlation matrix of the data and the sample size as sufficient statistic.
# Note: in the package there are different independence test for different type
#   of data:
#    - gaussCItest() gaussian data
#    - disCItest() discrete data
#    - binCItest() binary data

#### Estimation of the skeleton (undirected graph) ####
# skeleton(sufficient statistic,
#          conditional independence test,
#          number of features/variables,
#          significance level,
#          ...)

## Example of Gaussian test
data("gmG")
skeleton.fit.gaussian <- skeleton(suffStat = list(C = cor(gmG$x),
                                                  n = nrow(gmG$x)),
                                  indepTest = gaussCItest,
                                  p = ncol(gmG$x),
                                  alpha = 0.01)
par(mfrow = c(1,2))
plot(gmG$g, main = "")
plot(skeleton.fit.gaussian, main = "")
dev.off()

## Example: discrete variables
data("gmD")
data_matrix <- data.matrix(gmD$x)
skeleton.fit.discrete <- skeleton(suffStat = list(dm = data_matrix,
                                                  n = unlist(lapply(gmD$x, nlevels)), #nrow(data_matrix),
                                                  adaptDF = TRUE),
                                   indepTest = disCItest,
                                   p = ncol(data_matrix),
                                   alpha = 0.01)
par(mfrow = c(1,2))
plot(gmD$g, main = "")
plot(skeleton.fit.discrete, main = "")
dev.off()

#### Use pc(): ####
# pc(sufficient statistic,
#    conditional independence test,
#    number of features/variables,
#    significance level,
#    ...)
CPDAG.fit.gaussian <- pc(suffStat = list(C = cor(gmG$x),
                                         n = nrow(gmG$x)),
                         indepTest = gaussCItest,
                         p = ncol(gmG$x),
                         alpha = 0.01)
par(mfrow = c(1,2))
plot(gmG$g, main = "")
plot(CPDAG.fit.gaussian, main = "")
dev.off()


#### Use ida() and idaFast() ####
# To see the causal effect of an intervention we implement:
# ida("cause" variable,
#     "target" variable,
#      covariance matrix,
#      the estimated graph)
ida(1, 6, cov(gmG$x), CPDAG.fit.gaussian@graph)
# [1] 0.7536376 0.5487757
# This means that given the uncertainity of the estimated DAG, these two are two
#   possible causal structures between between V1 and V6.
#
# Now let's compute the effect of V1 on V4, V5 and V6: two possible way:
#ida(1, 4, cov(gmG$x), pc.fit@graph)
#ida(1, 5, cov(gmG$x), pc.fit@graph)
#ida(1, 6, cov(gmG$x), pc.fit@graph)
# or
idaFast(1, c(4,5,6), cov(gmG$x), CPDAG.fit.gaussian@graph)
#        [,1]       [,2]
# 4 0.0102699 0.01201369
# 5 0.2387461 0.01793453
# 6 0.7536376 0.54877566
# Note: the true values for the causal effects are 0, 0.2, 0.71 (since it's a
#   simulation we were able to compute them)

# Other example of ida()
data("gmI")
pc.fit.ida <- pc(suffStat = list(C = cor(gmI$x),
                                 n = nrow(gmI$x)),
                 indepTest = gaussCItest,
                 p = ncol(gmI$x),
                 alpha = 0.01)
par(mfrow = c(1,2))
plot(gmI$g, main = "")
plot(pc.fit.ida, main = "")
ida(2, 5, cov(gmI$x), pc.fit.ida@graph, method = "global", verbose = TRUE)
# by putting method = "global", return the value for each of the DAG in the
# (estimated) equivalence class; verbose = TRUE returns details on the
# regressions 
ida(2, 5, cov(gmI$x), pc.fit.ida@graph, method = "local")
dev.off()

# idaFast -> This function estimates the multiset of possible total causal
#            effects of one variable (x) on a several (i.e., a vector of) target
#            variables (y) from observational data.
idaFast(2, c(5,6,7), cov(gmI$x), pc.fit.ida@graph)

#### Other conditional independence test ####
# The user can define himself any conditional independence test to use inside
#   skeleton(), pc(), ida() or idaFast(), ecc..
# The function must be in the form:
#     indepTest(x, y, S, suffStat)
#   and it must return the p-value of the conditional independence test of V_x
#   and V_y given V_S, where x, y and S indicate column positions of the
#   original data matrix, and some information on the data in the form of a
#   sufficient statistic (like dataframe of the original dataset or the
#   cov.matrix).
# A way to get the partial correlation of V_x and V_y given V_S is by solving
#   the two linear regressione problem
#     lm(V_x ~ V_S) and
#     lm(V_y ~ V_S),
#   get the residuals and calculate the corr between the the residuals...
# suffStat is an object containing several pieces of info all used within the
#   the function to prooduce the p-value.
#
# Example:
#   how myCItest works:
#     if S is empty, then it execute a standard correlation test between x and
#     y
#     if S is NOT empty, then it computes residuals of the linear regressions
#     x ~ S and y ~ S, then it execute a standard correlation test between 
#     residuals of x and y
#     then it gives back the p-value of the test
myCItest <- function(x, y, S, suffStat){
  if(length(S) == 0){
    x. <- suffStat[, x]
    y. <- suffStat[, y]
  }else{
    rxy <- resid(lm.fit(y = suffStat[, c(x, y)],
                        x = cbind(1, suffStat[, S])))
    x. <- rxy[, 1]
    y. <- rxy[, 2]
  }
  cor.test(x., y.)$p.value
}
my.pc.fit <- pc(suffStat = gmG$x,
                indepTest = myCItest,
                p = ncol(gmG$x),
                alpha = 0.01)
par(mfrow = c(1,2))
plot(gmG$g, main = "")
plot(my.pc.fit, main = "")
dev.off()