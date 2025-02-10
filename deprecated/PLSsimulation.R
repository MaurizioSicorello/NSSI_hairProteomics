##############################################
# Simulation of PLS data with two highly correlated clusters
# Different values might hint that significance tests of loadings might often fail, because  
#  "component 1" might be "component 2" in other scenarios when they explain similar variance in the true model


# packages
library("MASS")
library("caret")

# simulation parameters
N = 68
nVars = 15
bootIt = 1000
noise = 0.3

predCovMat <- matrix(c(1, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0.5, 1, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0.5, 0.5, 1, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0.5, 0.5, 0.5, 1, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0.5, 0.5, 0.5, 0.5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 1, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0.5, 1, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0.5, 0.5, 1, 0.5, 0.5, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 1, 0.5, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 1, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1), 
                     nrow = 15, ncol = 15, byrow = TRUE)



# simulate data

data <- mvrnorm(n = N, mu = rep(0, nVars), Sigma = predCovMat)

comp1_weights <- c(rep(0.3, 3), rep(0, 15-3))
comp2_weights <- c(rep(0, 5), rep(0.3, 3), rep(0, 15-3-5))
comp3_weights <- c(rep(0, 10), rep(0.3, 3), rep(0, 15-10-3))

comp1 <- data %*% comp1_weights + rnorm(N, 0, noise)
comp2 <- data %*% comp2_weights + rnorm(N, 0, noise)
comp3 <- data %*% comp3_weights + rnorm(N, 0, noise)

outcome <- 0.3*comp1 + 0.3*comp2 + 0.3*comp3 + rnorm(N, 0, noise)

data <- data.frame(data)
names(data) <- paste0("x", c(1:nVars))

data$outcome <- outcome


# perform PLS

bootWeightOut <- data.frame(matrix(
      nrow = bootIt, ncol = (ncol(data)-1),
      dimnames = list(NULL, names((data))[-16])
    ))

for(i in 1:bootIt){
  
  dfBootWeights <- data[sample(x = c(1:nrow(data)), replace = TRUE), ]
  
  plsModel_boot <- plsr(data = dfBootWeights, outcome ~ ., ncomp = 3)
  
  bootWeightOut[i, ] <- loadings(plsModel_boot)[, 1]
  
}

boot_p_mat(bootWeightOut)

hist(bootWeightOut$x2)

plsModel <- plsr(data = data, outcome ~ ., ncomp = 3)
imp <- varImp(plsModel)$Overall

dfPlot <- data.frame("varImportance" = imp, "varNames" = factor(paste0("x", 1:15), levels = paste0("x", 1:15)))

ggplot(data = dfPlot, aes(x = varNames, y = varImportance)) +
  geom_bar(stat = "identity")





