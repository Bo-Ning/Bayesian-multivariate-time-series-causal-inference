# clean dataset
rm(list = ls())
setwd("~/BayesianCausalInference/")

# set seed
set.seed(2017030204)
library("MASS")

########################################################
################## Simulate dataset ####################
########################################################
time <- 100 # time length
n <- 5 # number of test
nc <- 50 # number of controls to be generated
n.cntl <- 10 # number of controls to be used for each response

# generate control datasets
cntl.data.pool <- matrix(0, time, nc)
for (i in 1:nc) {
  cntl.data.pool[, i] <- arima.sim(list(ar = c(0.6)), n = time)
}
cntl.data.pool <- cntl.data.pool + abs(min(cntl.data.pool))
# force observations to be positive

# generate test datasets
test.data <- matrix(0, time, n)
A <- matrix(0, n, n) # generate empty matrix, use to create AR(1) correlation
# matrix
graph.structure <- toeplitz(c(1, 1, rep(0, n-2)))
elem <- c(10,5,rep(0,n-2))
Sigma.inv <- toeplitz(elem)
Sigma <- solve(Sigma.inv)
mu <- tau <- matrix(0, n, time)
for (t in 1:time) {
  if (t == 1) {
    mu[, 1] <- 1
  } else {
    mu[, t] <- 0.8*mu[, t-1] + rnorm(n) * 0.1
  }
  test.data[t, ] <- mu[, t] + 
    1*cntl.data.pool[t, seq(1, n*2, by = 2)] + 
    2*cntl.data.pool[t, seq(2, n*2, by = 2)] + 
    mvrnorm(mu = rep(0, n), Sigma = Sigma)
  # add seasonality
  if (t == time) {
    test.data <- test.data + 
      0.1 * cos(2*pi/7*(1:time)) + 0.1 * sin(2*pi/7*(1:time))
  }
}
test.data <- t(t(test.data) + abs(colMins(test.data))+1)
# simulate causal impact
causal.period <- 81:100 # campaign runs 20 periods

# simulate causal impact
for (i in 1:n) {
  test.data[causal.period, i] <- test.data[causal.period, i] + (i-1)*log(1:20)/2
}

########################################################
################# Reorganize dataset ###################
########################################################
a <- 1:nc
index <- 1
cntl.data <- cntl.data.pool[, index:(index+1)]
for (i in 1:n) {
  if (i != 1) {
    cntl.data <- cbind(cntl.data, cntl.data.pool[, index:(index+1)])
  }
  random.sample <- sample(a[-c(1:(2*n))], 8)
  cntl.data.select <- cntl.data.pool[, random.sample]
  cntl.data <- cbind(cntl.data, cntl.data.select)
  index <- index + 2
}

########################################################
##################### Fit model ########################
########################################################
## Stage 1:
# EMVS for estimating beta
iterloop <- 10000
stationary <- TRUE
nseasons <- 7
graph <- TRUE
burnin <- 2000
num.sim.data <- 30
num.cores <- 1
graph.structure <- graph.structure
cntl.index <- rep(10, n) # cntl.index

# load function
source("two.stage.estimation.R")

MultivariateCausalInferenceRes <- 
  two.stage.estimation(test.data, cntl.index, cntl.data, 
                       graph = graph, graph.structure = graph.structure, 
                       circle = nseasons, causal.period = causal.period, 
                       s = 0.1,
                       emvs.iteration = 50, 
                       v0.value = seq(1e-6, 0.02, length.out = 5),
                       mcmc.iterloop = iterloop, burnin = burnin, 
                       stationary = TRUE, 
                       misspecification = FALSE,
                       num.sim.data = 30, num.cores = 1,
                       seed = 1, probs = 0.95,
                       plot.EMVS.figure = FALSE,
                       plot.title = "EMVS plot")

# collect results
beta.hat <- MultivariateCausalInferenceRes$beta.hat
mcmc.output <- MultivariateCausalInferenceRes$mcmc.output
threshold <- MultivariateCausalInferenceRes$threshold
ks.test.cntl <- MultivariateCausalInferenceRes$ks.test.cntl

# save results
write.csv(mcmc.output$prediction.sample, file = "prediction.sample.csv")
write.csv(mcmc.output$mu.sample, file = "mu.sample.csv")
write.csv(mcmc.output$Theta.sample, file= "theta.sample.csv")
write.csv(mcmc.output$D.sample, file = "D.sample.csv")
write.csv(mcmc.output$sigma.sample, file = "sigma.sample.csv")
write.csv(mcmc.output$sigma.U.sample, file = "sigma.U.sample.csv")
write.csv(mcmc.output$sigma.V.sample, file = "sigma.V.sample.csv")
write.csv(mcmc.output$sigma.W.sample, file = "sigma.W.sample.csv")
write.csv(beta.hat, file = "beta.hat")
write.csv(threshold, file = "threshold.csv")
write.csv(ks.test.cntl, file = "ks.test.csv")


