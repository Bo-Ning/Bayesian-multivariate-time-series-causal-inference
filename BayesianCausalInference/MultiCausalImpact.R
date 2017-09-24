# By Bo Ning
# Mar. 6. 2017

#' \code{MultiCausalImpact} the function to calculate KS distance and threshold
#' 
#' @param test.dataset: T by n dataset with time period T, number of datasets
#'         n
#' @param causal.period: Period of dataset has causal impact
#' @param nseasons: seasonality input for analysis
#' @param iterloop: iterations for MCMC
#' @param burnin: burn-in
#' @param stationary: adding constrain hidden state to be stationary or not
#' @param graph: impose graph restriction on the variance covariance matrix
#' @param graph.structural: if graph is TRUE, specify graphic structure 
#'        matrix
#' @param num.sim.data: number of counterfactual data to be simulated
#' @param probs: the percentile for deciding the threshold
#' @param num.cores: number of cores to run programming for parallel computing,
#'        the default is 1
MultiCausalImpact <- function(test.data, causal.period, cntl.term, seed = 1, nseasons = 12, 
                              iterloop = 1000, burnin = 100, stationary = TRUE,
                              graph = FALSE, graph.structure = NULL, num.sim.data = 30,
                              probs = 0.95, num.cores = NA) {
  ############## load packages #################
  library(Matrix)
  library(MCMCpack)
  library(matrixStats)
  library(parallel)
  if (graph == TRUE) {
    library(BDgraph)
  }
  
  ############## load functions ###############
  source("MCMC.multivariate.ssm.R")
  source("koopmanfilter.R")
  
  # set seed
  set.seed(seed)
  
  ############### organize dataset #################
  length <- dim(test.data)[1] # length of dataset
  d <- dim(test.data)[2] # dimension of dataset
  # seperate causal period dataset
  causal.period <- causal.period
  length.non.causal <- length - length(causal.period)
  # Fit deduct cntl.term from test.data
  test.data.tilde <- test.data - cntl.term
  
  #################################################
  # Step 1: Sample draws from posterior distributions of parameters
  #         and obtain the predicted distribution for causal period
  mcmc.model.output <- 
    MCMC.multivariate.ssm(test.data.tilde, causal.period,
                      nseasons = nseasons, iterloop = iterloop, 
                      burnin = burnin, stationary = stationary,
                      graph = graph, graph.structure = graph.structure)
  prediction.sample <- mcmc.model.output$prediction.sample
  a.last.sample <- mcmc.model.output$a.last.sample
  P.last.sample <- mcmc.model.output$P.last.sample
  if (stationary == TRUE) {
    Theta.sample <- mcmc.model.output$Theta.sample
  }
  sigma.sample <- mcmc.model.output$sigma.sample
  sigma.U.sample <- mcmc.model.output$sigma.U.sample
  sigma.V.sample <- mcmc.model.output$sigma.V.sample
  sigma.W.sample <- mcmc.model.output$sigma.W.sample
  D.sample <- mcmc.model.output$D.sample
  z <- mcmc.model.output$z
  R <- mcmc.model.output$R
  trans <- mcmc.model.output$trans
  
  ####################################################
  cat("\nEstimating trend for each simulated counterfactual: \n")
  # report progress
  
  # Step 2: Sample num.sim.data number of counterfactuals from predicted data
  num.sim.data <- num.sim.data # numbers of dataset to simulate
  # generate random numbers
  causal.length <- length(causal.period)
  counterfactual.data <- array(NA, dim = c(causal.length, d, num.sim.data))
  for (t in 1:causal.length) {
    for (dd in 1:d) {
      counterfactual.data[t, dd, ] <- 
        prediction.sample[causal.period[1]+t-1, dd, 
                          sample((burnin+1):iterloop, num.sim.data, replace = T)] 
    }
  }
  
  for (num in 1:num.sim.data) {
    counterfactual.data[, , num] <- 
      counterfactual.data[,,num] + cntl.term[causal.period, ]
  }
  
  ####################################################
  ## Step 3:
  # combine counterfactual dataset and observed dataset
  # and fit them into the model to obtain draws of trend
  combined.data <- abind(counterfactual.data, test.data[causal.period, ], 
                         along = 3)
  
  ####################################################
  ## Step 4: 
  # Fit dataset to calculate trend using parallel computing
  # library(parallel)
  combined.data.estimate.trend <- 
    array(NA, dim = c(causal.length, d, iterloop-burnin, num.sim.data+1))
  
  if (is.na(num.cores) == TRUE) {
    num.cores <- detectCores() - 1
  }
  
  pb  <- txtProgressBar(1, num.sim.data+1, style=3)    # report progress
  
  for (k in 1:(num.sim.data+1)) {
    # report progress
    setTxtProgressBar(pb, k)
    
    data.estimate.trend <- 
      mclapply(
        (burnin+1):iterloop,
        function(x) {
          if (stationary == TRUE) {
            trans[(d+1):(2*d), (d+1):(2*d)] <- Theta.sample[,,x]
          }
          
          # apply kalman-filter and simulation smoother
          alpha.plus <- Matrix(0, length(causal.period), 
                               (min(nseasons, length)+1)*d)
          Q <- bdiag(sigma.U.sample[,,x], sigma.V.sample[,,x],
                     sigma.W.sample[,,x])
          mu.ss <- rep(0, (min(nseasons, length)+1)*d)
          for (t in 1:length(causal.period)) {
            eta <- mvrnorm(1, mu = rep(0, 3*d), Q)
            mu.ss[(d+1):(2*d)] <- (diag(d) - Theta.sample[,,x]) %*% D.sample[,x]
            if (t == 1) {
              alpha.plus[t, ] <- mu.ss + trans %*% a.last.sample[, x] + R %*% eta
            }
            else {
              alpha.plus[t, ] <- mu.ss + trans %*% alpha.plus[t-1, ] + R %*% eta
            }
          }
          data.est.plus <- alpha.plus %*% z + 
            mvrnorm(n = length(causal.period), 
                    mu = rep(0, d), Sigma = sigma.sample[,,x])
          data.est.star <- combined.data[,,k] - data.est.plus 
          # Estimate alpha parameters
          sample.alpha.draws <- 
            koopmanfilter((min(nseasons, length)+1)*d,
                          data.est.star, trans, z, a.last.sample[, x],
                          P.last.sample[,,x], 2*sigma.sample[,,x], 2*Q, R)
          sample.alpha <- sample.alpha.draws + alpha.plus
          return(sample.alpha)
        }, mc.cores=num.cores)
    # convert list object to array
    for (i in 1:(iterloop - burnin)) {
      combined.data.estimate.trend[,,i,k] <- 
        as.matrix(data.estimate.trend[[i]][,1:d])
    }
  }
  
  ####################################################
  # Step 5:
  # compare two distributions:
  # \sum_T+1: T+n \mu_t | Y_obs and \sum_T+1: T+n \mu_t | Y_cf
  combined.data.estimate.culmulate.trend <- 
    array(NA, dim = c(causal.length, d, iterloop-burnin, num.sim.data+1))
  for (t in 1:length(causal.period)) {
    if (t==1) {
      combined.data.estimate.culmulate.trend[t, , , ] <- 
        combined.data.estimate.trend[1,,,]
    } else {
      combined.data.estimate.culmulate.trend[t, , , ] <- 
        apply(combined.data.estimate.trend[1:t,,,], c(2,3,4), mean)
    }}
  
  cat("\nCalculating ks distance...\n") # report progress
  ####################################################
  # Step 6: calculate the threshold for control variables
  ks.cntlsets <- array(NA, dim = c(length(causal.period), 
                                   num.sim.data*(num.sim.data-1), d))
  for (t in 1:length(causal.period)) {
    for (dd in 1:d) {
      a <- 1
      for (i in 1:(num.sim.data)) {
        for (j in 1:(num.sim.data)) {
          if (i != j) {
            ks.cntlsets[t, a, dd] <- 
              ks.test(combined.data.estimate.culmulate.trend[t, dd, , i], 
                      combined.data.estimate.culmulate.trend[t, dd, , j],
                      alternative = "less")$statistic
            a <- a + 1
          }}}}}
  threshold <- apply(ks.cntlsets, c(1,3), quantile, probs = probs)
  
  ####################################################
  # Step 7: calculate the ks distance between control and test variables 
  # Stack control trends given by simulated counterfactual datasets 
  stack.cntl.culmulate.trend <- 
    apply(combined.data.estimate.culmulate.trend[,,,1:num.sim.data], 
          c(1,2,3), mean)
  
  ks.test.cntl <- array(NA, dim = c(length(causal.period), d))
  for (dd in 1:d) {
    for (t in 1:length(causal.period)) {
      ks.test.cntl[t,dd] <- ks.test(
        combined.data.estimate.culmulate.trend[t, dd, , num.sim.data+1],
        stack.cntl.culmulate.trend[t, dd, ],
        alternative = "less")$statistic
    }}
  
  cat("Done! \n")
  # return result
  list(mcmc.output = mcmc.model.output, threshold = threshold,
       ks.test.cntl = ks.test.cntl)
}