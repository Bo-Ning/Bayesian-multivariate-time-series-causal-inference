# By Bo Ning
# Mar. 6. 2017
#' \code{two.stage.estimation} the function to calculate KS distance and threshold
#' 
#' @param test.data: T by n dataset with time period T, number of datasets
#'         n
#' @param cntl.index: index to indicate the location of cntl for each test data
#' @param cntl.data: T by p dataset with time perior T, number of total controls
#'                   should put control data for each test data accordingly
#' @param causal.period: Period of dataset has causal impact
#' @param s: the DAEMVS parameter
#' @param emvs.iteration: Iteration for EM algorithm; default is 50
#' @param v0.value: the small value for spike and slab prior
#' @param mcmc.iterloop: total iteration for MCMC; default is 10000
#' @param burnin: burnin for MCMC; default is 2000
#' @param stationary: adding constrain hidden state to be stationary or not
#' @param misspecification: EMVS use the model without time-dependent parameter
#' @param nseasons: seasonality input for analysis
#' @param iterloop: iterations for MCMC
#' @param burnin: burn-in
#' @param graph: impose graph restriction on the variance covariance matrix
#' @param graph.structural: if graph is TRUE, specify graphic structure 
#'        matrix
#' @param num.sim.data: number of counterfactual data to be simulated;
#'                      default is 30
#' @param probs: the percentile for deciding the threshold; default is 0.95
#' @param num.cores: number of cores to run programming for parallel computing,
#'        the default is 1
#' @param plot.EMVS.figure: output EMVS plot
#' @param plot.title: setting EMVS plot title, only useful if ask for EMVS plot

two.stage.estimation <- function(test.data, cntl.index, cntl.data, 
                                 graph = FALSE, graph.structure = FALSE, 
                                 circle, causal.period, s = 0.1,
                                 emvs.iteration = 50, 
                                 v0.value = seq(1e-6, 0.02, length.out = 5),
                                 mcmc.iterloop = 10000, burnin = 2000, 
                                 stationary = TRUE, 
                                 misspecification = FALSE,
                                 num.sim.data = 30, num.cores = 1,
                                 seed = 1, probs = 0.95,
                                 plot.EMVS.figure = FALSE,
                                 plot.title = NULL) {
  
  # load functions
  source("estimate.counterfactual.R")
  source("MultiCausalImpact.R")
  source("MungeMatrix.R")
  source("EMVS.R")
  
  # check if a package has been installed
  pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE)) {
      install.packages(x,dep=TRUE)
      if(!require(x,character.only = TRUE)) stop("Package not found")
    } 
  }
  pkgTest("MASS")
  pkgTest("BDgraph")
  pkgTest("matrixStats")
  pkgTest("parallel")
  pkgTest("matrixcalc")
  pkgTest("MCMCpack")
  
  # load library
  library(MASS)
  library(BDgraph)
  library(matrixStats)
  library(parallel)
  library(matrixcalc)
  
  T <- dim(test.data)[1]
  d <- dim(test.data)[2]
  
  if (T < d) {
    test.data <- t(test.data)
    T <- dim(test.data)[1]
    d <- dim(test.data)[2]
  }
  
  if (graph == TRUE) {
    graph.structure <- matrix(1, d, d)
  } else {
    if (graph.structure == FALSE) {
      stop("Graph structure must provde!")
    }
  }
  
  ########################################################
  ####################### Stage 1 ########################
  ########################################################
  ## Stage 1:
  # EMVS for estimating beta
  # for EMVS, s = 1; for DAEMVS, 0 <= s <= 1
  selection <- 
    estimate.counterfactual(test.data = test.data, cntl.index = cntl.index, 
                            cntl.data = cntl.data, graph.structure = graph.structure, 
                            circle = circle, causal.period = causal.period, s = s, 
                            iteration = emvs.iteration, 
                            v0.value = v0.value,
                            stationary = stationary, plot.figure = plot.EMVS.figure,
                            misspecification = misspecification, 
                            plot.title = plot.title)
  
  cntl.term <- selection$cntl.term
  EMVS.estimator <- selection$beta.hat
  
  ########################################################
  ####################### Stage 2 ########################
  ########################################################
  ## Stage 2:
  # MCMC for time varying parameters and covariance and variance matrices
  # Fit into timeseries model
  model.estimates <- 
    MultiCausalImpact(test.data = test.data, causal.period = causal.period, 
                      cntl.term = cntl.term, seed = seed, nseasons = circle, 
                      iterloop = iterloop, burnin = burnin, 
                      stationary = stationary, graph = graph, 
                      graph.structure = graph.structure,
                      num.sim.data = num.sim.data, probs = probs,
                      num.cores = num.cores)
  
  # collect result
  mcmc.output <- model.estimates$mcmc.output
  threshold <- model.estimates$threshold
  ks.test.cntl <- model.estimates$ks.test.cntl
  
  
  # return results
  return(list(beta.hat = EMVS.estimator, mcmc.output = mcmc.output,
              ks.test.cntl = ks.test.cntl, threshold = threshold))
  
}