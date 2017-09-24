# Bo Ning
# Mar, 6, 2017

# test.data <- test.data
# graph.structure <- graph.structure
# circle <- 7
# iteration <- 20
# stationary = TRUE
# misspecification = FALSE
# plot.figure <- TRUE
# s <- 0.1

estimate.counterfactual <- function(test.data, cntl.index, cntl.data, 
                                    graph.structure, circle = 7,
                                    causal.period, s = 0.1, iteration = 50,
                                    v0.value = seq(1e-6, 0.02, length.out = 5),
                                    stationary = TRUE, 
                                    misspecification = FALSE,
                                    plot.figure = FALSE, plot.title = NULL){
  
  cat("Starting Bayesian EM variable selection... \n")     # report progress
  
  iteration <- iteration
  test.data.non.causal <- test.data[-causal.period, ]
  cntl.data.non.causal <- cntl.data[-causal.period, ]
  source("EMVS.R")
  v0.value <- v0.value
  beta.v0 <- matrix(NA, nc, length(v0.value))
  v0 <- v1 <- theta <- rep(NA, length(v0.value))
  
  if (length(v0.value) == 1) {
    emvs.result <- EMVS(test.data.non.causal, cntl.index, cntl.data.non.causal,
                        graph.structure, circle, v0.value, s, 
                        iteration = iteration, stationary = stationary, 
                        misspecification = misspecification)
    beta.v0 <- emvs.result$beta[, 2]
    theta <- emvs.result$theta
    v1 <- emvs.result$v1
  } else {
    pb  <- txtProgressBar(1, length(v0.value), style=3) # creating progress bar
    
    for (i in 1:length(v0.value)) {
      # report progress
      setTxtProgressBar(pb, i)
   
      emvs.result <- EMVS(test.data.non.causal, cntl.index, cntl.data.non.causal,
                          graph.structure, circle, v0.value[i], s, 
                          iteration = iteration, stationary = stationary, 
                          misspecification = misspecification)
      beta.v0[, i] <- emvs.result$beta[, iteration+1]
      theta[i] <- emvs.result$theta[iteration]
      v1[i] <- emvs.result$v1
    }
  }
  
  # calculate the thresheld
  c <- sqrt(v1/v0.value)
  # beta.threshold <- sqrt( 2*v0.value * log( (1-theta)/theta * c ) * c^2 / (c^2-1))
  
  beta.threshold <- sqrt( (log(v0.value/v1) + 
                             2*log(theta/(1-theta) + 1e-10)) / (1/v1-1/v0.value))
                           # +
                           #   
                           #   2 * log(0.9/0.1)) / (1/v1-1/v0.value))
  dCntl <- sum(cntl.index)
  
  if (plot.figure == TRUE) {
    color <- rep("black", dCntl)
    color[seq(1,dCntl,by = 10)] <- "lightblue"
    color[seq(2,dCntl,by = 10)] <- "blue"
    matplot(v0.value, t(beta.v0), type = "l", col = color, xlab = expression(v[0]),
            ylab = expression(hat(beta)), ylim = c(-2.5, 2.5),
            cex.lab = 1,  mgp=c(2.2,1,0))
    lines(v0.value, beta.threshold, col = "red", lwd = 2)
    lines(v0.value, -beta.threshold, col = "red", lwd = 2)
    title(plot.title)
  }
  
  ###### Deduct the control variable part ######
  beta.star <- beta.threshold[length(v0.value)]
  if (length(v0.value) > 1) {
    beta.hat <- beta.v0[, dim(beta.v0)[2]]
  } else {
    beta.hat <- beta.v0
  }
  beta.hat[abs(beta.hat) < beta.star] <- 0
  cntl.term <- matrix(NA, dim(test.data)[1], dim(test.data)[2])
  index <- 1
  for (i in 1:dim(test.data)[2]) {
    cntl.term[, i] <- cntl.data[, index:(index+9)] %*% beta.hat[index:(index+9)]
    index <- index + 10
  }
  
  return(list(cntl.term = cntl.term, beta.hat = beta.hat))
}