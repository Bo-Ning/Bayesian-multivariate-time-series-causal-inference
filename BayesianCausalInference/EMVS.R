# By Bo Ning
# Mar. 06. 2017
# 
# test <- test.data[-causal.period, ]
# cntl <- cntl.data[-causal.period, ]
# cntl.index <- cntl.index
# graph.structure <- graph.structure
# circle <- 7
# v0 <- v0.value <- seq(1e-6, 0.02, length.out = 10)[1]
# s <- 1
# iterloop <- 50
# stationary <- T

EMVS <- function(test, cntl.index, cntl, graph.structure, circle,
                 v0 = 0.1, s = 0.1, iteration = 50, stationary = TRUE, 
                 misspecification = FALSE) {
  
  library(abind)
  source("koopmanfilter.R")
  
  ##################### EMVS ########################
  # organize dataset
  length <- dim(test)[1]
  n <- dim(test)[2]
  dCntl <- sum(cntl.index)
  
  cntl.input <- cntl
  # re-organize cntl dataset
  cntl <- NULL
  index <- 0
  for (dim in 1:n) {
    cntl <- abind(cntl, cntl.input[, (index+1):(index+cntl.index[dim])], 
                  along = 3)
    index <- index+cntl.index[dim]
  }
  
  ############################################################
  # rename y.est to y.matrix and cntl.est to x.matrix
  x.std <- array(0, c(length,cntl.index[1],n))
  for (i in 1:(dim(test)[2])) {
    x.std[,,i] <- t(t(cntl[,,i]) - colMeans(cntl[,,i]))
  }
  x.matrix <- matrix(0, n*length, dCntl)
  index <- 0
  for (dims in 1:n) {
    x.matrix[seq(dims,n*length,by=n), 
             (index+1):(index+cntl.index[dims])] <- x.std[, , dims]
    index <- index+cntl.index[dims]
  }
  y.std <- t((t(test) - colMeans(test)))
  y.matrix <- as.matrix(c(t(y.std)))
  
  # giving starting values for beta and sigma
  beta.hat <- rep(0, dCntl)
  y.tilde <- matrix(y.matrix - x.matrix %*% beta.hat, n, length)
  sigma.hat <- crossprod(y.std) / (length-1)
  # if sigma.hat is singular, we convert into a p.d. matrix
  if (is.singular.matrix(sigma.hat) == TRUE){
    sigma.hat <- MungeMatrix(sigma.hat)
  }
  sigma.hat.inv <- solve(sigma.hat)
  sigma.hat.inv[graph.structure == 0] <- 0
  B <- diag(n) # prior parameters for wishart priors
  delta <- n + 1 # prior parameter for wishart priors
  
  # giving staring points for v0, v1, a, b, theta
  v0 <- v0 # initial v0
  # v1 <- 100 # initial v1
  v1 <- 10
  a <- 1 # intial a, shape1 in beta distribution
  b <- 1  # intial b, shape2 in beta distribution
  theta <- rbeta(1, a, b) # intial theta from beta distribtuion
  
  if (misspecification == FALSE) {
    # giving starting values for Q
    k1 <- k2 <- k3 <- 0.1 # prior parameters for sigma.u, sigma.v, sigma.w
    sigma.u.hat <- chol2inv(rgwish(n = 1, adj = graph.structure,
                                   b = n+1, D = k1^2 * n * diag(n)))
    sigma.v.inv <- rgwish(n = 1, adj = graph.structure,
                          b = n+1, D = k2^2 * n * diag(n))
    sigma.v.hat <- chol2inv(sigma.v.inv)
    sigma.v.hat.inv <- solve(sigma.v.hat)
    sigma.w.hat <- chol2inv(rgwish(n = 1, adj = graph.structure,
                                   b = n+1, D = k3^2 * n * diag(n)))
    Q.hat <- bdiag(sigma.u.hat, sigma.v.hat, sigma.w.hat)
    Q.inv <- solve(Q.hat)
  
    # giving initial parameters for KFBS (Kalman filter and backward simulation)
    # initial a.int
    a.int <- rep(0, n*(circle+1)) 
    # initialize P.int
    P.int <- Matrix(0, n*(circle+1), n*(circle+1))
    
    if (stationary == TRUE) {
      P.int[1:(3*n), 1:(3*n)] <- diag(3*n)
    } else {
      P.int[1:(3*n), 1:(3*n)] <- diag(3*n)*10^6
    }
  
    # initial transition matrix
    trans <- Matrix(0, (circle+1)*n, (circle+1)*n)
    linear <- diag(2*n)
    linear[1:n, (n+1):(2*n)] <- diag(n)
    trans[1:(2*n), 1:(2*n)] <- linear
    # take initial variance of tau from the data
    if (stationary == TRUE) {
      data.yw <- ar.yw(as.matrix(test), aic = FALSE, order.max = 1,
                       demean = T, intercept = T)
      phi.hat <- matrix(cbind(data.yw$ar), n, n)
    } else {
      phi.hat <- diag(n)
    }
    trans[(n+1):(2*n), (n+1):(2*n)] <- phi.hat
    seasonal <- Matrix(0, circle-1, circle-1)
    seasonal[1, ] <- -1
    seasonal[2:(circle-1), 1:(circle-2)] <- diag(circle-2)
    for (dims in 1:n) {
      trans[seq((2*n+dims), (circle+1)*n, by=n), 
            seq((2*n+dims), (circle+1)*n, by=n)] <- seasonal
    }
  
    
    # define
    z <- Matrix(0, n*(circle+1), n) 
    z[1:n, ] <- diag(n)
    z[(2*n+1):(3*n), ] <- diag(n)
    
    # initialize R 
    R <- Matrix(0, n*(circle+1), n*3)
    R[1:(3*n), 1:(3*n)] <- diag(3*n)
    
    # initialize alpha.hat
    a.hat <- matrix(0, n*(circle+1), length)
  }
  
  # step 2: EM update parameters
  # create matrix to collect results
  iterloop <- iteration
  beta.update <- matrix(0, dCntl, iterloop+1)
  sigma.update <- array(NA, c(n,n,iterloop)) 
  v0.update <- v1.update <- theta.update <- rep(NA, iterloop)
  lp.update <- rep(NA, iterloop)
  if (misspecification == FALSE) {
    a.update <- array(NA, c(length, n*(circle+1),iterloop))
    phi.update <- array(NA, c(n,n,iterloop))
    sigma.u.update <- sigma.v.update <- sigma.w.update <- 
      array(NA, c(n,n,iterloop)) 
  }
  
  # begin for-loop
  for (iter in 1:iterloop) {

    # ----------------- E-step ------------------- #
    
    if (misspecification == FALSE) {
      # upsing kalman filter and backward smoother
      ## update alpha
      KFBS <- koopmanfilter(n*(circle+1), t(y.tilde), trans, z, a.int, P.int, 
                            sigma.hat, Q.hat, R, output.var.cov = TRUE)
      a.hat <- KFBS$a.sample
      P.hat <- KFBS$P.sample
      P.cov.hat <- KFBS$P.cov.sample
      
      # calculate expectation for E(alpha_t alpha_t') and E(alpha_t alpha_(t-1)')
      V.hat <- array(0, c(n*(circle+1),n*(circle+1),length))
      V.cov.hat <- array(0, c(n*(circle+1),n*(circle+1),length-1))
      for (i in 1:length) {
        V.hat[,,i] <- P.hat[,,i] + tcrossprod(a.hat[,i])
        if (i < length)
          V.cov.hat[,,i] <- P.cov.hat[,,i] + tcrossprod(a.hat[,i], a.hat[,i+1])
      }
      
    }
    # ----------------- M-step ------------------- #
    ## update A
    gamma1 <- dnorm(beta.hat, mean = 0, sd = sqrt(v1))
    gamma2 <- dnorm(beta.hat, mean = 0, sd = sqrt(v0))
    pstar <- (gamma1*theta)^s / ((gamma1*theta)^s + (gamma2*(1-theta))^s)
    A <- diag(dCntl)
    diag(A) <- (1-pstar)/v0 + pstar/v1
    
    ## update beta
    # organize dataset
    kron.sigma.inverse <- kronecker(diag(length), sigma.hat.inv)
    XcovX <- crossprod(x.matrix, kron.sigma.inverse) %*% x.matrix
    
    if (misspecification == FALSE) {
      a.z <- c(as.matrix(t(t(a.hat) %*% z)))
      XcovY <- crossprod(x.matrix, kron.sigma.inverse) %*% 
        (y.matrix - a.z)
    } else {
      XcovY <- crossprod(x.matrix, kron.sigma.inverse) %*% 
        (y.matrix)
    }
    
    beta.hat <- as.vector(solve(XcovX + A, XcovY))
    
    gamma1 <- dnorm(beta.hat, mean = 0, sd = sqrt(v1))
    gamma2 <- dnorm(beta.hat, mean = 0, sd = sqrt(v0))
    pstar <- (gamma1*theta)^s / ((gamma1*theta)^s + (gamma2*(1-theta))^s)
    
    # update theta
    theta <- (sum(pstar) + a - 1) / (a + b + dCntl - 2)
    
    if (misspecification == FALSE) {
      # update phi
      # vec(phi) ~ N(0, 0.01*I_{n^2}) prior of vec(phi)
      if (stationary == TRUE) {
        phi.term1 <- phi.term2 <- 0
        for (i in 1:(length-1)) {
          
          phi.term1 <- phi.term1 +
            kronecker(t(V.hat[(n+1):(2*n), (n+1):(2*n), i]), sigma.v.hat.inv)
          phi.term2 <- phi.term2 +
            kronecker(t(V.cov.hat[(n+1):(2*n),(n+1):(2*n), i]), sigma.v.hat.inv)
          
        }
        
        vec.phi <- solve(phi.term1 + 10*diag(n^2), phi.term2 %*% c(diag(n)))
        phi.hat <- matrix(vec.phi, n, n)
      } else {
        phi.hat <- diag(n)
      }
      
      # update trans
      trans[(n+1):(2*n), (n+1):(2*n)] <- phi.hat
      
      # update sigma
      y.tilde <- matrix(y.matrix - x.matrix %*% beta.hat, n, length)
      P.plus.aa <- 0
      for (i in 1:length) {
        P.plus.aa <- P.plus.aa + V.hat[,,i] 
      }
      sigma.mat <- y.tilde %*% t(y.tilde) - t(z) %*% a.hat %*% t(y.tilde) -
        y.tilde %*% t(a.hat) %*% z + t(z) %*% P.plus.aa %*% z
      sigma.hat <- (sigma.mat + B) / (length + delta - 2) 
      sigma.hat.inv <- solve(sigma.hat)
      #sigma.hat.inv[graph.structure == 0] <- 0
      if (min(abs(eigen(sigma.hat.inv)$values)) <= 0) {
        sigma.hat.inv <- MungeMatrix(sigma.hat.inv)
      }
      sigma.hat <- solve(sigma.hat.inv)
      
      # update Q
      Q.mat.term <- 0
      for (i in 1:(length-1)) {
        Q.mat.term <- Q.mat.term +
          V.hat[,,i+1] - trans %*% V.cov.hat[,,i] - 
          t(V.cov.hat[,,i]) %*% t(trans) + trans %*% V.hat[,,i] %*% t(trans)
      }
      Q.mat <- t(R) %*% Q.mat.term %*% R
      Q.hat <- (Q.mat + bdiag((n+1)*k1^2*B, (n+1)*k2^2*B, (n+1)*k3^2*B)) /
        (length + delta - 3)
      
      # update sigma.u
      sigma.u.hat <- Q.hat[1:n, 1:n]
      sigma.u.hat.inv <- solve(sigma.u.hat)
      # sigma.u.hat.inv[graph.structure == 0] <- 0
      if (min(abs(eigen(sigma.u.hat.inv)$values)) <= 0) {
        sigma.u.hat.inv <- MungeMatrix(sigma.v.hat.inv)
      }
      
      # update sigma.v
      sigma.v.hat <- Q.hat[(n+1):(2*n), (n+1):(2*n)]
      sigma.v.hat.inv <- solve(sigma.v.hat)
      # sigma.v.hat.inv[graph.structure == 0] <- 0
      if (min(abs(eigen(sigma.v.hat.inv)$values)) <= 0) {
        sigma.v.hat.inv <- MungeMatrix(sigma.v.hat.inv)
      }
      
      # update sigma.w
      sigma.w.hat.inv <- solve(Q.hat[(2*n+1):(3*n), (2*n+1):(3*n)])
      # sigma.w.hat.inv[graph.structure == 0] <- 0
      if (min(abs(eigen(sigma.w.hat.inv)$values)) <= 0) {
        sigma.w.hat.inv <- MungeMatrix(sigma.w.hat.inv)
      }
      
      # update Q.hat
      Q.inv <- bdiag(sigma.u.hat.inv, sigma.v.hat.inv, sigma.w.hat.inv)
      Q.hat <- solve(Q.inv)
    } else {
      
      # ------------- FOR MISSPECIFIED MODEL ---------------- #
      # update sigma
      y.tilde <- matrix(y.matrix - x.matrix %*% beta.hat, n, length)
      sigma.mat <- y.tilde %*% t(y.tilde) 
      sigma.hat <- (sigma.mat + B) / (length + delta - 2) 
      sigma.hat.inv <- solve(sigma.hat)
      #sigma.hat.inv[graph.structure == 0] <- 0
      if (min(abs(eigen(sigma.hat.inv)$values)) <= 0) {
        sigma.hat.inv <- MungeMatrix(sigma.hat.inv)
      }
      sigma.hat <- solve(sigma.hat.inv)
    }
    
    # collect result
    if (misspecification == FALSE) {
      a.update[,, iter] <- a.hat
      phi.update[,,iter] <- phi.hat
    }
    beta.update[, iter+1] <- beta.hat
    sigma.update[, , iter] <- as.matrix(sigma.hat)
    theta.update[iter] <- theta
  }
 
  return(list(beta = beta.update, theta = theta.update, v1 = v1))
   
}
