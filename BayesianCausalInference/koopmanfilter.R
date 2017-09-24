# By Bo Ning
# Mar. 6. 2017

### KOOPMANFILTER
#' \code{koopmanfilter} sample alpha from posterior distribution
#'
#' @param T time length 
#' @param n dimension of state equation
#' @param y observed data
#' @param trans coefficient matrix of state equation
#' @param z coefficient vector/matrix of hidden process
#' @param a starting point of mean of alpha
#' @param P starting point of variance of alpha
#' @param intcp intercept of hidden state 
#' @param sigma.epsln variance of state equation
#' @param Q covariance-variance matrix of hidden process
#' @param R rank deficient matrix 
#' @param output.mean.var.cov ask to give mean a_t, variance p_t and covariance 
#'                            p_{t, t-1}. This is used for EMVS algorithm
#' @return alpha.sigma
#' 
#' @seealso 
koopmanfilter <- function(n, y, trans, z, a.int, P.int,
                          sigma, Q, R, causal.period = NULL,
                          output.var.cov = FALSE) {
  
  # load function
  source("kalmflter.R")
  
  # get data length and dimension
  T <- dim(y)[1]
  d <- dim(y)[2]
  
  # create empty matrix to collect data
  alpha.sample <- matrix(0, n, T)
  v.sample <- matrix(0, T, d)
  F.sample <- array(0, c(d, d, T))
  L.sample <- array(0, c(T, n, n))
  alpha.plus <- matrix(0, T, n)
  r <- matrix(0, n, T+1)
  N <- array(0, c(n, n, T+1))
  a.sample <- matrix(0, n, T)
  a.ff <- matrix(0, n, T+1)
  a.ff[,1] <- as.vector(a.int)
  P.sample <- array(0, c(n, n, T))
  P.ff <- array(0, c(n, n, T+1))
  P.ff[,,1] <- as.matrix(P.int)
  P.cov.sample <- array(0, c(n, n, T-1))
  
  # kalman-filter
  if (is.null(causal.period) == TRUE) {
    a <- a.int
    P <- P.int 
    for (t in 1:T) {
      kalmanfilter <- kalmflter(y[t, ], trans, z, a, P, sigma, Q, R)
      v.sample[t, ] <- as.vector(kalmanfilter$v)
      F.sample[, , t] <- as.matrix(kalmanfilter$FF)
      L.sample[t, , ] <- as.matrix(kalmanfilter$L)
      a <- kalmanfilter$a
      P <- kalmanfilter$P
      a.ff[,t+1] <- as.vector(a)
      P.ff[,,t+1] <- as.matrix(P + t(P)) / 2 # force to be symmetric matrix
    }
  } else {
    a <- a.int
    P <- P.int
    for (t in 1:T) {
      if (t %in% causal.period) {
        kalmanfilter <- kalmflter(y[t, ], trans, z, a, P, sigma, Q, R, 
                                  missing = TRUE)
        a <- kalmanfilter$a
        P <- kalmanfilter$P
        a.ff[,t+1] <- as.vector(a)
        P.ff[,,t+1] <- as.matrix(P)
      } else {
        kalmanfilter <- kalmflter(y[t, ], trans, z, a, P, sigma, Q, R)
        v.sample[t, ] <- as.matrix(kalmanfilter$v)
        F.sample[, , t] <- as.matrix(kalmanfilter$FF)
        L.sample[t, , ] <- as.matrix(kalmanfilter$L)
        a <- kalmanfilter$a
        P <- kalmanfilter$P
        a.ff[,t+1] <- as.vector(a)
        P.ff[,,t+1] <- as.matrix(P + t(P)) / 2
        if (t == (causal.period[1]-1)) {
          a.last <- a
          P.last <- P
        }
      }
    }
  }
  
  # ------------ BACKWARD RECURSION ---------------- #
  
  # backward recursion to obtain draws of r and N
  for (t in T:1) {
    if (t %in% causal.period) {
      r[, t] <- as.vector(r[, t+1] %*% trans)
      N[,,t] <- as.matrix(t(trans) %*% N[,, t+1] %*% trans)
    } else {
      r[, t] <- as.matrix(z %*% solve(F.sample[, , t]) %*% v.sample[t, ] + 
                         t(L.sample[t, , ]) %*% r[, t+1])
      N[,,t] <- as.matrix(z %*% solve(F.sample[,,t]) %*% t(z)) + 
        as.matrix(t(L.sample[t,,]) %*% N[,,t+1] %*% L.sample[t,,])
    }
  }
  
  # obtain draws of alpha
  alpha.sample[, 1] <- as.matrix(a.int + P.int %*% r[, 1])
  a.sample[, 1] <- as.matrix(a.int + P.int %*% r[,1])
  P.sample[,,1] <- as.matrix(P.ff[,,1]) - 
    as.matrix(P.ff[,,1] %*% N[,,1] %*% P.ff[,,1])
  for (t in 2:T) {
    alpha.sample[, t] <- as.matrix(trans %*% alpha.sample[, t-1] + 
                                  R %*% Q %*% t(R) %*% r[, t])
    a.sample[,t] <- as.matrix(a.ff[, t] + P.ff[,,t] %*% r[, t])
    P.sample[,,t] <- as.matrix(P.ff[,,t]) - 
      as.matrix(P.ff[,,t] %*% N[,,t] %*% P.ff[,,t])
    P.sample[,,t] <- (P.sample[,,t] + t(P.sample[,,t])) / 2 
    P.cov.sample[,,t-1] <- as.matrix(P.ff[,,t-1] %*% t(L.sample[t-1,,])) %*% 
      as.matrix((diag(n) - N[,,t] %*% P.ff[,,t]))
  }
  alpha.sample <- t(alpha.sample)
  a.sample <- a.sample
  P.sample <- P.sample
  P.cov.sample <- P.cov.sample
  
  if (output.var.cov == TRUE) {
    list(alpha.sample = alpha.sample, a.sample = a.sample, P.sample = P.sample, 
         P.cov.sample = P.cov.sample)
  } else {
    # return value
    if (is.null(causal.period) == TRUE) {
      return(alpha.sample = alpha.sample)
    } else {
      list(alpha.sample = alpha.sample, a.last = a.last, P.last = P.last) 
    }
  }
}