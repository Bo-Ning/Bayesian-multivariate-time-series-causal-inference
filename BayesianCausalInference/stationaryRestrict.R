# Causal-Invertible VARMA(p,q) Estimation Code
#
# Peter Linton, Tucker McElroy, Anindya Roy
# 02/19/16
# Revised by Bo Ning
# March, 2016
#
#
# This Code is for fitting a Stationary VAR(1) model for tau_t

stationaryRestrict <- function(y, sigma, sigma.inv = NULL) {
  
  source("varp.R")
  
  # phi <- initial_var(phi)
  
  m <- dim(sigma)[1]
  # lower triangle have priors N(0, 5), diagonal elements have priors N(0, 3)
  #priorsig <- rep(c(rep(5,m*(m-1)/2), rep(3,m), rep(5,m*(m-1)/2)), p)
  priorsig <- c(rep(5,m*(m-1)/2), rep(5,m), rep(5,m*(m-1)/2))
  jumpsig = rep(0.05, m^2);
  priorlogor = 0; jump = .5
  y <- t(y)
  length <- dim(y)[2]
  # length <- dim(y)[1]
  if (is.null(sigma.inv) == TRUE) {
    sigma.inv <- solve(sigma)
  } 
  # find ols estimator of phiyw
  phi <- matrix(cbind(ar.ols(as.matrix(t(y)), aic = FALSE, order.max = 1,
                            demean = F, intercept = F)$ar), m, m)

  # if the not p.d, then constraint to p.d.
  phi <- initial_var(phi, Shrink = 1)
  preyw = par2pre_varp(phi)
  pre = preyw$pre
  # cat(pre, "\n")
  delta = preyw$delta
  U <- preyw$U
  # 
  if (det(U-diag(m)) <= 1e-10) {
     phi <- phi
   } else {
   
    length.v <- m*(m+1)/2
    prenew = pre
    prenew[1:length.v] = pre[1:length.v] + rnorm(length.v,0,jumpsig[1:length.v])
    logr = -.5*sum((prenew[1:length.v]^2-pre[1:length.v]^2)/(priorsig[1:length.v]^2)) -
      varppre_lkhd(y,prenew,delta,sigma,sigma.inv) +
      varppre_lkhd(y,pre,delta,sigma,sigma.inv)
    if(log(runif(1)) < logr){pre[1:length.v] = prenew[1:length.v]}

    # update q, given v and delta
    prenew[(length.v+1):(m^2)] = pre[(length.v+1):(m^2)] +
      rnorm(m^2-length.v,0,jumpsig[(length.v+1):(m^2)])
    logr = -.5*sum((prenew[(length.v+1):(m^2)]^2-
                      pre[(length.v+1):(m^2)]^2)/(priorsig[(length.v+1):(m^2)]^2)) -
      varppre_lkhd(y,prenew,delta,sigma,sigma.inv) +
      varppre_lkhd(y,pre,delta,sigma,sigma.inv)
    if(log(runif(1)) < logr){pre[(length.v+1):(m^2)] = prenew[(length.v+1):(m^2)]}
    
    # update delta given v, q
    deltanew = delta
    deltanew = 2*rbinom(1, 1, jump)-1
    logr = sum(sign(deltanew-delta)*priorlogor) -
      varppre_lkhd(y, pre, deltanew, sigma, sigma.inv) +
      varppre_lkhd(y, pre, delta, sigma, sigma.inv)
    if(log(runif(1))< logr){delta = deltanew}
    phi = pre2par_varp(pre, delta)
   }

  return(phi)
}
