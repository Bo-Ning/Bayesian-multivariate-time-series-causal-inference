# Peter Linton 2015
# Revised by Bo Ning 2016
#
# VARp Sims and fit
# var.bayes takes in a dataset and p integer (order of AR polynomial) and 
# returns a phi estimate of appropriate dimensions

library(mvtnorm)
library(Matrix)

###############################################

sqrtm <- function(A){
  return(eigen(A)$vectors%*%diag(sqrt(eigen(A)$values))%*%t(eigen(A)$vectors))  
}

################################
varp_lkhd<-function(y,phi,sigma,sigma.inv){
  m = nrow(y);n = ncol(y); 
  y1 = y[, 1:(n-1)]
  y2 = y[,2:n] - phi%*%y1  
  gp = matrix(solve((diag(m^2) - phi%x%phi),as.vector(sigma),tol=1e-40),m,m)
  q = sum(diag( solve(gp,as.vector(y[,1]),tol=1e-30)%*%as.vector(y[,1]))) +  
    sum(diag(crossprod(sigma.inv,y2)%*%t(y2) ))
  l = -m*n*log(2*pi)/2 - (log(det(gp)) + (n-1)*log(det(sigma)))/2 - q/2 
  return(l)
}

################################
varppre_lkhd <- function(y,pre,delta,sigma,sigma.inv){
  m = dim(y)[1]
  p = length(delta)
  phi = pre2par_varp(pre,delta)
  l = -varp_lkhd(y,phi,sigma,sigma.inv)
  if (is.nan(l)){l = 1e+5}
  return(l)
}

#############################################

pre2par_varp <- function(pre,delta){
  m = sqrt(length(pre))
  v = matrix(0,m,m)
  q = v
  l = diag(m)
  l[lower.tri(l)] = pre[1:choose(m,2)]
  d = diag(exp(pre[(choose(m,2)+1):choose(m+1,2)]))
  v = l%*%d%*%t(l)
  s= diag(0,m)
  s[lower.tri(s)] = pre[(choose(m+1,2)+1):m^2]
  s = s - t(s)
  q = diag(c(delta,rep(1,(m-1))))%*%(diag(m) - s)%*%solve(diag(m) + s) 
  phi = VQ2par(v,q)
  return(phi)
}

#############################################
# use v, q transform phi
VQ2par <- function(v,q){
  m = dim(v)[1]
  phi = matrix(0,m,m)
  # u1 = diag(m) + v
  u1 = diag(m)*0.9 + v
  u2 = sqrtm(v)%*%q%*%sqrtm(u1)
  txi = u2
  bigu = u1
  phi = txi%*%solve(bigu)
  return(phi)
}

#############################################
pre2VQ <- function(pre,delta){
  m = sqrt(length(pre)) 
  l = diag(m)
  l[lower.tri(l)] = pre[1:choose(m,2)] 
  d = diag(exp(pre[(choose(m,2)+1):choose(m+1,2)]))
  v = l%*%d%*%t(l)
  s = diag(0,m)
  s[lower.tri(s)] = pre[(choose(m+1,2)+1):(m^2)]
  s = s - t(s)
  # q = diag(c(delta,rep(1,(m-1))))%*%(diag(m) - s)%*%solve(diag(m) + s) %*%
  #  ((diag(m) - s)%*%solve(diag(m) + s))
  q <- diag(c(delta,rep(1,(m-1))))%*%(diag(m) - s)%*%solve(diag(m) + s)
  return(list(v = v, q=q))
}

#######################################
par2pre_varp <-function(phi){
  m = dim(phi)[1]
  pre <- rep(0, m^2)
  U = matrix(solve((diag(m^2) - phi%x%phi),as.vector(diag(m)),tol=1e-40),m,m)
  # v = U - diag(m)
  v = U - diag(m)*0.9
  q = solve(sqrtm(v))%*%phi%*%sqrtm(U)
  pre[1:choose(m+1,2)] = V2LDL(v)
  delta = sign(det(q))
  s = 2*solve(diag(m) + diag(c(delta,rep(1,(m-1))))%*%q) - diag(m)
  pre[(choose(m+1,2)+1):m^2] = s[lower.tri(s)]
  return(list(pre=pre,delta=delta,U=U))
}

#######################################

V2LDL <-function(v){
  m = dim(v)[1]
  pre = rep(0,choose(m+1,2))
  c = chol(v, tol = 1e-40)
  d = diag(c)
  l = t(c/d);
  pre[1:choose(m,2)] = l[lower.tri(l)]
  pre[(choose(m,2)+1):choose(m+1,2)] = log(d^2)
  return(pre)
}

##############################################
initial_var <- function(phi, Shrink = 1, Delta = 0.01){
  # Function initializes VARMA(p,q) fit using a VAR(L) and then from the residuals 
  # estimates a VMA(q) to give a VARMA(p,q) estimate
  # Shink allows projection to stationary space, set to 1
  # Inputs: 
  #    y: m x n data matrix
  #    L: integer for VAR(L) fit
  #    p: integer 
  #    q: integer
  #    Shrink: 1 if you want to constrain
  #    delta: small epsilon-style number (distance away from boundary)
  # Outputs:
  #    phi: m x m x p Array
  #    theta: m x m x q Array
  #    Sigma: m x m positive definite covariance matrix
  m <- dim(phi)[1] 
  if(max(abs(eigen(phi)$value))>=1){ 
    phi = phi/((max(abs(eigen(phi)$value)) + Delta))  
  }
  return(phi=phi)
}
