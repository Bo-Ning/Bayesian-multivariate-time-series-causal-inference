# Bo Ning
# Mar, 6, 2017

### KALMAN-FILTER FORWARD
#' \code{kalmflter} kalman-filter for forward draws of alphas
#'
#' @param y observed data
#' @param trans coefficient matrix of state equation
#' @param z coefficient vector/matrix of alpha
#' @param a starting point of mean of alpha
#' @param P starting point of variance of alpha
#' @param sigma.epsln variance of state equation
#' @param Q covariance-variance matrix of hidden process
#' @param R rank deficient matrix 
#' @param missing let kalman-filter do data imputation if missing is true
#' @return updated a, updated P, updated alpha
#' 
#' @seealso
kalmflter <- function(y, trans, z, a, P, sigma, Q, R, missing = FALSE) {
  if (missing == TRUE) {
    a.update <- trans %*% a
    P.update <- trans %*% P %*% t(trans) + R %*% Q %*% t(R)
    list(a = a.update, P = P.update)
  } else {
    v <- y - crossprod(z, a)  
    FF <- t(z) %*% P %*% z + sigma
    if (is.singular.matrix(as.matrix(FF)) == TRUE) {
      FF <- MungeMatrix(FF)
    }
    # print(det(FF))
    K <- trans %*% P %*% z %*% solve(as.matrix(FF))
    L <- trans - K %*% t(z)
    a.update <- trans %*% a + K %*% v
    P.update <- trans %*% P %*% t(L) + R %*% Q %*% t(R)
    # return updated values
    list(v = v, FF = FF, L = L, a = a.update, P = P.update)
  }
}