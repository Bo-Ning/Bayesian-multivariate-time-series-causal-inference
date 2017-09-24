# This function is used to convert singular matrix to positive definite
# symmetric matrix.

MungeMatrix <- function(v) {
  vScale <- sqrt(abs(diag(v)))
  vScale <- vScale %o% vScale
  vScale[vScale == 0] <- 1e-8
  standardizedV <- v / vScale
  m <- min(eigen(standardizedV)$values)
  if (m < 0) {
    # add 2 abs(m) to all eigenvalues:
    standardizedV <- standardizedV + 2 * (-m) * diag(nrow(v))
    # renormalize
    standardizedV <- standardizedV / (1 + 2 * (-m))
    v <- standardizedV * vScale
  }
  v
}