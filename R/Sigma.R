#' Function that computes the covariance matrix of ordinal patterns
#' from a time series, given embedding dimension
#'
#' @param TS time series
#' @param emb embedding dimension
#' @returns A covariance matrix

Sigma <- function(TS, emb){
  
  # Find OP probabilities
  q <- OPprob(TS, emb)
  qpos <- q[q>0]
  
  # Find sum of Q matrices
  k <- length(qpos)
  Q_lag <- matrix(0, nrow = k, ncol = k)
  
  for (l in 1:(emb - 1)){
    Qaux <- Qmatrix(TS, emb, l)
    Q_lag <- Q_lag + Qaux + t(Qaux)
  }
  
  return(diag(qpos) - (2 * emb - 1) * qpos %*% t(qpos) + Q_lag)
}



