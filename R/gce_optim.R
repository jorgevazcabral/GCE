# Functions for optimization

#' Objective function for primal optimization formulation using solnp
#'
#' Returns entropy value
#'
#' @param x0 Initial values for the parameters to be optimized over.
#' @param X .
#' @param n .
#' @param k .
#' @param m .
#' @param j .
#' @param s1 .
#' @param S .
#'
#' @return Entropy value
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

ObjFunGCE.primal.solnp <- function(x0,X,n,k,m,j,q1,q2,s1,S) {
  # p <- abs(matrix(x0[1:(k * m)], k, m))
  # w3 <- abs(matrix(x0[(k * m + 1):length(x0)], n, j))
  # return(sum(p * log(p)) + sum(w3 * log(w3)))
  p <- as.matrix(as.numeric(abs(matrix(x0[1:(k * m)], k, m))),
                 ncol = 1)
  w3 <- as.matrix(as.numeric(abs(matrix(x0[(k * m + 1):length(x0)], n, j))),
                  ncol = 1)
  q1rep <- rep(q1, k)
  q2rep <- rep(q2, n)

  return(as.numeric((t(p) %*% log(p) - t(p) %*% log(q1rep)) + (t(w3) %*% log(w3) - t(w3) %*% log(q2rep))))
}

#' Constraint function for primal optimization formulation using solnp
#'
#' Returns y value and probabilities
#'
#' @inheritParams ObjFunGCE.primal.solnp
#'
#' @return Entropy value
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

ConstFunGCE.primal.solnp <- function(x0,X,n,k,m,j,q1,q2,s1,S) {
  aux.x0.X <- matrix(x0[1:(k * m)], k, m)
  aux.x0.error <- matrix(x0[(k * m + 1):length(x0)], n, j)
  aux.y = as.matrix(X) %*% (rowSums(s1 * aux.x0.X)) +
    rowSums(S * aux.x0.error)
  aux.p = c(rowSums(aux.x0.X), rowSums(aux.x0.error))
  return(c(aux.y, aux.p))
}

#' Objective function for primal optimization formulation using solnl
#'
#' Returns entropy value
#'
#' @param x0 Initial values for the parameters to be optimized over.
#' @param dp .
#' @param k .
#' @param q1 .
#' @param n .
#' @param q2 .
#'
#' @return Entropy value
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

ObjFunGCE.primal.solnl <- function(x0, dp, k, q1, n, q2) {
  p <- as.matrix(x0[1:dp])
  w <- as.matrix(x0[(dp + 1):ncol(x0)])

  q1rep <- rep(q1, k)
  q2rep <- rep(q2, n)

  (t(p) %*% log(p) - t(p) %*% log(q1rep)) + (t(w) %*% log(w) - t(w) %*% log(q2rep))
}
