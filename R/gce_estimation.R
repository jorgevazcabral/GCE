#' Generalized Cross entropy estimation
#'
#' Internal function used to fit a linear regression model via generalized cross
#' entropy where initial support spaces can be provided or computed.
#'
#' @param y A vector of observations of length n, or a matrix with n rows.
#' @param X A design matrix of dimension n * k.
#' @param y.test A vector of observations of length n, or a matrix with n rows
#'  to be used to determine the out-of sample error. If \code{NULL} the error
#'  will be computed using \code{y}. The default is \code{y.test = NULL}.
#' @param X.test A design matrix of dimension n * k containing the same
#'  variables as \code{X} to be used to determine the out-of sample error.
#'  If \code{NULL} the error will be computed using \code{X}. The default is
#'  \code{X.test = NULL}.
#' @param residuals Boolean value. If \code{TRUE} residuals are returned.  The
#' default is \code{residuals = TRUE}.
#' @param fittedvalues Boolean value. If \code{TRUE} fitted values are returned.
#' The default is \code{fittedvalues = TRUE}.
#' @param int.one NULL or fixed positive upper limit for the supports on
#' standardized data; a pair (-a,a), with a positive, or a matrix ((k+1) x 2)
#' for the supports on original data. The default is \code{int.one = NULL}.
#' @param int.one.q1 A vector of prior weights for the signal. The sum of weights
#'  must be equal to 1. The length of \code{int.one.q1} defines M. The default
#'  is \code{int.one.q1 = c(1 / 3, 1 / 3, 1 / 3)}.
#' @param int.one.errormeasure Loss function to use for the selection of support.
#'  One of c("RMSE","MAE","MSE"). The default is
#'  \code{int.one.errormeasure = "RMSE"}.
#' @param int.two An interval centered around zero given in the form c(-a,a),
#' a positive. If NULL then a is calculated using the empirical three-sigma rule
#' Pukelsheim (1994). The default is \code{int.two = NULL}.
#' @param int.two.q2 A vector of prior weights for the noise. The sum of weights
#' must be equal to 1. The length of \code{int.two.q2} defines J. The default
#'  is \code{int.two.q2 = c(1 / 3, 1 / 3, 1 / 3)}.
#' @param method  Use \code{"primal.solnl"} (GCE using SQP) or
#'  \code{"primal.solnp"} (GME) for primal form of the optimization problem and
#'  \code{"dual"} (GME) for dual form. The default is \code{method = "dual"}.
#' @param method.maxfeval Maximum number of function evaluations for
#' \code{"primal.solnl"}. The default is \code{method.maxfeval = 1e+04}.
#' @param method.maxit Maximum number of minor (inner) iterations for
#'  \code{"primal.solnp"} or maximum number of iterations for
#'  \code{"primal.solnl"} and \code{"dual"}. The default is
#'  \code{method.maxit = 1e+05}.
#' @param method.tol Relative tolerance on feasibility and optimality for
#' \code{"method"}. The default is \code{method.tol = 1e-06}.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd
gce_estimation <- function(y,
                           X,
                           y.test = NULL,
                           X.test = NULL,
                           residuals = TRUE,
                           fittedvalues = TRUE,
                           int.one = NULL,
                           int.one.q1 = c(1 / 3, 1 / 3, 1 / 3),
                           int.one.errormeasure = "RMSE",
                           int.two = NULL,
                           int.two.q2 = c(1 / 3, 1 / 3, 1 / 3),
                           method = "dual",
                           method.maxfeval = 1e+04,
                           method.maxiter = 1e+05,
                           method.tol = 1e-06) {
  X_original <- X
  y_original <- y

  all.error.measures <-
    matrix(0, 1, 3, dimnames = list("value", c("MSE", "RMSE", "MAE")))

  n <- nrow(X)
  k <- ncol(X)

  m <- length(int.one.q1)
  j <- length(int.two.q2)


  ### Carefull with standardization

    if (any(c("(Intercept)", "X.Intercept.") %in% colnames(X_original))) {
      X <- scale(X[, -1])
      k <- k - 1
    } else {
      X <- scale(X)
    }

    y <- scale(y)

    intg <- matrix(c(-int.one, int.one), k, 2, byrow = length(int.one) == 1)

  if (is.null(int.two)) {
    int2 <- c(floor(-3 * sd(y[, 1])), ceiling(3 * sd(y[, 1])))
  } else {
    int2 <- int.two
  }

  s1 <- matrix(0, k, m)
  Z <- matrix(0, k, k * m)
  for (i in 1:k) {
    s1[i,] <- seq(intg[i, 1],
                  intg[i, 2],
                  by = (intg[i, 2] - intg[i, 1]) / (m - 1))
    Z[i, ((i - 1) * m + 1):(i * m)] <- s1[i, 1:m]
  }

  s2 <- seq(int2[1], int2[2], by = (int2[2] - int2[1]) / (j - 1))
  S <- matrix(rep(s2, n), ncol = length(s2), byrow = TRUE)
  V <- matrix(0, n, n * j)
  for (i in 1:n) {
    V[i, ((i - 1) * j + 1):(i * j)] <- s2
  }

  if (method == "primal.solnp") {
    res.opt = Rsolnp::solnp(
      pars = c(rep(1 / m, k * m), rep(1 / j, n * j)), #c(rep(int.one.q1,k),rep(int.two.q2, n)),
      fun = GCE:::ObjFunGCE.primal.solnp,
      eqfun = GCE:::ConstFunGCE.primal.solnp,
      eqB = c(y[, 1], rep(1, n + k)),
      LB = rep(1e-5, k * m + n * j),
      UB = rep(1, k * m + n * j),
      control = list(
        tol = method.tol,
        inner.iter = method.maxiter,
        trace = 0
      ),
      X = X,
      n = n,
      k = k,
      m = m,
      j = j,
      q1 = int.one.q1,
      q2 = int.two.q2,
      s1 = s1,
      S = S
    )

    p <- res.opt$pars[1:(k * m)]


    ### CONFIRMAR ESTE BLOCO ########
    aux2.p <- rep(0, k * m)
    for (i in 1:k) {
      aux2.p[((i - 1) * m + 1):(i * m)] <-
        p[seq(i, i + (m - 1) * k, by = k)]
    }
    p <- aux2.p
    #################################


    b <- Z %*% p

  } else if (method == "primal.solnl") {
    res.opt <- pracma::fmincon(
      x0 = t(c(rep(1 / m, k * m), rep(1 / j, n * j))),
      fn = ObjFunGCE.primal.solnl,
      Aeq = rbind(
        cbind(as.matrix(X) %*% Z, V),
        cbind(kronecker(diag(k), matrix(1, 1, m)), matrix(0, k, n *
                                                            j)),
        cbind(matrix(0, n, k * m), kronecker(diag(n), matrix(1, 1, j)))
      ),
      beq = c(y[, 1], rep(1, n + k)),
      lb = rep(1e-5, k * m + n * j),
      ub = rep(1, k * m + n * j),
      dp = k * m,
      k = k,
      q1 = int.one.q1,
      n = n,
      q2 = int.two.q2,
      tol = method.tol,
      maxfeval = method.maxfeval,
      maxiter = method.maxiter
    )

    p <- res.opt$par[1:(k * m)]
    b <- Z %*% p

  } else if (method == "dual") {
    dimZ <- ncol(Z) #i k*m
    t <- 1
    u <- 1 #m or m9
    lambda <- increm <- matrix(0, n, 1)
    iter <- 0
    while (u > method.tol & method.maxiter > iter) {
      iter <- iter + 1
      lambda <- lambda + increm
      newz <- exp(-t(X) %*% lambda %*% matrix(1, 1, dimZ) * Z)
      p_m2 <- newz / (newz %*% matrix(1, dimZ, dimZ))
      newv <- exp(-lambda %*% matrix(1, 1, j) * S)
      w9 <- newv / (newv %*% matrix(1, j, j))
      g9 <-
        y[, 1] - as.matrix(X) %*% ((Z * p_m2) %*% matrix(1, dimZ, 1)) -
        ((S * w9) %*% matrix(1, j, 1))
      inv_z <-
        diag((apply((p_m2 * (
          Z ^ 2
        )), 1, sum) -
          apply((p_m2 * Z), 1, sum) ^ 2) ^ (-1))
      inv_v <-
        diag((apply((w9 * (
          S ^ 2
        )), 1, sum) -
          apply((w9 * S), 1, sum) ^ 2) ^ (-1))
      temp <- inv_v %*% as.matrix(X)
      inv_H = -inv_v +
        ((temp %*% solve(inv_z + t(as.matrix(
          X
        )) %*% temp)) %*% t(temp))
      increm <- inv_H %*% g9
      t0 <- t
      t <- t(g9) %*% increm
      u <- abs(t - t0)
    }

    b <- apply(p_m2 * Z, 1, sum)

    p <- matrix(0, 1, k * m)

    for (i in 1:k) {
      pos <- (i - 1) * m + 1
      p[1,pos:(pos + m - 1)] <-
        p_m2[i, pos:(pos + m - 1)] + (sum(p_m2[i,-c(pos:(pos + m - 1))])/m)
    }
    p <- t(p)[,1] #k * t(p)
  }

  q1rep <- matrix(rep(int.one.q1, k))

  nepk <- matrix(0, k, 1)
  for (i in 1:k) {
    pos <- (i - 1) * m + 1
    nepk[i, 1] = -t(p[pos:(pos + m - 1)]) %*% log(p[pos:(pos + m - 1)]) /
      (-t(q1rep[pos:(pos + m - 1)]) %*% log(q1rep[pos:(pos + m - 1)]))
  }
  nep <- (-t(p) %*% log(p)) / (-t(q1rep) %*% (log(q1rep)))

  #if (length(int.one) == 1) {
    beta_hat <- b
    if (any(c("(Intercept)", "X.Intercept.") %in% colnames(X_original))) {
      b <- as.matrix(scale_back_coef(X, y, c(0, b), intercept = TRUE))
      nepk <- c(1, nepk) #### Review ####
    } else {
      b <-
        as.matrix(scale_back_coef(X, y, b, intercept = FALSE))
    }
  #}

  if (is.null(y.test) | is.null(X.test)) {
    y.fitted <- as.matrix(X_original) %*% b
    y.values <- as.matrix(y_original)
  } else {
    y.fitted <- as.matrix(X.test) %*% b
    y.values <- as.matrix(y.test)
  }

  all.error.measures[1] <- mean((y.values - y.fitted) ^ 2)
  all.error.measures[2] <- sqrt(all.error.measures[1])
  all.error.measures[3] <- mean(abs(y.values - y.fitted))

  res <- list(
    coefficients = b,
    standardized.coefficients = beta_hat,
    residuals = {
      if (isTRUE(residuals)) {
        y.values - y.fitted
      }
    },
    fitted.values = {
      if (isTRUE(fittedvalues)) {
        y.fitted
      }
    },
    nep = nep,
    nepk = nepk,
    all.error.measures = all.error.measures,
    error.measure = all.error.measures[
      which(colnames(all.error.measures) %in% int.one.errormeasure)],
    support = {if (length(int.one) == 1) {int.one} else {
        "given"
      }},
    support.matrix = #{
      #if (length(int.one) == 1) {
        as.matrix(
          cbind(
            c(NA,scale_back_coef(X, y, intg[,1], intercept = FALSE)),
            c(NA,scale_back_coef(X, y, intg[,2], intercept = FALSE))))
      #} else {intg}}
      ,
    support.min = 1,
    p = p,
    # p.pre.dual = {
    #   if (method == "dual") {
    #     p_m2
    #   } else {
    #     NULL
    #   }
    # },
    convergence = ifelse(
      method != "dual",
      res.opt$convergence,
      ifelse(iter > method.maxiter , 1, 0)
    )
  )

  names(res$error.measure) <- int.one.errormeasure

  class(res) <- "gce_family"

  return(res)

}
