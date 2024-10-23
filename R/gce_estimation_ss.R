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
#' @param int.one.method One of c("unique", "adaptive", "sweep"). The default is
#' \code{int.one.method = "sweep"}.
#' @param int.one.alpha A vector of positive upper limits for the support spaces
#' when \code{int.one.method = "sweep"} or a vector of values to multiply
#' estimated coefficients when \code{int.one.method = "adaptive"}. The default
#' is \code{int.one.alpha = NULL}.
#' @param int.one.alpha.n A positive value for the length of
#' \code{int.one.alpha} when \code{int.one.alpha = NULL}. The default is
#'  \code{int.one.alpha.n = 100}.
#' @param int.one.alpha.results Boolean value. If \code{TRUE} models for different
#' supports are returned. The default is \code{int.one.alpha.results = TRUE}.
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
#' @param verbose An integer to control how verbose the output is. For a value
#' of 0 no messages or output are shown and for a value of 3 all messages
#' are shown. The default is \code{verbose = 0}.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd
gce_estimation_ss <- function(y,
                              X,
                              y.test = NULL,
                              X.test = NULL,
                              residuals = TRUE,
                              fittedvalues = TRUE,
                              int.one = NULL,
                              int.one.q1 = c(1 / 3, 1 / 3, 1 / 3),
                              int.one.method = "sweep",
                              int.one.alpha = NULL,
                              int.one.alpha.results = TRUE,
                              int.one.errormeasure = "RMSE",
                              int.two = NULL,
                              int.two.q2 = c(1 / 3, 1 / 3, 1 / 3),
                              method = "dual",
                              method.maxfeval = 1e+04,
                              method.maxiter = 1e+05,
                              method.tol = 1e-06,
                              verbose = 0) {
  if (verbose >= 3) {
    cat("0% ", sep = "")
  }

  if (int.one.method == "specific") {
    names.supports <- colnames(int.one.alpha)
  } else {
    names.supports <- int.one.alpha
  }
  nsupports <- length(names.supports)

  res <- list(
    results = list(list()),
    error.measure = NULL,
    best = NULL,
    support = names.supports,
    support.min = NULL,
    model = NULL
  )

  for (su in 1:nsupports) {

    aux.int.one.alpha <- int.one.alpha[su]

    if (int.one.method == "specific") {
      aux.int.one.alpha <- int.one.alpha[,su]
    } else {
      aux.int.one.alpha <- int.one.alpha[su]
    }

    res$results[[su]] <-
      tryCatch({
        gce_estimation(
          y,
          X,
          y.test,
          X.test,
          residuals,
          fittedvalues,
          int.one = aux.int.one.alpha,
          int.one.q1,
          int.one.errormeasure,
          int.two,
          int.two.q2,
          method,
          method.maxfeval,
          method.maxiter,
          method.tol)
      #res$results[[su]]$support <- as.numeric(names.supports[su])
      },
      error = function(e) {
        return(list(error.measure = NA))
      })

    if ((length(res$results[[su]]) == 1) &
        (su != nsupports)) {
      for (h in (su + 1):(nsupports)) {
        res$results[[h]] <- list(error.measure = NA)
      }
      if (verbose >= 3) {
        cat("100% ", sep = "")
      }
      break
    }

    if (verbose >= 3) {
      cat(round((su) / nsupports * 100, 0), "% ", sep = "")
    }
  }

  names(res$results) <- names.supports

  measuresupports <- NULL
  for (ms in 1:nsupports) {
    measuresupports <- c(measuresupports,
                         res$results[[ms]]$error.measure)
  }

  res$error.measure <- measuresupports
  names(res$error.measure) <- names.supports

  res$best <-
    data.frame(support = names.supports[which.min(measuresupports)],
               measure = round(min(measuresupports,na.rm = TRUE), 4))
  res$support.min <- which.min(measuresupports)[[1]]
  res$model <- res$results[[which.min(measuresupports)]]

  if (!isTRUE(int.one.alpha.results)) {
    res$results <- NA
  }

  class(res) <- "gce_family"

  return(res)

}
