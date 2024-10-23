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
#' @param cv.nfolds number of folds used for cross-validation. The default is
#' \code{cv.nfolds = 5} and the smallest value allowable is \code{cv.nfolds = 3}.
#' @param cv.repeats number of repetitions for cross-validation. The default is
#' \code{cv.repeats = 1}.
#' @param cv.results Boolean value. If \code{TRUE} cross-validation models are
#' returned. The default is \code{cv.results = TRUE}.
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
#' @param seed A single value, interpreted as an integer, for reproducibility
#' or \code{NULL}. The default is \code{seed = NULL}.
#' @param verbose An integer to control how verbose the output is. For a value
#' of 0 no messages or output are shown and for a value of 3 all messages
#' are shown. The default is \code{verbose = 0}.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd
gce_estimation_cv <- function(y,
                              X,
                              y.test = NULL,
                              X.test = NULL,
                              residuals = TRUE,
                              fittedvalues = TRUE,
                              cv.nfolds = 5,
                              cv.repeats = 1,
                              cv.results = TRUE,
                              int.one = NULL,
                              int.one.q1 = c(1 / 3, 1 / 3, 1 / 3),
                              int.one.errormeasure = "RMSE",
                              int.two = NULL,
                              int.two.q2 = c(1 / 3, 1 / 3, 1 / 3),
                              method = "dual",
                              method.maxfeval = 1e+04,
                              method.maxiter = 1e+05,
                              method.tol = 1e-06,
                              seed = NULL,
                              verbose = 0) {
  if (verbose >= 2) {
    cat(0, "% ", sep = "")
  }

  mean.measure <- NULL
  res <- list(
    cvresults = list(),
    model = NULL,
    mean.error.measure.cv = NULL,
    sd.error.measure.cv = NULL,
    support = NULL,
    support.min = NULL
  )

  res$model <-
    gce_estimation(
      y,
      X,
      y.test,
      X.test,
      residuals,
      fittedvalues,
      int.one,
      int.one.q1,
      int.one.errormeasure,
      int.two,
      int.two.q2,
      method,
      method.maxfeval,
      method.maxiter,
      method.tol
    )

  if (verbose >= 2) {
    cat(round(1 / (cv.repeats * cv.nfolds + 1) * 100, 0), "% ", sep = "")
  }

  for (r in 1:cv.repeats) {
    res$cvresults[[r]] <- list()

    if (!is.null(seed)) {
      set.seed(seed + r)
    }

    #change_order <- sample(nrow(X))
    yX_CV <-
      groupdata2::fold(data.frame(y = as.numeric(y[, 1]), X),
                                      k = cv.nfolds)
      # data.frame(y = as.numeric(y[change_order, 1]),
      #            X[change_order, ],
      #            .folds = cut(seq(1, nrow(X)),
      #                         breaks = cv.nfolds,
      #                         labels = FALSE))

    for (cv in 1:cv.nfolds) {
      res$cvresults[[r]][[cv]] <- list()

      res_k <-
        gce_estimation(
          y = data.frame(V1 = yX_CV[yX_CV$.folds != cv, 1]),
          X = as.matrix(yX_CV[yX_CV$.folds != cv, -c(1, ncol(yX_CV))]),
          y.test = data.frame(V1 = yX_CV[yX_CV$.folds == cv, 1]),
          X.test = as.matrix(yX_CV[yX_CV$.folds == cv, -c(1, ncol(yX_CV))]),
          residuals,
          fittedvalues,
          int.one,
          int.one.q1,
          int.one.errormeasure,
          int.two,
          int.two.q2,
          method,
          method.maxfeval,
          method.maxiter,
          method.tol
        )

      mean.measure <- c(mean.measure, res_k$error.measure[[1]])
      res$cvresults[[r]][[cv]] <- res_k

      if (verbose >= 2) {
        cat(round(((r - 1) * cv.nfolds + cv + 1) / (cv.repeats * cv.nfolds + 1) *
                    100, 0), "% ", sep = "")
      }
    }
    names(res$cvresults[[r]]) <- paste0("fold", 1:cv.nfolds)
  }

  names(res$cvresults) <- paste0("repeats", 1:cv.repeats)

  res$mean.error.measure.cv <- mean(mean.measure)
  res$sd.error.measure.cv <- sd(mean.measure)

  res$support = int.one
  res$support.min = 1

  if (!isTRUE(cv.results)) {
    res$cvresults <- NA
  }

  class(res) <- "gce_family"

  return(res)

}
