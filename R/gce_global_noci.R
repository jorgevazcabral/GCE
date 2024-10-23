#' Generalized Cross entropy estimation
#'
#' Internal function to fit a linear regression model via generalized cross
#' entropy where initial support spaces can be provided or computed.
#'
#' @param formula An object of class \code{\link[stats]{formula}} (or one that
#'  can be coerced to that class): a symbolic description of the model to be
#'  fitted.
#' @param data A data frame (or object coercible by
#' \code{\link[stats]{as.data.frame}} to a data frame) containing the variables
#'  in the model.
#' @param data.test A data frame (or object coercible by
#' \code{\link[stats]{as.data.frame}} to a data frame) containing the same
#' variables as \code{data} to be used to determine the out-of sample error. If
#' \code{NULL} the error will be computed using \code{data}.
#' @param model Boolean value. if \code{TRUE}, the model frame used is returned.
#' The default is \code{model = FALSE}.
#' @param residuals Boolean value. If \code{TRUE} residuals are returned.  The
#' default is \code{residuals = TRUE}.
#' @param fittedvalues Boolean value. If \code{TRUE} fitted values are returned.
#' The default is \code{fittedvalues = TRUE}.
#' @param cv Boolean value. If \code{TRUE} cross-validation will be
#' performed. The default is \code{cv = TRUE}.
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
#' @param int.one.prior One of c("msc", "lm", "ridge"). The default is
#' \code{int.one.prior = "msc"}.
#' @param int.one.percentile A percentile value to be used in the computation of the
#'  maximum value of the absolute standardized coefficient.
#' @param int.one.factor A positive value to be used in the adjustment of the
#' computed maximum value of the absolute standardized coefficient.
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
#' @param seed A single value, interpreted as an integer, for reproducibility
#' or \code{NULL}. The default is \code{seed = NULL}.
#' @param verbose An integer to control how verbose the output is. For a value
#' of 0 no messages or output are shown and for a value of 3 all messages
#' are shown. The default is \code{verbose = 0}.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd
gce_global_noci <- function(formula,
                            data,
                            data.test = NULL,
                            model = FALSE,
                            residuals = TRUE,
                            fittedvalues = TRUE,
                            cv = TRUE,
                            cv.nfolds = 5,
                            cv.repeats = 1,
                            cv.results = TRUE,
                            int.one = NULL,
                            int.one.q1 = c(1 / 3, 1 / 3, 1 / 3),
                            int.one.prior = "msc",
                            int.one.percentile = 99.9,
                            int.one.factor = 5,
                            int.one.method = "sweep",
                            int.one.alpha = NULL,
                            int.one.alpha.n = 100,
                            int.one.alpha.results = TRUE,
                            int.one.errormeasure = "RMSE",
                            int.two = NULL,
                            int.two.q2 = c(1 / 3, 1 / 3, 1 / 3),
                            method = "dual",
                            method.maxfeval = 1e+04,
                            method.maxiter = 1e+05,
                            method.tol = 1e-06,
                            seed = NULL,
                            verbose = 0) {

  data <- as.data.frame(data)

  X <- as.data.frame(model.matrix(formula, data = data))
  y <- as.data.frame(data[row.names(data) %in% row.names(X),
                          colnames(data) == all.vars(formula)[1]])
  colnames(y) <- as.character(all.vars(formula)[1])

  if (!is.null(data.test)) {
    data.test <- as.data.frame(data.test)
    X.test <- as.data.frame(model.matrix(formula, data = data.test))
    y.test <-
      data.test[row.names(data.test) %in% row.names(X.test),
                colnames(data.test) == all.vars(formula)[1]]
    colnames(y) <- as.character(formula[2])
  } else {
    X.test <- NULL
    y.test <- NULL
  }

  # missing <- apply(apply(X, 1, is.na), 2, any)
  # if (any(missing)) {
  #   y <- data.frame(y[-which(missing),])
  #   X <- na.omit(X)
  # }

  if (any(c("(Intercept)", "X.Intercept.") %in% colnames(X))) {
    X.aux <- scale(X[, -1])
  } else {
    X.aux <- scale(X)
  }

  k <- ncol(X.aux)

  if (length(int.one) == 0) {
    if (k > 1) {
      if (int.one.prior == "ridge") {
        int.one <-
          glmnet::cv.glmnet(
            x = X.aux,
            y = scale(y) ,
            standardize = FALSE,
            intercept = FALSE,
            alpha = 0,
            lambda = 10^seq(log(10^(-10),base = 10),
                            log(10^(2),base = 10),
                            length.out = 100),
            nfolds = 5,
            type.measure = "default"
          )

        int.one <- round(int.one.factor * max(c(1,
                                          max(abs((int.one$glmnet.fit$beta))))),
                         3)
        } else
      if (int.one.prior == "lm") {

        int.one <-
          lm(data = data.frame(y = scale(y),
                               X.aux),
             formula = formula)

        int.one <-
          round(int.one.factor * max(c(1,
                                       abs(coef(
                                         int.one
                                       )))),
                3)
      } else if (int.one.prior == "msc") {

      max_cor_xx <- max(abs(cor(X.aux)[upper.tri(cor(X.aux))]))
      int.one <-
        round(
          int.one.factor * max_standardized_coef(k, max_cor_xx, int.one.percentile),3)
    }
      } else {
      int.one <- round(int.one.factor * 1,3)
    }
  }

  if (length(int.one) == 1 & length(int.one.alpha) == 0) {
    if (int.one.method == "sweep") {
      int.one.alpha <-
            round(exp(seq(log(int.one),log(1),length.out = int.one.alpha.n)),3)
    } else if (int.one.method %in% c("adaptive", "specific")) {
      prior_estimation <-
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

      if (int.one.method == "adaptive") {
      int.one.alpha <-
        c(int.one,
          round(
            max(
              abs(
                prior_estimation$standardized.coefficients)) * exp(
                  seq(
                    log(20),
                    log(2),
                    length.out = int.one.alpha.n - 1)),
            3))
      } else {
        if (method == "dual") {
          int.one.alpha <-
            cbind(rep(int.one,k),
                  round(
                    outer(abs(prior_estimation$standardized.coefficients),
                          exp(seq(log(20), log(2), length.out = int.one.alpha.n - 1)),
                          "*"),
                    3))
        } else {
      int.one.alpha <-
        cbind(rep(int.one,k),
              round(
                outer(abs(prior_estimation$standardized.coefficients)[,1],
                      exp(seq(log(20), log(2), length.out = int.one.alpha.n - 1)),
                      "*"),
                3))
        }
      colnames(int.one.alpha) <-
        c(int.one,
          round(exp(seq(log(20),log(2),length.out = int.one.alpha.n - 1)),3))
    }}

  } else if (length(int.one) == 1 & length(int.one.alpha) != 0) {
    int.one.alpha <- c(int.one,int.one.alpha)
    }

  if (isTRUE(cv)) {
    if (int.one.method == "unique") {
      res <-
        gce_estimation_cv(
          y,
          X,
          y.test,
          X.test,
          residuals,
          fittedvalues,
          cv.nfolds,
          cv.repeats,
          cv.results,
          int.one,
          int.one.q1,
          int.one.errormeasure,
          int.two,
          int.two.q2,
          method,
          method.maxfeval,
          method.maxiter,
          method.tol,
          seed,
          verbose)
    } else {
      res <-
        gce_estimation_sscv(
          y,
          X,
          y.test,
          X.test,
          residuals,
          fittedvalues,
          cv.nfolds,
          cv.repeats,
          cv.results,
          int.one,
          int.one.q1,
          int.one.method,
          int.one.alpha,
          int.one.alpha.results,
          int.one.errormeasure,
          int.two,
          int.two.q2,
          method,
          method.maxfeval,
          method.maxiter,
          method.tol,
          seed,
          verbose)
    }
  }
  else {
    if (int.one.method == "unique") {
      res <- list(model = list())
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
          method.tol)
    } else {
      res <-
        gce_estimation_ss(
          y,
          X,
          y.test,
          X.test,
          residuals,
          fittedvalues,
          int.one,
          int.one.q1,
          int.one.method,
          int.one.alpha,
          int.one.alpha.results,
          int.one.errormeasure,
          int.two,
          int.two.q2,
          method,
          method.maxfeval,
          method.maxiter,
          method.tol,
          verbose
        )
    }
  }

  res$model$coefficients <- data.frame(res$model$coefficients)

  rownames(res$model$coefficients) <- colnames(X)

  res$model$df.residual <- nrow(X) - ncol(X)
  if (isTRUE(model)) {
    res$model$model <- data.frame(y, X)
    res$model$model <-
      res$model$model[, colnames(res$model$model) != "X.Intercept."]
  } else {
    res$model$model <- NA
  }

  class(res) <- "gce_family"

  return(res)
}
