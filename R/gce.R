#' Generalized Cross entropy estimation
#'
#' This generic function fits a linear regression model via generalized cross
#' entropy. Initial support spaces can be provided or computed.
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
#' @param method Use \code{"primal.solnl"} (GCE using SQP) or
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
#' @param ci.B A single positive integer greater or equal to 10 for the number
#'  of bootstrap samples to be used for the computation of the boostrap
#'  confidence interval(s). Zero value will generate no sample and therefore no
#'  confidence interval(s). The default is \code{ci.B = 0}.
#' @param ci.verbose An integer to control how verbose the cross-validation
#' output is. For a value of 0 no messages or output are shown and for a
#' value of 3 all messages are shown. The default is \code{ci.verbose = 0}.
#' @param seed A single value, interpreted as an integer, for reproducibility
#' or \code{NULL}. The default is \code{seed = NULL}.
#' @param verbose An integer to control how verbose the output is. For a value
#' of 0 no messages or output are shown and for a value of 3 all messages
#' are shown. The default is \code{verbose = 0}.
#'
#' @details
#'
#' The \code{gce} function fits a linear regression model via generalized cross
#' entropy
#'
#' \eqn{Y =f(X) + \epsilon = \beta_0 + \beta_1 X_{1} + \dots + \beta_K X_{K} + \epsilon.}{}
#'
#' Models for \code{gce} are specified symbolically. A typical model has the form
#' response ~ terms where response is the (numeric) response vector and terms is
#'  a series of terms which specifies a linear predictor for response.
#'  \code{gce} calls the lower level functions \code{gce_global_noci},
#'  \code{gce_estimation_sscv}, \code{gce_estimation_ss},
#'  \code{gce_estimation_cv} and \code{gce_estimation}.
#'
#' @return
#' \code{gce} returns an object of class \code{"gce"}.
#'  The function \code{\link{summary.gce}} is used to obtain and print a summary
#'  of the results. The generic accessor functions \code{\link{coef.gce}},
#'  \code{\link{fitted.values.gce}}, \code{\link{residuals.gce}} and
#'  \code{\link{df.residual.gce}}, extract various useful features of the value
#'   returned by \code{object} of class \code{gce}.
#'
#'  An object of class "gce" is a list containing at least the following
#'  components:
#'
#' \item{coefficients}{a named data frame of coefficients}
#' \item{residuals}{the residuals, that is response minus fitted values.}
#' \item{fitted.values}{the fitted mean values.}
#' \item{call}{the matched call.}
#' \item{terms}{the terms object used.}
#' \item{df.residual}{the residual degrees of freedom.}
#' \item{model}{if requested (the default), the model frame used.}
#'
#' @seealso
#' \code{\link{summary.gce}} for more detailed summaries.
#' The generic functions \code{\link{plot.gce}}, \code{\link{print.gce}},
#'  \code{\link{coef.gce}} and \code{\link{confint.gce}}.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @references
#' Golan, A., Judge, G. G. and Miller, D. (1996)
#' \emph{Maximum entropy econometrics : robust estimation with limited data.},
#' Wiley.\cr
#' Jaynes, E.T. (1957)
#' \emph{Information Theory and Statistical Mechanics, Physical Review,
#'  106, 620-630},
#' \doi{/10.1103/PhysRev.106.620}.
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual")
#'
#' res_gce_package
#'
#' @export

gce <- function(formula,
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
                ci.B = 0,
                ci.verbose = 0,
                seed = NULL,
                verbose = 0) {

    gce_validate(
      formula,
      data,
      data.test,
      model,
      residuals,
      fittedvalues,
      cv,
      cv.nfolds,
      cv.repeats,
      cv.results,
      int.one,
      int.one.q1,
      int.one.prior,
      int.one.percentile,
      int.one.factor,
      int.one.method,
      int.one.alpha,
      int.one.alpha.n,
      int.one.alpha.results,
      int.one.errormeasure,
      int.two,
      int.two.q2,
      method,
      method.maxfeval,
      method.maxiter,
      method.tol,
      ci.B,
      ci.verbose,
      seed,
      verbose
    )

  if (verbose >= 1) {
    cat("Estimation\n", sep = "")
  }

  cl <- match.call()

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  res <-
    GCE:::gce_global_noci(
      formula,
      data,
      data.test,
      model,
      residuals,
      fittedvalues,
      cv,
      cv.nfolds,
      cv.repeats,
      cv.results,
      int.one,
      int.one.q1,
      int.one.prior,
      int.one.percentile,
      int.one.factor,
      int.one.method,
      int.one.alpha,
      int.one.alpha.n,
      int.one.alpha.results,
      int.one.errormeasure,
      int.two,
      int.two.q2,
      method,
      method.maxfeval,
      method.maxiter,
      method.tol,
      seed,
      verbose
    )

  res$call <- cl
  res$formula <- formula
  res$terms <- mt
  res$int.one.method <- int.one.method

  if (ci.B != 0 & ci.verbose >= 1) {
    cat("\nConfidence intervals", sep = "")
  }

  if (ci.B != 0) {
    res$model$conf.int <-
      data.frame(matrix(
        NA,
        ncol = ci.B,
        nrow = nrow(res$model$coefficients)
      ))
    row.names(res$model$conf.int) <-
      row.names(res$model$coefficients)

    #for (b in 1:ci.B) {
    b = 0
    attempts.b <- 0
    complete_b <- NULL

    while (b < ci.B) {
      attempts.b <- attempts.b + 1
      b <- b + 1
      if (is.null(complete_b)) {b <- 1} else {
        if (!((b-1) %in% complete_b)) {b <- b-1}}

      if (ci.verbose >= 1) {
      cat("\nb = ",b,"\n", sep = "")
      }

      tryCatch({
      res$model$conf.int[, b] <-
        GCE:::gce_global_noci(
          formula,
          data = {
            set.seed(ifelse(is.null(seed),0,seed) + attempts.b)
            data[sample(nrow(data), nrow(data), replace = TRUE), ]
          },
          data.test,
          model,
          residuals,
          fittedvalues,
          cv,
          cv.nfolds,
          cv.repeats,
          cv.results,
          int.one,
          int.one.q1,
          int.one.prior,
          int.one.percentile,
          int.one.factor,
          int.one.method,
          int.one.alpha,
          int.one.alpha.n,
          int.one.alpha.results,
          int.one.errormeasure,
          int.two,
          int.two.q2,
          method,
          method.maxfeval,
          method.maxiter,
          method.tol,
          seed = ifelse(is.null(seed),0,seed) + attempts.b,
          ci.verbose
        )$model$coefficients

      complete_b <- c(complete_b, b)
      },
      error = function(e) {
        b <<- b-1
      })

    }
  }

  class(res) <- "gce"

  return(res)

}
