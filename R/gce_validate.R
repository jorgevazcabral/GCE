#' Validation of gce arguments
#'
#' This function validates the arguments for \code{\link[gce]}
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
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @noRd

gce_validate <- function(formula,
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

  if (missing(formula)) stop("argument `formula` is required.",call. = FALSE)
  if (missing(data)) stop("argument `data` is required.",call. = FALSE)

  if (is.null(int.one.prior)) {
    stop('argument `int.one.prior` must be one of c("msc", "lm", "ridge")',
         call. = FALSE)
  } else if (!int.one.prior %in% c("msc", "lm", "ridge")) {
    stop('argument `int.one.prior` must be one of c("msc", "lm", "ridge")',
         call. = FALSE)
  }

  if (is.null(int.one.method)) {
    stop('argument `int.one.method` must be one of c("unique", "adaptive", "sweep", "specific")',
         call. = FALSE)
  } else if (!int.one.method %in% c("unique", "adaptive", "sweep", "specific")) {
    stop('argument `int.one.method` must be one of c("unique", "adaptive", "sweep", "specific")',
         call. = FALSE)
  }

  if (int.one.method == "sweep" & length(int.one) > 1) {
    stop('when using `int.one.method == "sweep"` argument `int.one` must be one of NULL or a scalar.',
         call. = FALSE)
  }

  if (is.null(method)) {
    stop('argument `method` must be one of c("primal.solnl", "primal.solnp", "dual")',
         call. = FALSE)
  } else if (!method %in% c("primal.solnl", "primal.solnp", "dual")) {
    stop('argument `method` must be one of c("primal.solnl", "primal.solnp", "dual")',
         call. = FALSE)
  }

  if (is.null(int.one.errormeasure)) {
    stop('argument `int.one.errormeasure` must be one of c("RMSE","MAE","MSE")',
         call. = FALSE)
  } else if (!int.one.errormeasure %in% c("RMSE","MAE","MSE")) {
    stop('argument `int.one.errormeasure` must be one of c("RMSE","MAE","MSE")',
         call. = FALSE)
  }

  logical.param <- c(
    model = model,
    cv = cv,
    residuals = residuals,
    fittedvalues = fittedvalues,
    cv.results = cv.results,
    int.one.alpha.results = int.one.alpha.results
  )

  if (!is.logical(logical.param) |
      length(logical.param) != 6) {
    stop(
      "The following arguments must be TRUE or FALSE:
          `model`,`cv`, `residuals`, `fittedvalues`,
          `cv.results`, `int.one.alpha.results`",
      call. = FALSE
    )
  }

  if (is.null(verbose) |
      is.null(ci.verbose)) {
    stop("arguments `verbose` and `ci.verbose` must be one of 0, 1, 2 or 3.",
         call. = FALSE)
  } else if (!(verbose %in% c(0, 1, 2, 3)) |
             !(ci.verbose %in% c(0, 1, 2, 3))) {
    stop("arguments `verbose` and `ci.verbose` must be one of 0, 1, 2 or 3.",
         call. = FALSE)
  }

  if (is.null(cv.nfolds)) {
    stop("argument `cv.nfolds` must be an integer equal or greater than 3.",
         call. = FALSE)
  } else if (cv.nfolds < 3 | cv.nfolds %% 1 != 0) {
    stop("argument `cv.nfolds` must be an integer equal or greater than 3.",
         call. = FALSE)
  }

  if (is.null(cv.repeats)) {
    stop("argument `cv.repeats` must be a positive integer.", call. = FALSE)
  } else if (cv.repeats < 1 | cv.repeats %% 1 != 0) {
    stop("argument `cv.repeats` must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(int.one.q1) | !is.numeric(int.two.q2)) {
    stop("arguments `int.one.q1` and `int.two.q2` must be numeric vectors.", call. = FALSE)
  } else if (sum(int.one.q1) != 1 | sum(int.two.q2) != 1) {
    stop("the sum of the elements of arguments `int.one.q1` or `int.two.q2` must be equal to 1.",
         call. = FALSE)
  } else if (length(int.one.q1) == 1 | length(int.two.q2) == 1) {
    stop("the number of elements of arguments `int.one.q1` or `int.two.q2` must be greater than 1.",
         call. = FALSE)
  }

  if (!is.numeric(ci.B)) {
    stop("argument `ci.B` must be 0 or a positive integer greater or equal to 10.",
         call. = FALSE)
    if (ci.B < 10) {
      stop("argument `ci.B` must be 0 or a positive integer greater or equal to 10.",
           call. = FALSE)
  }}
}
