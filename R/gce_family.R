# Generic accessor functions of the gce class

#' Print gce
#'
#' Print gce
#'
#' @param object fitted \code{gce} object.
#' @param digits  significant digits in printout.
#' @param ... additional print arguments.
#'
#'
#' @return A small summary of a \code{gce} object is returned.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     verbose = FALSE)
#'
#' res_gce_package
#' print(res_gce_package)
#'
#' @export

print.gce <-
  function(object, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    if (length(coef(object))) {
      cat("Coefficients:\n")
      print.default(format(coef(object), digits = digits), print.gap = 2L,
                    quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(object)
  }

#' Extract Model Formula from gce object
#'
#' Returns the model used to fit gce object.
#'
#' @param object fitted \code{gce} object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     verbose = FALSE)
#'
#' formula(res_gce_package)
#' # y ~ .
#'
#' @export
formula.gce <- function(object) {
  object$formula
}

#' Extract gce Model Coefficients
#'
#' Extract coefficients from a gce object
#'
#' @param object Fitted \code{gce} model object.
#'
#' @return Returns the coefficients from a gce object
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     verbose = FALSE)
#'
#' coef(res_gce_package)
#' @export

coef.gce <- function(object) {
  aux.coef <- object$model$coefficients$res.model.coefficients
  names(aux.coef) <- rownames(object$model$coefficients)
  aux.coef
}

#' @rdname coefficients.gce
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     verbose = FALSE)
#'
#' coefficients(res_gce_package)
#' @export
coefficients.gce <- coef.gce

#' Extract gce Model Residuals
#'
#' \code{residuals} is a function which extracts model residuals from gce
#' objects.
#' The abbreviated form \code{resid} is an alias for \code{residuals}.
#'
#' @param object Fitted \code{gce} model object.
#'
#' @return Returns the residuals from a gce object
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     verbose = FALSE)
#'
#' residuals(res_gce_package)
#'
#' @export

residuals.gce <- function(object,
                          type = c("working", "response", "deviance",
                                          "pearson")) {
  type <- match.arg(type)
  r <- as.numeric(object$model$residuals)
  if (any(is.na(object$model$model))) {
    names(r) <- 1:length(r) } else {
      names(r) <-  rownames(object$model$model)}
  res <- switch(type, working = , response = r, deviance = ,
                pearson = if (is.null(object$weights)) r else r * sqrt(object$weights))

  res
}

#' @rdname resid.gce
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     verbose = FALSE)
#'
#' resid(res_gce_package)
#'
#' @export
resid.gce <- residuals.gce

#' Model Deviance
#'
#' Returns the deviance of a fitted model \code{gce} object.
#'
#' @param object Fitted \code{gce} model object.
#'
#' @return The value of the deviance extracted from a \code{gce} object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     verbose = FALSE)
#'
#' deviance(res_gce_package)
#'
#' @export
deviance.gce <- function (object)
  sum(residuals.gce(object)^2, na.rm = TRUE)

#' Variable Names of gce Fitted Models
#'
#' Simple utility returning variable names.
#'
#' @param object Fitted \code{gce} model object.
#'
#' @return A character vector.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     verbose = FALSE)
#'
#' variable.names(res_gce_package)
#'
#' @export
variable.names.gce <- function(object) {
  rownames(object$model$coefficients)
}

#' Case Names of gce Fitted Models
#'
#' Simple utility returning case names.
#'
#' @param object Fitted \code{gce} model object.
#'
#' @return A character vector.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     verbose = FALSE)
#'
#' case.names(res_gce_package)
#'
#' @export

case.names.gce <- function(object) {
  rownames(object$model$model)
}

#' Calculate gce Fitted Values
#'
#' The fitted values for the linear model represented by a gce object are
#' extracted.
#'
#' @param object Fitted \code{gce} model object.
#'
#' @return Returns a vector with the fitted values for the linear model
#' represented by a gce object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     verbose = FALSE)
#'
#' fitted(res_gce_package)
#'
#' @export

fitted.gce <- function(object) {
  aux.fitted <- as.numeric(object$model$fitted.values)
  if (any(is.na(object$model$model))) {
    names(aux.fitted) <- 1:length(aux.fitted) } else {
      names(aux.fitted) <-  rownames(object$model$model)
    }
  aux.fitted
}

#' @rdname fitted.values.gce
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     verbose = FALSE)
#'
#' fitted.values(res_gce_package)
#'
#' @export
fitted.values.gce <- fitted.gce

#' Residual Degrees-of-Freedom
#'
#' Returns the residual degrees-of-freedom extracted from a fitted model gce
#' object.
#'
#' @param object Fitted \code{gce} model object.
#'
#' @return The value of the residual degrees-of-freedom extracted from a gce
#' object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     verbose = FALSE)
#'
#' df.residual(res_gce_package)
#'
#' @export

df.residual.gce <- function(object) {
  df.residual(object$model)
}

#' Confidence Intervals for gce Model Parameters
#'
#' Computes bootstrap confidence intervals for one or more parameters in a gce
#' fitted model.
#'
#' @param object Fitted \code{gce} model object.
#' @param parm a specification of which parameters are to be given confidence
#' intervals, either a vector of numbers or a vector of names. If missing,
#' all parameters are considered.
#' @param level the confidence level required.
#' @param method .
#'
#' @return A matrix (or vector) with columns giving lower and upper confidence
#' limits for each parameter. These will be labelled as (1-level)/2 and
#' 1 - (1-level)/2 in % (by default 2.5% and 97.5%).
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     ci = 95,
#'     verbose = FALSE)
#'
#' confint(res_gce_package)
#' confint(res_gce_package, parm = c("X001"), level = 0.99)
#'
#' @export

confint.gce <- function(object,
                        parm,
                        level = 0.95,
                        method = "basic") {
  if (is.null(object$model$conf.int)) {
    return(NULL)
  } else {
    pnames <- row.names(object$model$coefficients)
    if (missing(parm)) {
      parm <- pnames
    } else
      if (is.numeric(parm)) {
        parm <- pnames[parm]
      }
    if (!all(parm %in% pnames)) {
      stop("Invalid choice of parameter")
    }
    if (!is.numeric(level)) {
      stop("Non numeric level!")
    } else {
      if (level > 1 | level < 0) {
        stop("level must be greater than 0 smaller than 1")
      }
    }
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
    if (method == "basic") {
      return(t(apply(
        object$model$conf.int[parm, ], 1, quantile, a
      )))
    }
  }
}

#' Summarise a linear regression model via generalized cross entropy fit
#'
#' Generic function used to produce summary information from a fitted linear
#' regression model via generalized cross entropy as represented by
#' \code{object} of class \code{gce}.
#'
#' @param object Fitted \code{gce} model object.
#' @param ci.level the confidence level required to compute a bootstrap
#' confidence interval.
#' @param ci.method the method used to compute a bootstrap a confidence interval.
#'
#' @return The function computes and returns a list of summary statistics of the
#'  fitted linear model given in \code{object}, using the components
#'  (list elements) "call" and "terms" from its argument, plus...
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     ci = 95,
#'     verbose = FALSE)
#'
#' sm_res_gce_package <- summary(res_gce_package)
#' str(sm_res_gce_package)
#' sm_res_gce_package_99 <- summary(res_gce_package, ci.level = 0.99,
#'  ci.method = "basic")
#' sm_res_gce_package_99$coefficients
#' @export
summary.gce <- function(object,
                        ci.level = 0.95,
                        ci.method = "basic") {
  if (is.null(object$terms))
    stop("invalid 'gce' object:  no 'terms' component")
  if (!inherits(object, "gce"))
    warning("calling summary.gce(<fake-gce-object>) ...")

  rdf <- object$model$df.residual
  r <- as.numeric(object$model$residuals)
  n <- length(r)
  p <- nrow(object$model$coefficients)
  f <- object$model$fitted.values

  mss <-
    if ("(Intercept)" %in% rownames(object$model$coefficients))
      sum((f - mean(f)) ^ 2)
  else
    sum(f ^ 2)
  rss <- sum(r ^ 2)

  resvar <- rss / rdf

  p1 <- 1L:p
  ans <- object[c("call", "terms")]
  ans$residuals <- r

  if ("conf.int" %in% names(object$model)) {
    LB <- confint.gce(object, level = ci.level, method = ci.method)[, 1]
    UB <- confint.gce(object, level = ci.level, method = ci.method)[, 2]
  } else {
    LB <- NA
    UB <- NA
  }

  ans$coefficients <- cbind(
    Estimate = object$model$coefficients$res.model.coefficients,
    `LB` = LB,
    `UB` = UB
  )

  colnames(ans$coefficients) <- c(
    "Estimate",
    paste("LB ", ci.level * 100, "% CI", sep =
            ""),
    paste("UB ", ci.level * 100, "% CI", sep =
            "")
  )

  row.names(ans$coefficients) <-
    row.names(object$model$coefficients)
  ans$aliased <-
    is.na(object$model$coefficients$res.model.coefficients)
  names(ans$aliased) <- row.names(object$model$coefficients)
  ans$sigma <- sqrt(resvar)
  ans$df <- c(p, rdf)
  if (p != attr(object$terms, "intercept")) {
    df.int <- if (attr(object$terms, "intercept"))
      1L
    else
      0L
    ans$r.squared <- mss / (mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n -
                                                       df.int) / rdf)
  }
  else
    ans$r.squared <- ans$adj.r.squared <- 0

  ans$support <- object$best$support
  ans$cv.error.measure <- object$best$measure
  names(ans$cv.error.measure) <-
    names(object$results$model$error.measure)
  class(ans) <- "summary.gce"
  ans
}

#' Print Summary of gce Model Fits
#'
#' \code{print.summary} method for class "gce".
#'
#' @param object Fitted \code{gce} model object.
#' @param digits The number of significant digits to use when printing..
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The function returns a list of summary statistics of the
#'  fitted linear model given in object, using the components (list elements)
#'  "call" and "terms" from its argument, plus...
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     ci = 95,
#'     verbose = FALSE)
#'
#' summary(res_gce_package)
#' summary(res_gce_package, ci.level = 0.99, ci.method = "basic")
#'
#' @export
print.summary.gce <-
  function(object, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\nCall:\n",
        paste(deparse(object$call), sep = "\n", collapse = "\n"),
        "\n\n",
        sep = "")
    resid <- object$residuals
    df <- object$df
    rdf <- df[2L]
    cat("Residuals:\n", sep = "")

    if (rdf > 5L) {
      nam <- c("Min", "1Q", "Median", "3Q", "Max")
      rq <- if (length(dim(resid)) == 2L)
        structure(
          apply(t(resid), 1L, quantile), dimnames = list(nam,
                                                         dimnames(resid)[[2L]]))
      else {
        zz <- zapsmall(quantile(resid), digits + 1L)
        structure(zz, names = nam)
      }
      print(rq, digits = digits, ...)
    }
    else if (rdf > 0L) {
      print(resid, digits = digits, ...)
    }
    else {
      cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
      cat("\n")
    }

    if (length(object$aliased) == 0L) {
      cat("\nNo Coefficients\n")
    }
    else {
      cat("\nCoefficients:\n")
      coefs <- object$coefficients
      if (any(aliased <- object$aliased)) {
        cn <- names(aliased)
        coefs <- matrix(
          NA,
          length(aliased),
          1,
          dimnames = list(cn,colnames(coefs)))
        coefs[!aliased,] <- object$coefficients
      }
      print(coefs)
    }
    cat("\nResidual standard error:",
        format(signif(object$sigma,
                      digits)),
        "on",
        rdf,
        "degrees of freedom")
    cat("\n")
    cat("Multiple R-squared: ",
        formatC(object$r.squared, digits = digits))
    cat(",\tAdjusted R-squared: ",
        formatC(object$adj.r.squared,
                digits = digits))
    cat("\n")
    cat("Choosen support: ", object$support)
    cat("\n")
    cat(
      "CV-",
      names(object$cv.error.measure),
      ":  ",
      formatC(object$cv.error.measure,
              digits = digits),
      sep = ""
    )
    cat("\n")
  }

#' Plot Diagnostics for a gce Object
#'
#' Three plots are currently available: a plot of observed against fitted
#' values, a Scale-Location plot of residuals against fitted values and a Q-Q
#' plot of residuals.
#'
#' @param object Fitted \code{gce} model object.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' res_gce_package <-
#'   gce(y~.,
#'     data = gce_k3_yk2_cn1_df,
#'     model = TRUE,
#'     method = "dual",
#'     ci = 95,
#'     verbose = FALSE)
#'
#' plot(res_gce_package)
#'
#' @export
plot.gce <-
  function (object,
            which = c(1, 2, 3, 4),
            prior = TRUE,
            caption = list("Residuals vs Fitted",
                           "Q-Q Residuals",
                           "Scale-Location",
                           "Error vs supports"),
            panel = {if (add.smooth) {function(x, y, ...){panel.smooth(x,
                                                                     y,
                                                                     iter = iter.smooth,
                                                                     ...)}} else points},
            sub.caption = NULL,
            main = "",
            ask = prod(par("mfcol")) < length(which) && dev.interactive(),
            id.n = 3,
            labels.id = names(residuals(object)),
            cex.id = 0.75,
            qqline = TRUE,
            add.smooth = getOption("add.smooth"),
            iter.smooth = 3,
            label.pos = c(4, 2),
            cex.caption = 1,
            cex.oma.main = 1.25,
            extend.ylim.f = 0.08,
            ...) {

    dropInf <- function(x, h) {
      if (any(isInf <- h >= 1)) {
        warning(gettextf("not plotting observations with leverage one:\n  %s",
                         paste(which(isInf), collapse = ", ")), call. = FALSE,
                domain = NA)
        x[isInf] <- NaN
      }
      x
    }

  if (!inherits(object, "gce"))
    stop("use only with \"gce\" objects")
  if (!is.numeric(which) || any(which < 1) || any(which > 4))
    stop("'which' must be in 1:6")

  show <- rep(FALSE, 4)
  show[which] <- TRUE

  r <- residuals(object)
  yh <- fitted(object)
  n <- length(r)

  #rs <- r

  hii <- hat(object$model$model[,-1])
  rs <- dropInf(r,hii)

  rds <- rs
  l.fit <- "Fitted values"
  if (is.null(id.n))
    id.n <- 0L
  else {
    id.n <- as.integer(id.n)
    if (id.n < 0L || id.n > n)
      stop(gettextf("'id.n' must be in {1,..,%d}", n),
           domain = NA)
  }
  if (id.n > 0L) {
    if (is.null(labels.id))
      labels.id <- paste(1L:n)
    iid <- 1L:id.n
    show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
    if (any(show[2L:3L])) {
      show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid]
      show.rds <- sort.list(abs(rds), decreasing = TRUE)[iid]
    }
    text.id <- function(x, y, ind, adj.x = TRUE, usr = par("usr")) {
      labpos <- if (adj.x)
        label.pos[(x > mean(usr[1:2])) + 1L]
      else 3
      text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE,
           pos = labpos, offset = 0.25)
    }
  }
  getCaption <- function(k) if (length(caption) < k)
    NA_character_
  else as.graphicsAnnot(caption[[k]])
  if (is.null(sub.caption)) {
    cal <- object$call
    if (!is.na(m.f <- match("formula", names(cal)))) {
      cal <- cal[c(1, m.f)]
      names(cal)[2L] <- ""
    }
    cc <- deparse(cal, 80)
    nc <- nchar(cc[1L], "c")
    abbr <- length(cc) > 1 || nc > 75
    sub.caption <- if (abbr)
      paste(substr(cc[1L], 1L, min(75L, nc)), "...")
    else cc[1L]
  }
  one.fig <- prod(par("mfcol")) == 1
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }


  if (show[1L]) {
    ylim <- range(r, na.rm = TRUE)
    if (id.n > 0)
      ylim <- extendrange(r = ylim, f = extend.ylim.f)
    dev.hold()
    ylab1 <- "Residuals"
    plot(yh, r, xlab = l.fit, ylab = ylab1, main = main,
         ylim = ylim, type = "n", ...)
    panel(yh, r, ...)
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(getCaption(1), 3, 0.25, cex = cex.caption)
    if (id.n > 0) {
      y.id <- r[show.r]
      y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
      text.id(yh[show.r], y.id, show.r)
    }
    abline(h = 0, lty = 3, col = "gray")
    dev.flush()
  }

  if (show[2L]) {

    ylim <- range(rds, na.rm = TRUE)
    ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
    dev.hold()
      qq <- qqnorm(rds, main = main, ylab = "Standardized residuals", ylim = ylim,
                   ...)
      if (qqline)
        qqline(rds, lty = 3, col = "gray50")

    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(getCaption(2), 3, 0.25, cex = cex.caption)
    if (id.n > 0)
      text.id(qq$x[show.rds], qq$y[show.rds], show.rds)
    dev.flush()


    # qqnorm(r,
    #        main = "Q-Q Residuals")
    # qqline(r, lty = 3, col = "gray50")
  }

  if (show[3L]) {
    sqrtabsr <- sqrt(abs(rs))
    ylim <- c(0, max(sqrtabsr, na.rm = TRUE))
    yl <- as.expression(substitute(sqrt(abs(YL)),
                                   list(YL = as.name("Standardized residuals"))))
    yhn0 <- yh
    dev.hold()
    plot(yhn0, sqrtabsr, xlab = l.fit, ylab = yl, main = main,
         ylim = ylim, type = "n", ...)
    panel(yhn0, sqrtabsr, ...)
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(getCaption(3), 3, 0.25, cex = cex.caption)
    if (id.n > 0)
      text.id(yhn0[show.rs], sqrtabsr[show.rs], show.rs)
    dev.flush()

    # plot(yh,
    #      sqrtabsr,
    #      xlab = "",
    #      ylab = "",
    #      main = "Scale-location",
    #      #ylim = ylim,
    #      type = "n")
  }

  if (show[4L] & !(object$int.one.method %in% c("unique"))) {
    x <- round(as.numeric(object$ok.support), 3)
    y <- object$mean.error.measure.cv[1:length(x)]
    sdev <- object$sd.error.measure.cv[1:length(x)]
    nep <- NULL
    nnep <- length(object$ok.support)

    for (i in 1:nnep) {
      nep <- c(nep, object$results$results[[i]]$nep)
    }

    if (!prior & object$int.one.method %in% c("adaptive", "specific")) {
      x <- x[-1]
      y <- y[-1]
      sdev <- sdev[-1]
      nep <- nep[-1]
    }

    # par(mar = c(5,5,2,5))
    plot(
      x,
      y,
      ylim = range(c(y - sdev, y + sdev)),
      xlab = "Upper limit of the support spaces",
      ylab = names(object$model$error.measure),
      type = "n",
      ...
    )

    arrows(
      x,
      y - sdev,
      x,
      y + sdev,
      length = 0.05,
      angle = 90,
      code = 3,
      ...
    )

    abline(v = round(as.numeric(
      gsub("support.", "", object$best$support)
    ), 3),
    lty = 3,
    ...)

    points(
      x,
      y,
      pch = 20,
      #cex = 1.5,
      col = {if (!prior & object$int.one.method %in% c("adaptive", "specific")) {
        rep("red",length(x) - 1)
      } else {c("blue",rep("red",length(x) - 1))}},
      ...)

    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(getCaption(4), 3, 0.25, cex = cex.caption)

    # par(new = T)
    # plot(
    #   x,
    #   nep,
    #   axes = FALSE,
    #   xlab = NA,
    #   ylab = NA,
    #   pch=16,
    #   col = "blue"
    # )
    # axis(side = 4)
    # mtext(side = 4, line = 3, 'normalized entropy')
  }
  if (!one.fig && par("oma")[3L] >= 1)
    mtext(sub.caption, outer = TRUE, cex = cex.oma.main)
  invisible()
}
