#' Estimated maximum absolute regression standardized coefficient
#'
#' Returns the estimated maximum absolute regression standardized coefficient
#'
#' @param k Number of parameters
#' @param maxcor Maximum absolute correlation
#' @param percentile Percentile
#'
#' @return Returns an \code{percentile} based estimated maximum absolute value
#' for a standardized regression coefficient given a linear regression model
#' with \code{k} parameters and a maximum correlation value of \code{maxcor}.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#' max_standardized_coef(k = 3, maxcor = .70, percentile = 90)
#' # [1] 1.024147
#'
#' @export

max_standardized_coef <- function(k, maxcor, percentile) {
  percentile <- percentile / 100

  b <-
    -4.753 + 2.451 * k + 15.938 * maxcor + 6.513 * percentile - 9.198 * k *
    maxcor - 3.424 * k * percentile - 20.898 * maxcor * percentile + 11.710 *
    k * maxcor * percentile

  if (exp(b) < 1) {return(1)} else {return(exp(b))}

}
