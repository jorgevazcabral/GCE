#' Scale coefficients back
#'
#' Returns the unscaled regression coefficients
#'
#' @param X.scaled A matrix.
#' @param y.scaled A vector.
#' @param betas.scaled A vector of given scaled coefficients.
#' @param intercept logical indicating if
#'
#' @return Returns a vector of unscaled numeric regression coefficients.
#'
#' @author Jorge Cabral, \email{jorgecabral@@ua.pt}
#'
#' @examples
#'
#' @export

scale_back_coef <-
  function (X.scaled,
            y.scaled,
            betas.scaled,
            intercept = TRUE) {
    k <- attributes(X.scaled)[[1]][2]

    sigma_X <- attributes(X.scaled)$`scaled:scale`
    sigma_y <- attributes(y.scaled)$`scaled:scale`

    aux <-
      sapply(1:k, function(feat) {
        attributes(X.scaled)$`scaled:center`[feat] * sigma_y * betas.scaled[feat + 1] /
          sigma_X[feat]
      })

    if (isTRUE(intercept)) {
      betas <- rep(0, k + 1)
      betas[1] <-
        betas.scaled[1] * sigma_y + attributes(y.scaled)$`scaled:center` - sum(aux)
      betas[-1] <-
        sapply(1:k, function(feat) {
          betas.scaled[feat + 1] * sigma_y / sigma_X[feat]
        })
    } else {
      betas <- rep(0, k)
      betas <-
        sapply(1:k, function(feat) {
          betas.scaled[feat] * sigma_y / sigma_X[feat]
        })
    }

    return(betas)
  }

