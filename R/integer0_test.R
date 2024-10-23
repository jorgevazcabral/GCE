#' Return 0
#'
#' Returns the value 0 when numeric(0) or identical(0) is introduced
#'
#' @param data value
#'
#' @return Returns 0 or data.
#'
#' @examples
#' x = numeric(0)
#' y = 2
#' x*y
#' integer0_test(x*y)
#' # [1] 0
#' integer0_test(y)
#' # [1] 2
#'
#' @noRd

integer0_test <- function(data) {
  if (identical(data, integer(0)) | identical(data, numeric(0))) {
    return(0)
  }

  else {
    return(data)
  }
}
