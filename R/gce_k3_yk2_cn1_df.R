#' Simulated data set generated with gdata_GCE
#'
#' Simulated data, used to demonstrate the features of gce
#'
#' @format A \code{data.frame} containing:
#' \describe{
#'   \item{X001}{A N(0,1) independent variable.}
#'   \item{X002}{A N(0,1) independent variable.}
#'   \item{X003}{A N(0,1) independent variable.}
#'   \item{X004}{A N(0,1) independent variable.}
#'   \item{X005}{A N(0,1) independent variable.}
#'   \item{y}{A Dependent variable}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' # gce_k3_yk2_cn1_df <-
#' #   gdata_GCE(
#' #     n = 100,
#' #     cont.k = 3,
#' #     y.gen.cont.k = 2,
#' #     y.gen.cont.beta = c(5, 10),
#' #     intercept.beta = 0,
#' #     Xgenerator = "svd",
#' #     condnumber = 1,
#' #     mu = 0,
#' #     sd = 1,
#' #     error_dist = "normal",
#' #     error_normal_mean = 0,
#' #     error.dist.snr = 5,
#' #     dataframe = TRUE,
#' #     seed = 230676
#' #   )
#' data(gce_k3_yk2_cn1_df)

"gce_k3_yk2_cn1_df"
