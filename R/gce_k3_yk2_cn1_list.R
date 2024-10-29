#' Simulated data set generated with gdata_GCE
#'
#' Simulated data, used to demonstrate the features of gce
#'
#' @format A list containing:
#' \describe{
#'   \item{X}{A (100x5)-matrix of N(0,1) independent variables.}
#'   \item{y}{A Dependent variable}
#'   \item{coefficients}{The coefficients used to generate y - corrigir}
#'   \item{y.coefficients}{The coefficients used to generate y  - corrigir}
#'   \item{epsilon}{The Error added to the product of X and the coefficients.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' # gce_k3_yk2_cn1_list <-
#' #   gdata_GCE(
#' #     n = 100,
#' #     cont.k = 3,
#' #     y.gen.cont.k = 2,
#' #     y.gen.cont.beta = c(5, 10),
#' #     intercept.beta = 0,
#' #     Xgenerator.method = "svd",
#' #     condnumber = 1,
#' #     mu = 0,
#' #     sd = 1,
#' #     error.dist = "normal",
#' #     error.dist.mean = 0,
#' #     error.dist.snr = 5,
#' #     seed = 230676
#' #   )
#' data(gce_k3_yk2_cn1_list)
#' data.frame(gce_k3_yk2_cn1_list$X, y = gce_k3_yk2_cn1_list$y)

"gce_k3_yk2_cn1_list"
