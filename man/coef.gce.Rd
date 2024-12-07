% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gce_family.R
\name{coef.gce}
\alias{coef.gce}
\title{Extract gce Model Coefficients}
\usage{
\method{coef}{gce}(object)
}
\arguments{
\item{object}{Fitted \code{gce} model object.}
}
\value{
Returns the coefficients from a gce object
}
\description{
Extract coefficients from a gce object
}
\examples{
res_gce_package <-
  gce(y~.,
    data = gce_k3_yk2_cn1_df,
    model = TRUE,
    method = "dual",
    verbose = FALSE)

coef(res_gce_package)
}
\author{
Jorge Cabral, \email{jorgecabral@ua.pt}
}
