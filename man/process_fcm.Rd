% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funreg_utils.R
\name{process_fcm}
\alias{process_fcm}
\title{process_fcm}
\usage{
process_fcm(
  mxfundata,
  id,
  r = "r",
  value = "fundiff",
  analysis_vars,
  quantile_transform = FALSE,
  multilevel = FALSE
)
}
\arguments{
\item{mxfundata}{Dataframe of spatial summary functions from multiplex imaging data, in long format. Can be estimated using the function \code{extract_summary_functions} or provided separately.}

\item{id}{Character string, the name of the variable that identifies each unique subject.}

\item{r}{Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".}

\item{value}{Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".}

\item{analysis_vars}{Optional list of variables to be retained for downstream analysis.}

\item{quantile_transform}{If TRUE, a quantile transformation is applied to the functional predictor before modeling}
}
\value{
A \code{dataframe} with matrix-valued covariates \code{l_int}, \code{t_int}, and \code{func} for use in a linear or additive functional Cox model.
}
\description{
Internal function called by \code{run_fcm} that transforms long format functional data for use in a linear or additive functional Cox model.
}
\examples{
# simulate data
set.seed(1001)

}
\author{
Julia Wrobel \email{julia.wrobel@emory.edu}

Alex Soupir \email{alex.soupir@moffitt.org}
}
\keyword{internal}
