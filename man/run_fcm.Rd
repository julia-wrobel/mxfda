% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_fcm.R
\name{run_fcm}
\alias{run_fcm}
\title{process_fcm}
\usage{
run_fcm(
  mxfundata,
  form,
  id,
  r = "r",
  value = "fundiff",
  afcm = FALSE,
  analysis_vars,
  quantile_transform = FALSE,
  ...
)
}
\arguments{
\item{mxfundata}{Dataframe of spatial summary functions from multiplex imaging data, in long format. Can be estimated using the function \code{extract_summary_functions} or provided separately.}

\item{form}{Formula to be fed to mgcv in the form of survival_time ~ x1 + x2. Does not contain functional predictor. Character valued. Data must contain censoring variable called "event".}

\item{id}{Character string, the name of the variable that identifies each unique subject.}

\item{r}{Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".}

\item{value}{Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".}

\item{afcm}{If TRUE, runs additive functional Cox model. If FALSE, runs linear functional cox model. Defaults to linear functional cox model.}

\item{analysis_vars}{List of variables to be retained for downstream analysis, including variables from formula.}

\item{quantile_transform}{Defaults to FALSE. If TRUE, a quantile transformation is applied to the functional predictor before modeling using the \code{process_fcm} function.}

\item{...}{Optional other arguments to be passed to \code{fpca.face}}

\item{knots}{Number of knots for defining spline basis.}
}
\value{
A \code{list} which is a linear or additive functional Cox model fit. See \code{mgcv::gam} for more details.
}
\description{
Fit a functional Cox regression model.
}
\examples{
# simulate data
set.seed(1001)


}
