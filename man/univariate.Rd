% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_univariate.R
\name{univariate}
\alias{univariate}
\title{univariate}
\usage{
univariate(
  mximg,
  markvar,
  mark1,
  mark2,
  r_vec,
  func = c(Kest, Lest, Gest),
  edge_correction,
  breaks = NULL
)
}
\arguments{
\item{mximg}{Dataframe of cell-level multiplex imaging data for a single image.
Should have variables \code{x} and \code{y} to denote x and y spatial locations of each cell.}

\item{markvar}{The name of the variable that denotes cell type(s) of interest. Character.}

\item{mark1}{dummy filler, unused}

\item{mark2}{dummy filler, unused}

\item{r_vec}{Numeric vector of radii over which to evaluate spatial summary functions. Must begin at 0.}

\item{func}{Spatial summary function to calculate. Options are c(Kest, Lest, Gest) which denote Ripley's K, Besag's L, and nearest neighbor G function, respectively.}

\item{edge_correction}{Character string that denotes the edge correction method for spatial summary function. For Kest and Lest choose one of c("border", "isotropic", "Ripley", "translate", "none"). For Gest choose one of c("rs", "km", "han")}

\item{breaks}{ignored}
}
\value{
A \code{data.frame} containing:
\item{r}{the radius of values over which the spatial summary function is evaluated}
\item{sumfun}{the values of the spatial summary function}
\item{csr}{the values of the spatial summary function under complete spatial randomness}
\item{fundiff}{sumfun - csr, positive values indicate clustering and negative values repulsion}
}
\description{
Internal function called by \code{\link[=extract_summary_functions]{extract_summary_functions()}} to calculate a univariate spatial summary function for a single image.
}
\details{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\author{
Julia Wrobel \email{julia.wrobel@emory.edu}
}
