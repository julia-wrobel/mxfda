% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_summary_function.R
\name{add_summary_function}
\alias{add_summary_function}
\title{Add Summary Function}
\usage{
add_summary_function(mxFDAobject, summary_function_data, metric)
}
\arguments{
\item{mxFDAobject}{object of class \code{mxFDA}}

\item{summary_function_data}{data frame with \code{summary_key} from \code{mxFDA} object as key column for summary function}

\item{metric}{character vector with either 'uni' or 'bi' and 'k', 'l', or 'g'; e.g. 'uni g'}
}
\value{
an updated \code{mxFDA} object with a derived value added. See \code{\link[=make_mxfda]{make_mxfda()}} for more details.
}
\description{
Sometimes other ways of calculating summary functions is wanted and is done in other packages,
in this instance the data can be loaded into the \code{mxFDA} object.
}
\author{
Alex Soupir \email{alex.soupir@moffitt.org}
}
