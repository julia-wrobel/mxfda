#' impute_fpca
#'
#' Internal function called by \code{TITLE: regression function} that imputes missing data in functional predictors using FPCA.
#'
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
#' @importFrom refund fpca.face
#' @importFrom tidyr pivot_wider pivot_longer
#' @import dplyr
#'
#' @return A \code{dataframe} where the missing function values (NA) for the \code{value} variable have been replaced with estimates from FPCA.
#'
#' @examples
#' # simulate data
#' set.seed(1001)
#'
#'
#' @param mxfundata Dataframe of spatial summary functions from multiplex imaging data, in long format. Can be estimated using the function \code{extract_summary_functions} or provided separately.
#' @param id Character string, the name of the variable that identifies each unique subject.
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param knots Number of knots for defining spline basis.
#' @param analysis_vars Optional list of variables to be retained for downstream analysis.
#' @param ... Optional other arguments to be passed to \code{fpca.face}
impute_fpca = function(mxfundata,
                    id,
                    r = "r",
                    value = "fundiff",
                    knots = NULL,
                    analysis_vars,
                    ...){

  mxfundata <- mxfundata %>%
    select(all_of(c(id, r, value, analysis_vars))) %>%
    pivot_wider(names_from = r,
                names_glue = "r_{round(r, 2)}",
                values_from = value)

  mat <- mxfundata %>%
    select(starts_with("r_")) %>%
    as.matrix()

  if(is.null(knots)) knots = floor(ncol(mat)/3)
  if(knots > ncol(mat) - 5) knots = floor(ncol(mat)/3)

  # run fpca
  mx_fpc <- fpca.face(Y = mat, knots = knots, ...)

  # impute missing values using FPCA fitted values
  mxfundata = mxfundata %>%
    # convert to long format
    pivot_longer(contains("r_"), names_to = r, values_to = value) %>%
    mutate(Yhat = as.vector(t(mx_fpc$Yhat)),
           imputed = ifelse(is.na(get(value)), Yhat, get(value)))

  mxfundata[[value]] <- mxfundata$imputed
  #mxfundata[[value]] <- mxfundata$value


  select(mxfundata, -imputed, -Yhat)
}



