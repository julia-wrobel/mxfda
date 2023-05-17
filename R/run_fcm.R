#' process_fcm
#'
#' Internal function called by \code{TITLE: regression function} that transforms long format functional data for use in a linear or additive functional Cox model.
#'
#' @importFrom tidyr pivot_wider
#' @importFrom mgcv gam
#' @import dplyr
#'
#' @return A \code{dataframe} with matrix-valued covariates \code{l_int}, \code{t_int}, and \code{func} for use in a linear or additive functional Cox model.
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
#' @param quantile_transform Defaults to FALSE. If TRUE, a quantile transformation is applied to the functional predictor before modeling using the \code{process_fcm} function.
#' @param ... Optional other arguments to be passed to \code{fpca.face}
#' @export
run_fcm <- function(mxfundata,
                        form,
                        survival_status,
                        id,
                        r = "r",
                        value = "fundiff",
                        analysis_vars,
                        quantile_transform = FALSE,
                        ...){

  # check for missing values in the functional predictor
  if(anyNA(mxfundata[[value]])){
    mxfundata <- impute_fpca(mxfundata, id = id, r = r, value = value,
                               analysis_vars = analysis_vars)
    message("Functional predictor contains NA values that were imputed using FPCA")
  }

  mxfundata <- process_fcm(mxfundata, id, r, value, analysis_vars, quantile_transform)

  # transform to matrix
  # fit specified model

  # TO DO: try providing formula in the form of a Surv object, so that you don't
  # have so many arguments about specifying variables. BUT for now do what is fastest.

  fit_afc <- gam(survival_days ~ stage + adjuvant_therapy +
                   ti(t_int, G, by=l_int, bs=c("cr","cr"), k=c(10,10), mc=c(FALSE,TRUE)),
                 weights=survival_status,
                 data=Gdf, family=cox.ph())


  if(quantile_transform){
    # quantile transformed version of the data
    mxfundata$func <- apply(func, 2, function(y) ecdf(y)(y))
  }else{
    mxfundata$func <- func
  }


}

