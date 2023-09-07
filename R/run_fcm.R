#' process_fcm
#'
#' Fit a functional Cox regression model.
#'
#' @importFrom tidyr pivot_wider
#' @importFrom mgcv gam cox.ph
#' @import dplyr
#'
#' @return A \code{list} which is a linear or additive functional Cox model fit. See \code{mgcv::gam} for more details.
#'
#' @examples
#' # simulate data
#' set.seed(1001)
#'
#'
#' @param mxfundata Dataframe of spatial summary functions from multiplex imaging data, in long format. Can be estimated using the function \code{extract_summary_functions} or provided separately.
#' @param form Formula to be fed to mgcv in the form of survival_time ~ x1 + x2. Does not contain functional predictor. Character valued. Data must contain censoring variable called "event".
#' @param id Character string, the name of the variable that identifies each unique subject.
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param knots Number of knots for defining spline basis.
#' @param afcm If TRUE, runs additive functional Cox model. If FALSE, runs linear functional cox model. Defaults to linear functional cox model.
#' @param analysis_vars List of variables to be retained for downstream analysis, including variables from formula.
#' @param quantile_transform Defaults to FALSE. If TRUE, a quantile transformation is applied to the functional predictor before modeling using the \code{process_fcm} function.
#' @param smooth Option to smooth data using FPCA. Defaults to FALSE.
#' @param ... Optional other arguments to be passed to \code{fpca.face}
#' @export
run_fcm <- function(mxfundata,
                        form,
                        id,
                        r = "r",
                        value = "fundiff",
                    afcm = FALSE,
                        analysis_vars,
                        quantile_transform = FALSE,
                    smooth = FALSE,
                        ...){

  # check for missing values in the functional predictor
  if(smooth){
    mxfundata <- impute_fpca(mxfundata, id = id, r = r, value = value,
                             analysis_vars = analysis_vars, smooth = TRUE)
    message("Functional predictor contains NA values that were imputed using FPCA")
  }
  if(anyNA(mxfundata[[value]])){
    mxfundata <- impute_fpca(mxfundata, id = id, r = r, value = value,
                               analysis_vars = analysis_vars)
    message("Functional predictor contains NA values that were imputed using FPCA")
  }

  mxfundata <- process_fcm(mxfundata, id, r, value, analysis_vars, quantile_transform)

  # fit linear or additive functional Cox model
  if(afcm){
    form =  paste0(form, '+ ti(t_int, func, by=l_int, bs=c("cr","cr"), k=c(10,10), mc=c(FALSE,TRUE))')

    fit_fcm <- gam(formula = as.formula(form),
                   weights = event,
                   data = mxfundata, family = cox.ph())

    class(fit_fcm) <- append("afcm", class(fit_fcm))
  }else{
    form =  paste0(form, '+ s(t_int, by=l_int*func, bs="cr", k=20)')

    fit_fcm <- gam(formula = as.formula(form),
                   weights = event,
                   data = mxfundata, family = cox.ph())

    class(fit_fcm) <- append("lfcm", class(fit_fcm))
  }

 fit_fcm
}

