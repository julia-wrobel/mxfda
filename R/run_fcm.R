#' process_fcm
#'
#' Fit a functional Cox regression model.
#'
#'
#' @param mxFDAobject Dataframe of spatial summary functions from multiplex imaging data, in long format. Can be estimated using the function \code{extract_summary_functions} or provided separately.
#' @param model_name character string to give the fit model in the functional cox slot
#' @param formula Formula to be fed to mgcv in the form of survival_time ~ x1 + x2. Does not contain functional predictor. Character valued. Data must contain censoring variable called "event".
#' @param event character string for the column in Metadata that contains 1/0 for the survival event
#' @param metric name of calculated spatial metric to use
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param knots Number of knots for defining spline basis.
#' @param afcm If TRUE, runs additive functional Cox model. If FALSE, runs linear functional cox model. Defaults to linear functional cox model.
#' @param analysis_vars List of variables to be retained for downstream analysis, including variables from formula.
#' @param quantile_transform Defaults to FALSE. If TRUE, a quantile transformation is applied to the functional predictor before modeling using the \code{process_fcm} function.
#' @param smooth Option to smooth data using FPCA. Defaults to FALSE.
#' @param filter_cols a named vector of factors to filter summary functions to in `c(Derived_Column = "Level_to_Filter")` format
#' @param ... Optional other arguments to be passed to \code{fpca.face}
#'
#' @return A \code{list} which is a linear or additive functional Cox model fit. See \code{mgcv::gam} for more details.
#'
#' @author Julia Wrobel \email{julia.wrobel@@emory.edu}
#' @author Alex Soupir \email{alex.soupir@@moffitt.org}
#'
#' @examples
#' #load ovarian mxFDA object
#' data('ovarian_FDA')
#'
#' #run ghe lfcm model
#' ovarian_FDA = run_fcm(ovarian_FDA, model_name = "fit_lfcm",
#'                       formula = survival_time ~ age, event = "event",
#'                       metric = "uni g", r = "r", value = "fundiff",
#'                       analysis_vars = c("age", "survival_time"),
#'                       afcm = FALSE)
#'
#' @export
run_fcm <- function(mxFDAobject,
                    model_name,
                    formula,
                    event = "event",
                    metric = "uni k",
                    r = "r",
                    value = "fundiff",
                    afcm = FALSE,
                    analysis_vars = NULL,
                    quantile_transform = FALSE,
                    smooth = FALSE,
                    filter_cols = NULL,
                    ...,
                    knots = NULL){
  #get the right data from the object
  if(length(metric) != 1)
    stop("Please provide a single spatial metric to calculate functional cox models with")
  metric = unlist(strsplit(metric, split = " "))
  #check for slot in summaries
  metric.exists(mxFDAobject, metric)

  mxfundata = get_data(mxFDAobject, metric, 'summaries') %>%
    filter_data(filter_cols) %>%
    dplyr::full_join(mxFDAobject@Metadata, by = mxFDAobject@sample_key)
  #join everything needed to fit the model into a vector for analysis vars
  analysis_vars = unique(c(analysis_vars,
                           grep("~", paste0(formula), invert = TRUE, value = TRUE),
                           event))

  if(!(event %in% colnames(mxFDAobject@Metadata)))
    stop("`event` needs to be a column name in the metadata")
  if(!one_zero(mxFDAobject@Metadata[[event]]))
    stop("The event column needs to be in format 0/1")
  # check for missing values in the functional predictor
  if(smooth){
    mxfundata <- impute_fpca(mxfundata, id = mxFDAobject@sample_key, r = r, value = value,
                             analysis_vars = analysis_vars, smooth = TRUE)
    message("Functional predictor contains NA values that were imputed using FPCA")
  }
  if(anyNA(mxfundata[[value]])){
    mxfundata <- impute_fpca(mxfundata, id = mxFDAobject@sample_key , r = r, value = value, knots = knots,
                               analysis_vars = analysis_vars, smooth)
    message("Functional predictor contains NA values that were imputed using FPCA")
  }
  form = deparse(stats::formula(formula))
  mxfundata <- process_fcm(mxfundata, mxFDAobject@sample_key, r, value, analysis_vars, quantile_transform)

  # fit linear or additive functional Cox model
  if(afcm){
    form =  paste0(form, '+ ti(t_int, func, by=l_int, bs=c("cr","cr"), k=c(10,10), mc=c(FALSE,TRUE))')

    fit_fcm <- mgcv::gam(formula = stats::as.formula(form),
                         weights = mxfundata[[event]],
                         data = mxfundata,
                         family = mgcv::cox.ph())

    class(fit_fcm) <- append("afcm", class(fit_fcm))
  }else{
    form =  paste0(form, '+ s(t_int, by=l_int*func, bs="cr", k=20)')

    fit_fcm <- mgcv::gam(formula = as.formula(form),
                   weights = mxfundata[[event]],
                   data = mxfundata,
                   family = mgcv::cox.ph())

    class(fit_fcm) <- append("lfcm", class(fit_fcm))
  }

  if(grepl("[B|b]", metric[1]) & grepl("[K|k]", metric[2])) mxFDAobject@`functional_cox`$Kcross[[model_name]] = fit_fcm
  if(grepl("[B|b]", metric[1]) & grepl("[G|g]", metric[2])) mxFDAobject@`functional_cox`$Gcross[[model_name]] = fit_fcm
  if(grepl("[B|b]", metric[1]) & grepl("[L|l]", metric[2])) mxFDAobject@`functional_cox`$Lcross[[model_name]] = fit_fcm
  if(grepl("[U|u]", metric[1]) & grepl("[K|k]", metric[2])) mxFDAobject@`functional_cox`$Kest[[model_name]] = fit_fcm
  if(grepl("[U|u]", metric[1]) & grepl("[G|g]", metric[2])) mxFDAobject@`functional_cox`$Gest[[model_name]] = fit_fcm
  if(grepl("[U|u]", metric[1]) & grepl("[L|l]", metric[2])) mxFDAobject@`functional_cox`$Lest[[model_name]] = fit_fcm

 return(mxFDAobject)
}

