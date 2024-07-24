#' Run function Cox models for data with multiple samples per subject
#'
#' Fit a functional Cox regression model when there are multiple functions per subject, which arise from multiple samples per subject.
#' It is not necessary for all subjects to have the same number of samples.The function first performs a multilevel functional principal components
#' analysis (MFPCA) decomposition to the spatial summary function. Then, the average curve for each subject is used in a functional Cox model (FCM).
#' Variation around each subject's mean is captured by calculating the standard deviation of the level 2 scores from MFPCA, then including this as
#' a scalar variable in the FCM called "level2_score_sd".
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
#' @param filter_cols a named vector of factors to filter summary functions to in `c(Derived_Column = "Level_to_Filter")` format
#' @param pve Proportion of variance explained by multilevel functional principal components analysis in mfpca step
#' @param ... Optional other arguments to be passed to \code{fpca.face}
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @return A \code{list} which is a linear or additive functional Cox model fit. See \code{mgcv::gam} for more details.
#'
#' @author Julia Wrobel \email{`r juliawrobel_email`}
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @examples
#' #load ovarian mxFDA object
#' data('lung_FDA')
#'
#' # run the lfcm model
#' lung_FDA = run_mfcm(lung_FDA, model_name = "fit_mlfcm",
#'                       formula = survival_days ~ age,
#'                       event = "survival_status",
#'                       metric = "uni g", r = "r", value = "fundiff",
#'                       pve = 0.99,
#'                       afcm = FALSE)
#' @export
run_mfcm <- function(mxFDAobject,
                    model_name,
                    formula,
                    event = "event",
                    metric = "uni k",
                    r = "r",
                    value = "fundiff",
                    afcm = FALSE,
                    filter_cols = NULL,
                    pve = 0.99,
                    ...,
                    knots = NULL){

  if(!(event %in% colnames(mxFDAobject@Metadata)))
    stop("`event` needs to be a column name in the metadata")
  if(!one_zero(mxFDAobject@Metadata[[event]]))
    stop("The event column needs to be in format 0/1")
  # check for missing values in the functional predictor

#### FOR ALEX TO ADD!
  ## first check if mfda has been run. Right now it only works off of my specific example (K function)
  mfpca_obj <- run_mfpca(mxFDAobject,
                         metric = metric,
                         r = r,
                         value = value,
                         pve = pve,
                         lightweight = TRUE)


  analysis_vars <- unique(c(all.vars(formula), event, mxFDAobject@subject_key))
  meta_vars <- mxFDAobject@Metadata %>% select(all_of(analysis_vars)) %>% distinct()
  metric = unlist(strsplit(metric, split = " "))
  #get summary func name
  if(grepl("[B|b]", metric[1]) & grepl("[K|k]", metric[2])) func = "Kcross"
  if(grepl("[B|b]", metric[1]) & grepl("[G|g]", metric[2])) func = "Gcross"
  if(grepl("[B|b]", metric[1]) & grepl("[L|l]", metric[2])) func = "Lcross"
  if(grepl("[U|u]", metric[1]) & grepl("[K|k]", metric[2])) func = "Kest"
  if(grepl("[U|u]", metric[1]) & grepl("[G|g]", metric[2])) func = "Gest"
  if(grepl("[U|u]", metric[1]) & grepl("[L|l]", metric[2])) func = "Lest"

  ### this is specific to the K example I have been working with
  mxfundata = mfpca_obj@functional_mpca[[func]]$Yi_hat %>%
    dplyr::left_join(meta_vars)

  analysis_vars <- c(analysis_vars, "level2_score_sd")


  form = deparse(stats::formula(formula))
  mxfundata <- process_fcm(mxfundata, mxFDAobject@subject_key, r, value, analysis_vars, multilevel = TRUE)


  # fit linear or additive functional Cox model
  if(afcm){
    form =  paste0(form, '+ level2_score_sd + ti(t_int, func, by=l_int, bs=c("cr","cr"), k=c(10,10), mc=c(FALSE,TRUE))')
    weights = mxfundata[[event]]
    fit_fcm <- mgcv::gam(formula = stats::as.formula(form),
                         weights = weights,
                         data = mxfundata,
                         family = mgcv::cox.ph())

    class(fit_fcm) <- append("afcm", class(fit_fcm))
  }else{
    form =  paste0(form, '+ level2_score_sd + s(t_int, by=l_int*func, bs="cr", k=20)')
    weights = mxfundata[[event]]
    fit_fcm <- mgcv::gam(formula = stats::as.formula(form),
                   weights = weights,
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

