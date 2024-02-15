#' Extract Model
#'
#' Currently only extracts functional cox models not mixed functional cox models.
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param metric metric functional PCA data to extract, e.g. 'uni k'
#' @param type one of "cox", "mcox", or "sofr" to specify the type of model to extract
#' @param model_name character string of the model name to retrieve
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @return fit functional model
#'
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @examples
#' #load ovarian mxFDA object
#' data('ovarian_FDA')
#'
#' #run the lfcm model
#' ovarian_FDA = run_fcm(ovarian_FDA, model_name = "fit_lfcm",
#'                       formula = survival_time ~ age, event = "event",
#'                       metric = "uni g", r = "r", value = "fundiff",
#'                       analysis_vars = c("age", "survival_time"),
#'                       afcm = FALSE)
#'
#' #extract model
#' mod = extract_model(ovarian_FDA, 'uni g', 'cox', 'fit_lfcm')
#'
#' @export
extract_model = function(mxFDAobject, metric, type, model_name){
  #data preparation
  if(is.null(model_name))
    stop("Must provide a model")
  #get the right metric
  if(length(metric) != 1)
    stop("Please provide a single spatial metric to extract surface for")
  #get the right model type
  if(length(type) != 1 | !(type %in% c('cox', 'mcox', 'sofr')))
    stop("Must provide appropriate model type")
  metric = unlist(strsplit(metric, split = " "))

  if(type == 'cox'){
    #get linear regression
    if(grepl("[B|b]", metric[1]) & grepl("[K|k]", metric[2])) fit = mxFDAobject@functional_cox$Kcross[[model_name]]
    if(grepl("[B|b]", metric[1]) & grepl("[G|g]", metric[2])) fit = mxFDAobject@functional_cox$Gcross[[model_name]]
    if(grepl("[B|b]", metric[1]) & grepl("[L|l]", metric[2])) fit = mxFDAobject@functional_cox$Lcross[[model_name]]
    if(grepl("[U|u]", metric[1]) & grepl("[K|k]", metric[2])) fit = mxFDAobject@functional_cox$Kest[[model_name]]
    if(grepl("[U|u]", metric[1]) & grepl("[G|g]", metric[2])) fit = mxFDAobject@functional_cox$Gest[[model_name]]
    if(grepl("[U|u]", metric[1]) & grepl("[L|l]", metric[2])) fit = mxFDAobject@functional_cox$Lest[[model_name]]
  } else if(type == 'mcox'){
    #get mixed effects
    if(grepl("[B|b]", metric[1]) & grepl("[K|k]", metric[2])) fit = mxFDAobject@functional_mcox$Kcross[[model_name]]
    if(grepl("[B|b]", metric[1]) & grepl("[G|g]", metric[2])) fit = mxFDAobject@functional_mcox$Gcross[[model_name]]
    if(grepl("[B|b]", metric[1]) & grepl("[L|l]", metric[2])) fit = mxFDAobject@functional_mcox$Lcross[[model_name]]
    if(grepl("[U|u]", metric[1]) & grepl("[K|k]", metric[2])) fit = mxFDAobject@functional_mcox$Kest[[model_name]]
    if(grepl("[U|u]", metric[1]) & grepl("[G|g]", metric[2])) fit = mxFDAobject@functional_mcox$Gest[[model_name]]
    if(grepl("[U|u]", metric[1]) & grepl("[L|l]", metric[2])) fit = mxFDAobject@functional_mcox$Lest[[model_name]]
  } else if(type == 'sofr'){
    #get scalar on functional regression
    if(grepl("[B|b]", metric[1]) & grepl("[K|k]", metric[2])) fit = mxFDAobject@scalar_on_functional$Kcross[[model_name]]
    if(grepl("[B|b]", metric[1]) & grepl("[G|g]", metric[2])) fit = mxFDAobject@scalar_on_functional$Gcross[[model_name]]
    if(grepl("[B|b]", metric[1]) & grepl("[L|l]", metric[2])) fit = mxFDAobject@scalar_on_functional$Lcross[[model_name]]
    if(grepl("[U|u]", metric[1]) & grepl("[K|k]", metric[2])) fit = mxFDAobject@scalar_on_functional$Kest[[model_name]]
    if(grepl("[U|u]", metric[1]) & grepl("[G|g]", metric[2])) fit = mxFDAobject@scalar_on_functional$Gest[[model_name]]
    if(grepl("[U|u]", metric[1]) & grepl("[L|l]", metric[2])) fit = mxFDAobject@scalar_on_functional$Lest[[model_name]]
  }

  if(is.null(fit))
    stop("Model ", model_name, " doesn't exist")

  return(fit)
}
