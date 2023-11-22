#' Extract Model
#'
#' Currently only extracts functional cox models not mixed functional cox models.
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param what what functional PCA data to extract, e.g. 'uni k'
#' @param model_name character string of the model name to retrieve
#'
#' @return fit functional model
#'
#' @author Alex Soupir \email{alex.soupir@@moffitt.org}
#'
#' @examples
#' #set seed
#' set.seed(333)
#'
#' @export
extract_model = function(mxFDAobject, what, model_name){
  #data preparation
  if(is.null(model_name))
    stop("Must provide a model")
  #get the right data
  if(length(what) != 1)
    stop("Please provide a single spatial metric to extract surface for")
  what = unlist(strsplit(what, split = " "))

  #get the model
  if(grepl("[B|b]", what[1]) & grepl("[K|k]", what[2])) fit = mxFDAobject@functional_cox$Kcross[[model_name]]
  if(grepl("[B|b]", what[1]) & grepl("[G|g]", what[2])) fit = mxFDAobject@functional_cox$Gcross[[model_name]]
  if(grepl("[B|b]", what[1]) & grepl("[L|l]", what[2])) fit = mxFDAobject@functional_cox$Lcross[[model_name]]
  if(grepl("[U|u]", what[1]) & grepl("[K|k]", what[2])) fit = mxFDAobject@functional_cox$Kest[[model_name]]
  if(grepl("[U|u]", what[1]) & grepl("[G|g]", what[2])) fit = mxFDAobject@functional_cox$Gest[[model_name]]
  if(grepl("[U|u]", what[1]) & grepl("[L|l]", what[2])) fit = mxFDAobject@functional_cox$Lest[[model_name]]

  return(fit)
}
