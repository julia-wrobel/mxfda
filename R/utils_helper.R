#' Get Data
#'
#' Returns data from slots within the mxfDA object
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param what character string of length 2 for what to return. produced by unlisting "uni k" or "bi g fpca" type strings
#' @param type the type of data to get for the what summary function
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @return slot data
#' @keywords internal
#'
#' @export
get_data = function(mxFDAobject, what, type){

  if(grepl("mfpca", type, ignore.case = TRUE)){
    if(grepl("[B|b]", what[1]) & grepl("[K|k]", what[2])) dat = mxFDAobject@`functional_mpca`$Kcross
    if(grepl("[B|b]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`functional_mpca`$Gcross
    if(grepl("[B|b]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobject@`functional_mpca`$Lcross
    if(grepl("[U|u]", what[1]) & grepl("[K|k]", what[2])) dat = mxFDAobject@`functional_mpca`$Kest
    if(grepl("[U|u]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`functional_mpca`$Gest
    if(grepl("[U|u]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobject@`functional_mpca`$Lest
  } else if(grepl("fpca", type, ignore.case = TRUE)){
    if(grepl("[B|b]", what[1]) & grepl("[K|k]", what[2])) dat = mxFDAobject@`functional_pca`$Kcross
    if(grepl("[B|b]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`functional_pca`$Gcross
    if(grepl("[B|b]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobject@`functional_pca`$Lcross
    if(grepl("[U|u]", what[1]) & grepl("[K|k]", what[2])) dat = mxFDAobject@`functional_pca`$Kest
    if(grepl("[U|u]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`functional_pca`$Gest
    if(grepl("[U|u]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobject@`functional_pca`$Lest
  } else if(grepl("summ", type, ignore.case = TRUE)){
    if(grepl("[B|b]", what[1]) & grepl("[K|k]", what[2])) dat = mxFDAobject@`bivariate_summaries`$Kcross
    if(grepl("[B|b]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`bivariate_summaries`$Gcross
    if(grepl("[B|b]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobject@`bivariate_summaries`$Lcross
    if(grepl("[U|u]", what[1]) & grepl("[K|k]", what[2])) dat = mxFDAobject@`univariate_summaries`$Kest
    if(grepl("[U|u]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`univariate_summaries`$Gest
    if(grepl("[U|u]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobject@`univariate_summaries`$Lest
  }
  return(dat)
}

#' Is Empty
#'
#' checks whether or not a slot is empty in the mxfDA objecct
#'
#' @param obj object with slots
#' @param slot.name character string of for the name of a slot
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @return `TRUE`, `FALSE`, or `stop` command
#' @keywords internal
#'
#' @export
is.empty <- function(obj, slot.name) {
  if (slot.name %in% methods::slotNames(class(obj))) {
    return(is.null(methods::slot(obj, slot.name)) || length(methods::slot(obj, slot.name)) == 0)
  } else {
    stop(paste0("Object does not have a '", slot.name, "' slot."))
  }
}

#' Filter Data
#'
#' not used in an important want currently
#'
#' @param dat a data frame which to filter
#' @param filter_columns named vector where names is the column and element is what to filter to
#'
#' @details `r lifecycle::badge('experimental')`
#'
#' @return data frame filtered
#' @keywords internal
#'
#' @export
#'
filter_data = function(dat, filter_columns){
  if(is.null(filter_columns)) return(dat)
  for(col in seq(filter_columns)){
    dat = dat %>%
      dplyr::filter(get(names(filter_columns[col])) == !!filter_columns[col])
  }
  return(dat)
}

#' Check if Spatial Summary Exists
#'
#' Makes sure in processing functions that the spatial summary function has been calculated before running FPCA or models
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param metric a character vector split with and passed by other functions
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @return a stopped function if slot not in summary functions
#' @keywords internal
#'
#' @export
#'
metric.exists = function(mxFDAobject, metric){
  if(grepl("[U|u]", metric[1])){
    if(!grepl(metric[2], names(mxFDAobject@univariate_summaries), ignore.case = TRUE))
      stop("Missing summary function provided")
  }
  if(grepl("[B|b]", metric[1])){
    if(!grepl(metric[2], names(mxFDAobject@bivariate_summaries), ignore.case = TRUE))
      stop("Missing summary function provided")
  }
}

#' 1/0 Checker
#'
#' @param vec vector of hopefully 1s and 0s but anything
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @return `TRUE` or `FALSE` depending on if something other than 0/1 is in the vector
#' @keywords internal
#'
#' @export
#'
one_zero = function(vec){
  vec = unique(vec)
  length(setdiff(vec, c(0, 1))) == 0
}
