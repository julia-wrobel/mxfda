get_data = function(mxFDAobject, what, type){

  if(grepl("pca", type, ignore.case = TRUE)){
    if(grepl("[B|b]", what[1]) & grepl("[K|k]", what[2])) dat = mxFDAobject@`functional_pca`$Kcross
    if(grepl("[B|b]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`functional_pca`$Gcross
    if(grepl("[B|b]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobject@`functional_pca`$Lcross
    if(grepl("[U|u]", what[1]) & grepl("[K|k]", what[2])) dat = mxFDAobject@`functional_pca`$Kest
    if(grepl("[U|u]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`functional_pca`$Gest
    if(grepl("[U|u]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobject@`functional_pca`$Lest
  }
  if(grepl("summ", type, ignore.case = TRUE)){
    if(grepl("[B|b]", what[1]) & grepl("[K|k]", what[2])) dat = mxFDAobject@`bivariate_summaries`$Kcross
    if(grepl("[B|b]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`bivariate_summaries`$Gcross
    if(grepl("[B|b]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobject@`bivariate_summaries`$Lcross
    if(grepl("[U|u]", what[1]) & grepl("[K|k]", what[2])) dat = mxFDAobject@`univariate_summaries`$Kest
    if(grepl("[U|u]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`univariate_summaries`$Gest
    if(grepl("[U|u]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobjectx@`univariate_summaries`$Lest
  }
  return(dat)
}

is.empty <- function(obj, slot.name) {
  if (slot.name %in% methods::slotNames(class(obj))) {
    return(is.null(methods::slot(obj, slot.name)) || length(methods::slot(obj, slot.name)) == 0)
  } else {
    stop(paste0("Object does not have a '", slot.name, "' slot."))
  }
}
