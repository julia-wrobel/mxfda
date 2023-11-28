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

is.empty <- function(obj, slot.name) {
  if (slot.name %in% methods::slotNames(class(obj))) {
    return(is.null(methods::slot(obj, slot.name)) || length(methods::slot(obj, slot.name)) == 0)
  } else {
    stop(paste0("Object does not have a '", slot.name, "' slot."))
  }
}

filter_data = function(dat, filter_columns){
  if(is.null(filter_columns)) return(dat)
  for(col in seq(filter_columns)){
    dat = dat %>%
      dplyr::filter(get(names(filter_columns[col])) == !!filter_columns[col])
  }
  return(dat)
}

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

one_zero = function(vec){
  vec = unique(vec)
  length(setdiff(vec, c(0, 1))) == 0
}
