get_data = function(mxFDAobject, what, type){

  if(grepl("func", type, ignore.case = TRUE)){
    if(grepl("[B|b]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`Functional Data`$Gcross
    if(grepl("[B|b]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobject@`Functional Data`$Lcross
    if(grepl("[U|u]", what[1]) & grepl("[K|k]", what[2])) dat = mxFDAobject@`Functional Data`$Kest
    if(grepl("[U|u]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`Functional Data`$Gest
    if(grepl("[U|u]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobject@`Functional Data`$Lest
  }
  if(grepl("summ", type, ignore.case = TRUE)){
    if(grepl("[B|b]", what[1]) & grepl("[K|k]", what[2])) dat = mxFDAobject@`Bivariate Summaries`$Kcross
    if(grepl("[B|b]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`Bivariate Summaries`$Gcross
    if(grepl("[B|b]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobject@`Bivariate Summaries`$Lcross
    if(grepl("[U|u]", what[1]) & grepl("[K|k]", what[2])) dat = mxFDAobject@`Univariate Summaries`$Kest
    if(grepl("[U|u]", what[1]) & grepl("[G|g]", what[2])) dat = mxFDAobject@`Univariate Summaries`$Gest
    if(grepl("[U|u]", what[1]) & grepl("[L|l]", what[2])) dat = mxFDAobjectx@`Univariate Summaries`$Lest
  }
  return(dat)
}
