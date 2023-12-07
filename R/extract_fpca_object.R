#' Extract FPCA object
#'
#' Function that extracts the FPCA object created either by [run_fpca()] or [run_mfpca()] from the `mxFDA` object
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param what what functional PCA data to extract, e.g. 'uni k'
#'
#' @details
#' `r lifecycle::badge('stable')`
#'
#' Output object can be visualized with [refund.shiny::plot_shiny()]
#'
#' @return `fpca` object created with [run_fcm()]
#'
#' @author Alex Soupir \email{alex.soupir@@moffitt.org}
#'
#' @examples
#' #set seed
#' set.seed(333)
#'
#' @export
extract_fpca_object = function(mxFDAobject, what){
  if(!inherits(mxFDAobject, "mxFDA")) stop("supply object of class `mxFDA`")
  #check if object is of class mxFDA
  what = unlist(strsplit(what, split = " "))
  if(length(what) < 3) stop("need to specify which FPCA")

  dat = get_data(mxFDAobject, what[1:2], type = what[3])

  if(grepl("m", what[3]))
    return(dat$mfpc_object)
  return(dat$fpc_object)
}
