#' Extract FPCA object
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param what what functional PCA data to extract, e.g. 'uni k'
#' @param r column id for radii
#' @param value value used to calculate fpcs
#'
#' @return fpca object
#'
#' @author Alex Soupir \email{alex.soupir@@moffitt.org}
#'
#' @examples
#' #set seed
#' set.seed(333)
#'
#' @export
extract_fpca_object = function(mxFDAobject, what, r, value){
  if(!inherits(mxFDAobject, "mxFDA")) stop("supply object of class `mxFDA`")
  #check if object is of class mxFDA
  what = unlist(strsplit(what, split = " "))

  dat = get_data(mxFDAobject, what, type = "fpca")

  return(dat$fpc_object)
}
