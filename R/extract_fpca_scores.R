#' Extract FPCA scores
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param what what functional PCA data to extract, e.g. 'uni k'
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
extract_fpca_scores = function(mxFDAobject, what){
  if(!inherits(mxFDAobject, "mxFDA")) stop("supply object of class `mxFDA`")
  #check if object is of class mxFDA
  what = unlist(strsplit(what, split = " "))

  dat = get_data(mxFDAobject, what[1:2], type = what[3])

  if(grepl("m", what[3]))
    return(dat[1:2])

  #combine metadata with the scores
  #need to have the same number of metadata as number of derived samples
  #pretty sure if the number of samples in metadata don't match the number of spatial samples this will not work
  #no alternatiave because score_df doesn't report sample ID with it.
  new_df = mxFDAobject@Metadata %>%
    dplyr::bind_cols(dat$score_df)
  return(new_df)
}
