#' Extract FPCA object
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param what what functional PCA data to extract, e.g. 'uni k'
#' @param id sample id column used in `run_fpca` function
#' @param r column id for radii
#' @param value value used to calculate fpcs
#'
#' @return fpca object
#' @export
#'
extract_fpca_object = function(mxFDAobject, what, id, r, value){
  if(!inherits(mxFDAobject_oneImage, "mxFDA")) stop("supply object of class `mxFDA`")
  #check if object is of class mxFDA
  what = unlist(strsplit(what, split = " "))

  dat = get_data(mxFDAobject, what, type = "functional")
  # mxfundata = get_data(mxFDAobject, what, type = "summaries") %>%
  #   select(all_of(c(id, r, value))) %>%
  #   pivot_wider(names_from = r,
  #               names_prefix = "r_",
  #               values_from = value)
  return(dat$fpc_object)

  #fpca object already has yhat and y, need to be replaced?
  # fpc_object = dat$fpc_object
  # fpc_object$Y = dplyr::bind_cols(mxfundata, dat$score_df) %>%
  #   dplyr::select(dplyr::contains("r_"))
  # Yhat_dims = dim(fpc_object$Yhat)
  # fpc_object$Yhat = matrix(rnorm(Yhat_dims[1], Yhat_dims[2]),
  #                          nrow = Yhat_dims[1],
  #                          ncol = Yhat_dims[2])
  # return(fpc_object)
}
