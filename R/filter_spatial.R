#' Filter Spatial data
#'
#' function to filter the spatial data slot of the `mxFDA` object.
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param ... expressions that return a logical TRUE/FALSE value when evaluated on columns of the spatial data slot. These expressions get passed to [dplyr::filter()] so must be compatible.
#'
#' @returns mxFDAobject with updated spatial slot
#'
#' @author Alex Soupir \email{alex.soupir@@moffitt.org}
#'
#' @references [dplyr::filter()]
#'
#' @examples
#' # example code
#' set.seed(333)
#'
#' @export
#'
filter_spatial = function(mxFDAobject, ...){
  d = dplyr::quos(...)
  spat = mxFDAobject@Spatial
  nspat = spat %>%
    dplyr::filter(!!!d)
  mxFDAobject@Spatial = nspat
  return(mxFDAobject)
}
