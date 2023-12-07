#' Filter Spatial data
#'
#' function to filter the spatial data slot of the `mxFDA` object.
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param ... expressions that return a logical TRUE/FALSE value when evaluated on columns of the spatial data slot. These expressions get passed to [dplyr::filter()] so must be compatible.
#' @param force logical whether or not to return empty spatial data *if* filtering results in 0 rows
#'
#' @returns object of class `mxFDA` with the spatial slot filtered. See [make_mxfda()] for more details on object
#'
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @references [dplyr::filter()]
#'
#' @examples
#' #load ovarian mxFDA object
#' data(ovarian_FDA)
#'
#' #filter ages greater than 50
#' ovarian_FDA_age50 = filter_spatial(ovarian_FDA, age >= 50)
#'
#' @export
#'
filter_spatial = function(mxFDAobject, ..., force = FALSE){
  d = dplyr::quos(...)
  spat = mxFDAobject@Spatial
  nspat = spat %>%
    dplyr::filter(!!!d)
  if(!force & nrow(nspat) == 0)
    stop("No rows left in the spatial data after filtering")
  mxFDAobject@Spatial = nspat
  return(mxFDAobject)
}
