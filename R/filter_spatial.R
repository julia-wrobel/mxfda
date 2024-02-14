#' Filter Spatial data
#'
#' function to filter the spatial data slot of the `mxFDA` object.
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param ... expressions that return a logical TRUE/FALSE value when evaluated on columns of the meta data slot. These expressions get passed to [dplyr::filter()] so must be compatible.
#' @param based_on character for which data slot to use for filtering, either 'meta', or 'spatial'. Default to 'meta'.
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
#' ovarian_FDA_age50 = filter_spatial(ovarian_FDA, age >= 50, based_on = 'meta')
#'
#' @export
#'
filter_spatial = function(mxFDAobject, ..., based_on = 'meta', force = FALSE){
  d = dplyr::quos(...)
  spat = mxFDAobject@Spatial
  #filter spatial data by samples with meta filtering
  if(based_on == "meta"){
    clin = mxFDAobject@Metadata %>%
      dplyr::filter(!!!d)
    nspat = spat %>%
      dplyr::filter(get(mxFDAobject@sample_key) %in% clin[[mxFDAobject@sample_key]])
  } else if(based_on == "spatial"){
    nspat = spat %>%
      dplyr::filter(!!!d)
  } else {
    stop("`based_on` needs to be 'meta' or 'spatial'")
  }
  if(!force & nrow(nspat) == 0)
    stop("No rows left in the spatial data after filtering")
  mxFDAobject@Spatial = nspat
  return(mxFDAobject)
}

