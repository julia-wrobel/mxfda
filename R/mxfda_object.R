#' Make mxFDA class object
#'
#' @param metadata metadata with column `key` for samples
#' @param spatial spatial information, either list or df, with column `key`
#' @param key column name to link metadata to the spatial data
#'
#' @return object of class mxFDA
#' @export
#'
make_mxfda = function(metadata,
                      spatial,
                      key){
  #check that input data is in right format
  if(!inherits(clinical, "data.frame"))
    stop("Clinical as a data frame")
  if(!inherits(spatial, "list"))
    spatial = do.call(dplyr::bind_rows, spatial)
  if(!(key %in% colnames(metadata) & key %in% colnames(spatial)))
    stop("Key doesn't connect metadata and spatial data")

  datClass = methods::new("mxFDA",
                          Metadata = metadata,
                          Spatial = spatial,
                          key = key)
  return(datClass)
}

setClass("mxFDA", slots = list(
  Metadata = "data.frame",
  Spatial = "data.frame",
  key = "character",
  `Univariate Summaries` = "list", #G, K, L
  `Bivariate Summaries` = "list", #G, K, L
  `Functional Data` = "list" #maybe subset for FDA
))
