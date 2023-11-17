#' Make mxFDA class object
#'
#' @param metadata metadata with column `key` for samples
#' @param spatial spatial information, either list or df, with column `key`. `spatial` can be empty if inputting data already derived
#' @param subject_key column name for subject ID
#' @param sample_key column linking metadata to spatial data
#'
#' @return object of class mxFDA
#' @export
#'
#' @author Alex Soupir \email{alex.soupir@@moffitt.org}
#'
#' @examples
#' #set seed
#' set.seed(333)
#'
make_mxfda = function(metadata,
                      spatial = NULL,
                      subject_key,
                      sample_key){
  #check that input data is in right format
  if(!inherits(metadata, "data.frame"))
    stop("Clinical as a data frame")
  if(inherits(spatial, "list"))
    spatial = do.call(dplyr::bind_rows, spatial)
  #make dummy spatial if null
  if(is.null(spatial)){
    spatial = data.frame()
  } else if(!(sample_key %in% colnames(metadata) & sample_key %in% colnames(spatial))){
    stop("Sample key doesn't connect metadata and spatial data")
  }

  if(!(subject_key %in% colnames(metadata)))
    stop("Subject key not in metadata")

  datClass = methods::new("mxFDA",
                          Metadata = metadata,
                          Spatial = spatial,
                          subject_key = subject_key,
                          sample_key = sample_key)
  return(datClass)
}

setClass("mxFDA", slots = list(
  Metadata = "data.frame",
  Spatial = "data.frame",
  subject_key = "character",
  sample_key = "character",
  `univariate_summaries` = "list", #G, K, L
  `bivariate_summaries` = "list", #G, K, L
  `functional_pca` = "list", #maybe subset for FDA
  `functional_cox` = "list"
))
