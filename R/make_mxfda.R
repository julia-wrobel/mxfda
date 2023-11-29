#' Make mxFDA class object
#'
#' Used to create an object of class `mxFDA` that can be used with the [mxfda] package for functional data analysis.
#'
#' @param metadata metadata with column `key` for samples
#' @param spatial spatial information, either list or df, with column `key`. `spatial` can be empty if inputting data already derived
#' @param subject_key column name for subject ID
#' @param sample_key column linking metadata to spatial data
#'
#' @return object of class mxFDA
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @author Alex Soupir \email{alex.soupir@@moffitt.org}
#'
#' @examples
#' #select sample metadata
#' clinical = lung_df %>%
#'   dplyr::select(image_id, patient_id, patientImage_id,
#'                 gender, age, survival_days, survival_status, stage) %>%
#'   dplyr::distinct()
#' #select the spatial information
#' spatial = lung_df %>%
#'   dplyr::select(-image_id, -gender, -age, -survival_days, -survival_status, -stage)
#' sample_id_column = "patientImage_id"
#' #create the mxFDA object
#' mxFDAobject = make_mxfda(metadata = clinical,
#'                          spatial = spatial,
#'                          subject_key = "patient_id",
#'                          sample_key = sample_id_column)
#'
#' @export
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

methods::setClass("mxFDA", slots = list(
  Metadata = "data.frame",
  Spatial = "data.frame",
  subject_key = "character",
  sample_key = "character",
  `univariate_summaries` = "list", #G, K, L
  `bivariate_summaries` = "list", #G, K, L
  `functional_pca` = "list", #one for 'regular'
  `functional_mpca` = "list", #one for mixed effects (?)
  `functional_cox` = "list",
  `functional_mcox` = "list" #one for mixed effects (?)
))
