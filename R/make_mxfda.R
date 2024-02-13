#' Make mxFDA class object
#'
#' Used to create an object of class `mxFDA` that can be used with the [mxfda] package for functional data analysis.
#'
#' @param metadata metadata with columns `subject_key` and `sample_key`
#' @param spatial spatial information, either list or df, with column `sample_key`. `Spatial` can be empty if inputting data already derived. See [add_summary_function()] for more details.
#' @param subject_key column name in `Metadata` for subject ID
#' @param sample_key column linking `Metadata` to `Spatial` data
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @return S4 object of class `mxFDA`
#' \item{Metadata}{slot of class `data.frame` that contains sample and subject level information}
#' \item{Spatial}{slot of class `data.frame` that contains point level information within samples. An example would be cells belonging to TMA cores}
#' \item{subject_key}{slot of class `character` that corresponds to a column in the `Metadata` slot that groups samples at a subject level. An example would be "*patient_id*"}
#' \item{sample_key}{slot of class `character` that corresponds to a column both in the `Metadata` and `Spatial` slots that links samples to characteristics}
#' \item{univariate_summaries}{slot of class `list` where univariate summary functions calculated on `Spatial` would be stored}
#' \item{bivariate_summaries}{slot of class `list` where bivariate summary functions calculated on `Spatial` would be stored}
#' \item{functional_pca}{slot of class `list` where FPCA results are stored}
#' \item{functional_mpca}{slot of class `list` where MFPCA results are stored}
#' \item{functional_cox}{slot of class `list` where functional cox model results are stored}
#' \item{functional_mcox}{slot of class `list` where mixed functional cox model results are stored}
#' \item{scalar_on_function}{slot of class `list` where functional models are fit to scalar responses}
#'
#' @author Alex Soupir \email{`r alexsoupir_email`}
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
  `functional_mcox` = "list", #one for mixed effects (?)
  `scalar_on_functional` = "list"
))
