#' Multiplex imaging data from a non-small cell lung cancer study.
#'
#' This data is adapted from the VectraPolarisData Bioconductor package. There are multiple ROIs for each patient.
#'
#' @format ## `lung_df`
#' A data frame with 879,694 rows and 19 columns:
#' \describe{
#'   \item{image_id}{Image id for a given patient}
#'   \item{patient_id}{Unique patient id}
#'   \item{age}{Patient age at time of cancer diagnosis}
#'   \item{survival_days}{Survival time from diagnosis, in days}
#'   \item{survival_status}{Censoring variable, 1 = death, 0 = censor}
#'   \item{x}{Cell x position}
#'   \item{y}{Cell y position}
#'  ...
#' }
#' @source <https://bioconductor.org/packages/release/data/experiment/html/VectraPolarisData.html>
"lung_df"

#' Multiplex imaging data from an ovarian cancer tumor microarray
#'
#' This data is adapted from the VectraPolarisData Bioconductor package and comes from a tumor-microarray of tissue samples from 128 patients with ovarian cancer.
#' There is one patient per subject.
#'

#' Spatial summary functions of ovarian cancer multiplex imaging data.
#'
#' This data is adapted from the VectraPolarisData Bioconductor package.
#' Signal between the survival outcome and spatial summary functions has been augmented for teaching purposes.
#' Spatial relationship is summarized using the nearest neighbor G function.
#'
#' @format ## `ovarian_FDA`
#' An mxFDA object with augmented ovarian cancer multiplex immunofluorescence data, and NN G(r) calculated:
#' \describe{
#'   \item{Metadata}{information about the spatial samples with column `sample_key` column in both}
#'   \item{Spatial}{cell-level information with `x` and `y` columns along with `sample_key` to link to `Metadata`}
#'   \item{subject_key}{column in `Metadata` that may have multiple `sample_key` values for each, akin to patient IDs}
#'   \item{sample_key}{column in both `Metadata` and `Spatial` that is a 1:1 with the samples (unique per sample)}
#'   \item{univariate_summaries}{univariate summary slot with nearest neighbor G calculared}
#'   \item{bivariate_summaries}{empty slot available for bivariate summaries}
#'   \item{functional_pca}{empty slot for functional PCA data of summaries}
#'   \item{functional_cox}{empty slot for functional models}
#' }
#' @source <https://bioconductor.org/packages/release/data/experiment/html/VectraPolarisData.html>
"ovarian_FDA"

#' Multiplex imaging data from a non-small cell lung cancer study
#'
#' This data is adapted from the VectraPolarisData Bioconductor package. There are multiple ROIs for each patient.
#'
#' Spatial summary functions of lung cancer multiplex imaging data.
#'
#' This data is adapted from the VectraPolarisData Bioconductor package.
#' Signal between the survival outcome and spatial summary functions has been augmented for teaching purposes.
#' Spatial relationship is summarized using the nearest neighbor G function.
#'
#' Includes only spatial samples that had 10 or more radii with calculable G function
#'
#' @format ## `lung_FDA`
#' An mxFDA object with augmented non-small cel lung cancer multiplex immunofluorescence data, and NN G(r) calculated:
#' \describe{
#'   \item{Metadata}{information about the spatial samples with column `sample_key` column in both}
#'   \item{Spatial}{cell-level information with `x` and `y` columns along with `sample_key` to link to `Metadata`}
#'   \item{subject_key}{column in `Metadata` that may have multiple `sample_key` values for each, akin to patient IDs}
#'   \item{sample_key}{column in both `Metadata` and `Spatial` that is a 1:1 with the samples (unique per sample)}
#'   \item{univariate_summaries}{univariate summary slot with nearest neighbor G calculared}
#'   \item{bivariate_summaries}{empty slot available for bivariate summaries}
#'   \item{functional_pca}{empty slot for functional PCA data of summaries}
#'   \item{functional_cox}{empty slot for functional models}
#' }
#' @source <https://bioconductor.org/packages/release/data/experiment/html/VectraPolarisData.html>
"lung_FDA"
