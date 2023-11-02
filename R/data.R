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
#' A data frame with 12,726 rows and 8 columns:
#' \describe{
#'   \item{patient_id}{Unique patient id}
#'   \item{age}{Patient age at time of cancer diagnosis}
#'   \item{survival_time}{Survival time from diagnosis, in days}
#'   \item{event}{Censoring variable, 1 = death, 0 = censor}
#'   \item{r}{radius}
#'   \item{sumfun}{Value of G function}
#'   \item{csr}{Value of G function under assumption of complete spatial randomness}
#'   \item{fundiff}{sumfun-csr}
#' }
#' @source <https://bioconductor.org/packages/release/data/experiment/html/VectraPolarisData.html>
"ovarian_FDA"
