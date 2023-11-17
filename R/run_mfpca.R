#' run_fpca
#'
#' This is a wrapper for the function \code{mfpca.face} from the \code{refund} package. EXPAND
#'
#'
#' @param mxFDAobject object of class \code{mxFDA} created by `make_mxfda` with metrics derived with `extract_summary_functions`
#' @param id Character string, the name of the variable that identifies each unique subject.
#' @param image_id Character string, the name of the variable that identifies each unique image sample
#' @param metric name of calculated spatial metric to use
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param knots Number of knots for defining spline basis.Defaults to the number of measurements per function divided by 2.
#' @param analysis_vars Optional list of variables to be retained for downstream analysis.
#' @param lightweight Default is FALSE. If TRUE, removes Y and Yhat from returned mFPCA object. A good option to select for large datasets.
#' @param ... Optional other arguments to be passed to \code{mfpca.face}
#'
#' @return A \code{mxFDA} object containing:
#' \item{mxfundata}{The original dataframe of spatial summary functions, with scores from FPCA appended for downstream modeling}
#' \item{fpc_object}{A list of class "fpca" with elements described in the documentation for  \code{refund::fpca.face}}
#'
#' @author unknown \email{first.last@@domain.extension}
#'
#' @references Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
#' Fast covariance estimation for high-dimensional functional data.
#' \emph{Statistics and Computing}, 26, 409-421.
#' DOI: 10.1007/s11222-014-9485-x.
#'
#' @importFrom refund mfpca.face
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_wider
#' @import dplyr
#'
#' @examples
#' # simulate data
#' set.seed(1001)
#'
#' @export
run_mfpca = function(mxFDAobject,
                    id,
                    image_id,
                    metric = "uni k",
                    r = "r",
                    value = "fundiff",
                    knots = NULL,
                    analysis_vars = NULL,
                    lightweight = FALSE,
                    ...){
  #get the right data
  if(length(metric) != 1) stop("Please provide a single spatial metric to calculate functional PCA with")
  metric = unlist(strsplit(metric, split = " "))
  mxfundata = get_data(mxFDAobject, metric, 'summaries')

  ## add stopper here to check if they provided a metric that doesn't exist
  # right now we can split the patientImage_id but usually that won't work. Need to update S4 object to take both patient and image id
  mxfundata = mxfundata %>%
    separate(patientImage_id, into = c("id", "image_id"), sep = "_")

  # NEED TO CHANGE THIS BEHAVIOR- WHAT DO WE DO IF WE HAVE AN NA FOR THE IMAGE SUMMARY FUNCTION
  mxfundata = mxfundata %>% filter(!is.na(r))
  index_range <- range(mxfundata[[r]])

  # this seems to break when there are NA values, what behavior do I want for that?
  mxfundata <- mxfundata %>%
    select(all_of(c(id, image_id, r, value))) %>%
    pivot_wider(names_from = r,
                names_prefix = "r_",
                values_from = value) %>%
    arrange(id, image_id)

  mat <- mxfundata %>%
    select(starts_with("r_")) %>%
    as.matrix()

  if(is.null(knots)) knots = floor(ncol(mat)/3)
  if(knots > ncol(mat) - 5) knots = floor(ncol(mat)/3)

  # run fpca
  mx_mfpc <- mfpca.face(Y = mat, id = mxfundata[[id]],
                        visit = mxfundata[[image_id]],
                        twoway = FALSE, # argument, but default is false
                        knots = knots, ...)

  if(lightweight){
    mx_mfpc$Y <- NULL
    mx_mfpc$Xhat <- NULL
  }
  mx_mfpc$index_range <- index_range

  # how do we want to return the scores? Right now returning the level 1 scores and the standard deviation of the level 2 scores
  score_df <- setNames(as.data.frame(mx_mfpc$scores$level1), paste0("level1_", 1:mx_mfpc$npc$level1))
  score_df$id <- unique(mxfundata[[id]])

  scores_level2 <- setNames(as.data.frame(mx_mfpc$scores$level2), paste0("level2_", 1:mx_mfpc$npc$level2))
  scores_level2 <- scores_level2 %>%
    mutate(id = mxfundata[[id]],
           image_id = mxfundata[[image_id]])



  fpca_dat = list(score_df = score_df,
                  scores_level2 = scores_level2,
                  mfpc_object = mx_mfpc)

  return(fpca_dat)

}
