#' run_fpca
#'
#' This is a wrapper for the function \code{mfpca.face} from the \code{refund} package. EXPAND
#'
#'
#' @param mxFDAobject object of class \code{mxFDA} created by [make_mxfda()] with metrics derived with [extract_summary_functions()]
#' @param metric name of calculated spatial metric to use
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param knots Number of knots for defining spline basis.Defaults to the number of measurements per function divided by 2.
#' @param lightweight Default is FALSE. If TRUE, removes Y and Yhat from returned mFPCA object. A good option to select for large datasets.
#' @param ... Optional other arguments to be passed to \code{mfpca.face}
#' @param twoway whether to model within within sample variability (such as )
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @return A \code{mxFDA} object  with the `functional_mpca` slot for the respective spatial summary function containing:
#' \item{mxfundata}{The original dataframe of spatial summary functions, with scores from FPCA appended for downstream modeling}
#' \item{fpc_object}{A list of class "fpca" with elements described in the documentation for  \code{refund::fpca.face}}
#'
#' @author unknown \email{first.last@@domain.extension}
#' @author Julia Wrobel \email{`r juliawrobel_email`}
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @references Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
#' Fast covariance estimation for high-dimensional functional data.
#' \emph{Statistics and Computing}, 26, 409-421.
#' DOI: 10.1007/s11222-014-9485-x.
#'
#' @importFrom refund mfpca.face
#' @importFrom tidyr pivot_wider
#' @importFrom rlang `:=`
#' @import dplyr
#'
#' @examples
#' #load data
#' data(lung_FDA)
#'
#' #run mixed fpca
#' lung_FDA = run_mfpca(lung_FDA, metric = 'uni g')
#'
#' @export
run_mfpca = function(mxFDAobject,
                    metric = "uni k",
                    r = "r",
                    value = "fundiff",
                    knots = NULL,
                    lightweight = FALSE,
                    ...,
                    twoway = FALSE){
  #get the right data
  if(length(metric) != 1) stop("Please provide a single spatial metric to calculate functional PCA with")
  metric = unlist(strsplit(metric, split = " "))
  #check for slot in summaries
  metric.exists(mxFDAobject, metric)
  #get data
  mxfundata = get_data(mxFDAobject, metric, 'summaries')

  #add subject column
  mxfundata = mxFDAobject@Metadata %>%
    dplyr::select(dplyr::all_of(c(mxFDAobject@subject_key, mxFDAobject@sample_key))) %>%
    dplyr::right_join(mxfundata, by = mxFDAobject@sample_key)

  # NEED TO CHANGE THIS BEHAVIOR- WHAT DO WE DO IF WE HAVE AN NA FOR THE IMAGE SUMMARY FUNCTION
  mxfundata = mxfundata %>% filter(!is.na(r))
  index_range <- range(mxfundata[[r]])

  # this seems to break when there are NA values, what behavior do I want for that?
  mxfundata <- mxfundata %>%
    select(all_of(c(mxFDAobject@subject_key, mxFDAobject@sample_key, r, value))) %>%
    pivot_wider(names_from = r,
                names_prefix = "r_",
                values_from = value) %>%
    arrange(!!mxFDAobject@subject_key, !!mxFDAobject@sample_key)

  mat <- mxfundata %>%
    select(starts_with("r_")) %>%
    as.matrix()

  if(is.null(knots)) knots = floor(ncol(mat)/3)
  if(knots > ncol(mat) - 5) knots = floor(ncol(mat)/3)

  # run fpca
  mx_mfpc <- mfpca.face(Y = mat, id = mxfundata[[mxFDAobject@subject_key]],
                        visit = mxfundata[[mxFDAobject@sample_key]],
                        knots = knots, ...,
                        twoway = FALSE) # argument, but default is false

  if(lightweight){
    mx_mfpc$Y <- NULL
    mx_mfpc$Xhat <- NULL
  }
  mx_mfpc$index_range <- index_range

  # how do we want to return the scores? Right now returning the level 1 scores and the standard deviation of the level 2 scores
  score_df <- stats::setNames(as.data.frame(mx_mfpc$scores$level1), paste0("level1_", 1:mx_mfpc$npc$level1))
  score_df[[mxFDAobject@subject_key]] <- unique(mxfundata[[mxFDAobject@subject_key]])

  scores_level2 <- stats::setNames(as.data.frame(mx_mfpc$scores$level2), paste0("level2_", 1:mx_mfpc$npc$level2))
  scores_level2 <- scores_level2 %>%
    mutate(!!mxFDAobject@subject_key := mxfundata[[mxFDAobject@subject_key]],
           !!mxFDAobject@sample_key := mxfundata[[mxFDAobject@sample_key]])



  fpca_dat = list(score_df = score_df,
                  scores_level2 = scores_level2,
                  mfpc_object = mx_mfpc)

  if(grepl("[B|b]", metric[1]) & grepl("[K|k]", metric[2])) mxFDAobject@`functional_mpca`$Kcross = fpca_dat
  if(grepl("[B|b]", metric[1]) & grepl("[G|g]", metric[2])) mxFDAobject@`functional_mpca`$Gcross = fpca_dat
  if(grepl("[B|b]", metric[1]) & grepl("[L|l]", metric[2])) mxFDAobject@`functional_mpca`$Lcross = fpca_dat
  if(grepl("[U|u]", metric[1]) & grepl("[K|k]", metric[2])) mxFDAobject@`functional_mpca`$Kest = fpca_dat
  if(grepl("[U|u]", metric[1]) & grepl("[G|g]", metric[2])) mxFDAobject@`functional_mpca`$Gest = fpca_dat
  if(grepl("[U|u]", metric[1]) & grepl("[L|l]", metric[2])) mxFDAobject@`functional_mpca`$Lest = fpca_dat

  return(mxFDAobject)

}
