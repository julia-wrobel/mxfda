#' run_fpca
#'
#' This is a wrapper for the function \code{fpca.face} from the \code{refund} package. EXPAND
#'
#' @importFrom refund fpca.face
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_wider
#' @import dplyr
#'
#' @return A \code{list} containing:
#' \item{mxfundata}{The original dataframe of spatial summary functions, with scores from FPCA appended for downstream modeling}
#' \item{fpc_object}{A list of class "fpca" with elements described in the documentation for  \code{refund::fpca.face}}
#'
#' @references Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
#' Fast covariance estimation for high-dimensional functional data.
#' \emph{Statistics and Computing}, 26, 409-421.
#' DOI: 10.1007/s11222-014-9485-x.

#' @examples
#' # simulate data
#' set.seed(1001)
#'
#'
#' @param mxfundata Dataframe of spatial summary functions from multiplex imaging data, in long format. Can be estimated using the function \code{extract_summary_functions} or provided separately. For FPCA there should be one function per image, and only one image per subject.
#' @param id Character string, the name of the variable that identifies each unique subject.
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param knots Number of knots for defining spline basis.Defaults to the number of measurements per function divided by 2.
#' @param analysis_vars Optional list of variables to be retained for downstream analysis.
#' @param lightweight Default is FALSE. If TRUE, removes Y and Yhat from returned FPCA object. A good option to select for large datasets.
#' @param ... Optional other arguments to be passed to \code{fpca.face}
#' @export
run_fpca = function(mxfundata,
                    id,
                    r = "r",
                    value = "fundiff",
                    knots = NULL,
                    analysis_vars = NULL,
                    lightweight = FALSE,
                    ...){

  index_range <- range(mxfundata[[r]])

  mxfundata <- mxfundata %>%
    select(all_of(c(id, r, value, analysis_vars))) %>%
    pivot_wider(names_from = r,
                names_glue = "r_{round(r, 2)}",
                values_from = value)

  mat <- mxfundata %>%
    select(starts_with("r_")) %>%
    as.matrix()

  if(is.null(knots)) knots = floor(ncol(mat)/3)
  if(knots > ncol(mat) - 5) knots = floor(ncol(mat)/3)

  # run fpca
  mx_fpc <- fpca.face(Y = mat, knots = knots, ...)

  if(lightweight){
    mx_fpc$Y <- NULL
    mx_fpc$Yhat <- NULL
  }
  mx_fpc$index_range <- index_range

  score_df = setNames(as.data.frame(mx_fpc$scores), paste0("fpc_", 1:mx_fpc$npc))
  # append all FPCA scores to dataframe that has one row per subject
  mxfundata = bind_cols(mxfundata, score_df)

  list(mxfundata = mxfundata,
       fpc_object = mx_fpc)

}
