#' run_fpca
#'
#' This is a wrapper for the function \code{fpca.face} from the \code{refund} package. EXPAND
#'
#'
#' @param mxFDAobject object of class \code{mxFDA} created by `make_mxfda` with metrics derived with `extract_summary_functions`
#' @param metric name of calculated spatial metric to use
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param knots Number of knots for defining spline basis.Defaults to the number of measurements per function divided by 2.
#' @param analysis_vars Optional list of variables to be retained for downstream analysis.
#' @param lightweight Default is FALSE. If TRUE, removes Y and Yhat from returned FPCA object. A good option to select for large datasets.
#' @param filter_cols a named vector of factors to filter summary functions to in `c(Derived_Column = "Level_to_Filter")` format
#' @param ... Optional other arguments to be passed to \code{fpca.face}
#'
#' @details
#' `r lifecycle::badge('stable')`
#'
#' The `filter_cols` parameter is useful when the summary function was input by the user using [add_summary_function()] and the multiple marks were assessed; a column called "Markers" with tumor infiltrating lymphocytes as well as cytotoxic T cells. This parameter allows for filtering down to include only one or the other.
#'
#' @return A \code{mxFDA} object with the `functional_pca` slot filled for the respective spatial summary function containing:
#' \item{mxfundata}{The original dataframe of spatial summary functions, with scores from FPCA appended for downstream modeling}
#' \item{fpc_object}{A list of class "fpca" with elements described in the documentation for  \code{refund::fpca.face}}
#'
#' @author Julia Wrobel \email{`r juliawrobel_email`}
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @references Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
#' Fast covariance estimation for high-dimensional functional data.
#' \emph{Statistics and Computing}, 26, 409-421.
#' DOI: 10.1007/s11222-014-9485-x.
#'
#' @importFrom refund fpca.face
#'
#' @examples
#' #load ovarian mxFDA object
#' data('ovarian_FDA')
#'
#' #run the FPCA
#' ovarian_FDA = run_fpca(ovarian_FDA, metric = "uni g", r = "r", value = "fundiff",
#'                        lightweight = TRUE,
#'                        pve = .99)
#'
#' @export
run_fpca = function(mxFDAobject,
                    metric = "uni k",
                    r = "r",
                    value = "fundiff",
                    knots = NULL,
                    analysis_vars = NULL,
                    lightweight = FALSE,
                    filter_cols = NULL,
                    ...){
  #get the right data
  if(length(metric) != 1) stop("Please provide a single spatial metric to calculate functional PCA with")
  metric = unlist(strsplit(metric, split = " "))
  #check for slot in summaries
  metric.exists(mxFDAobject, metric)
  #get data
  mxfundata = get_data(mxFDAobject, metric, 'summaries') %>%
    filter_data(filter_cols)

  computed_vals = mxfundata %>%
    dplyr::group_by(dplyr::across(!!mxFDAobject@sample_key)) %>%
    dplyr::summarise(number_computable = sum(!is.na(get(value)))) %>% #number of radii with value
    dplyr::mutate(Keep = ifelse(number_computable < 4, FALSE, TRUE),
                  Keep = factor(Keep, levels = c(TRUE, FALSE))) #true means keep %>%
  cvs = table(computed_vals$Keep) %>% data.frame() #calculated values summed
  #let user know what is kept/removed
  message(paste0(cvs[cvs$Var1 == "TRUE", "Freq"],
                 " sample have >= 4 values for FPCA; removing ",
                 cvs[cvs$Var1 == "FALSE", "Freq"], " samples"))

  #make sure that the selected columns are present to be used for fpca
  if(!(r %in% colnames(mxfundata)))
    stop("'r' not in summary data")
  if(!(value %in% colnames(mxfundata)))
    stop("'value' not in summary data")

  index_range <- range(mxfundata[[r]], na.rm = TRUE)

  mxfundata <- mxfundata %>%
    dplyr::filter(get(mxFDAobject@sample_key) %in% #filter out the samples that don't have enough values
                    (computed_vals %>%
                       dplyr::filter(Keep == "TRUE") %>%
                       pull(!!mxFDAobject@sample_key))) %>%
    dplyr::select(dplyr::all_of(c(mxFDAobject@sample_key, r, value))) %>%
    tidyr::pivot_wider(names_from =  dplyr::all_of(r),
                names_prefix = "r_",
                values_from =  dplyr::all_of(value))

  mat <- mxfundata %>%
    dplyr::select(dplyr::starts_with("r_")) %>%
    as.matrix()

  if(is.null(knots)) knots = floor(ncol(mat)/3)
  if(knots > ncol(mat) - 5) knots = floor(ncol(mat)/3)

  # run fpca
  mx_fpc <- refund::fpca.face(Y = mat, knots = knots, ...)

  if(lightweight){
    mx_fpc$Y <- NULL
    mx_fpc$Yhat <- NULL
  }
  mx_fpc$index_range <- index_range

  score_df <- stats::setNames(as.data.frame(mx_fpc$scores), paste0("fpc", 1:mx_fpc$npc))

  # append all FPCA scores to dataframe that has one row per subject, then convert to long format
  mxfundata = dplyr::bind_cols(mxfundata, score_df) %>%
    dplyr::select(-dplyr::starts_with("r_"))

  fpca_dat <- list(score_df = mxfundata,
                   fpc_object = mx_fpc)

  if(grepl("[B|b]", metric[1]) & grepl("[K|k]", metric[2])) mxFDAobject@`functional_pca`$Kcross = fpca_dat
  if(grepl("[B|b]", metric[1]) & grepl("[G|g]", metric[2])) mxFDAobject@`functional_pca`$Gcross = fpca_dat
  if(grepl("[B|b]", metric[1]) & grepl("[L|l]", metric[2])) mxFDAobject@`functional_pca`$Lcross = fpca_dat
  if(grepl("[U|u]", metric[1]) & grepl("[K|k]", metric[2])) mxFDAobject@`functional_pca`$Kest = fpca_dat
  if(grepl("[U|u]", metric[1]) & grepl("[G|g]", metric[2])) mxFDAobject@`functional_pca`$Gest = fpca_dat
  if(grepl("[U|u]", metric[1]) & grepl("[L|l]", metric[2])) mxFDAobject@`functional_pca`$Lest = fpca_dat

  return(mxFDAobject)

}
