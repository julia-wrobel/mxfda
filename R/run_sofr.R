#' Run Scalar on Function Regression
#'
#' Fit a scalar-on-function regression model. Uses refund::pfr under the hood for computations, and stores results in the mxfda object.
#'
#' @param mxFDAobject Dataframe of spatial summary functions from multiplex imaging data, in long format. Can be estimated using the function \code{extract_summary_functions} or provided separately.
#' @param model_name character string to give the fit model
#' @param formula Formula to be fed to mgcv in the form of outcome ~ x1 + x2. Does not contain functional predictor. Character valued.
#' @param family Exponential family distribution to be passed to \code{mgcv::gam}. Defaults to "gaussian". Select "binomial" for binary outcome.
#' @param metric Name of calculated spatial metric to use
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param knots Number of knots for defining spline basis.
#' @param smooth Option to smooth data using FPCA. Defaults to FALSE.
#' @param filter_cols a named vector of factors to filter summary functions to in `c(Derived_Column = "Level_to_Filter")` format
#' @param ... Optional other arguments to be passed to \code{fpca.face}
#'
#' @importFrom refund pfr lf
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @return A \code{list} which is a linear or additive functional Cox model fit. See \code{mgcv::gam} for more details.
#'
#' @author Julia Wrobel \email{`r juliawrobel_email`}
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @examples
#' #load ovarian mxFDA object
#' data('ovarian_FDA')
#'
#' # run scalar on function regression model with a continuous outcome (age)
#' ovarian_FDA = run_sofr(ovarian_FDA,
#'                        model_name = "fit_sofr",
#'                        formula = age~stage,
#'                        metric = "uni g", r = "r", value = "fundiff")
#'
#' # run scalar on function regression model with a binary outcome (stage)
#' # also known as functional logistic regression
#' ovarian_FDA = run_sofr(ovarian_FDA,
#'                        model_name = "fit_sofr",
#'                        formula = stage~age,
#'                        family = "binomial",
#'                        metric = "uni g", r = "r", value = "fundiff")
#'
#' @export
run_sofr <- function(mxFDAobject,
                    model_name,
                    formula,
                    family = "gaussian",
                    metric = "uni k",
                    r = "r",
                    value = "fundiff",
                    smooth = FALSE,
                    filter_cols = NULL,
                    ...,
                    knots = NULL){
  #get the right data from the object
  if(length(metric) != 1)
    stop("Please provide a single spatial metric to calculate functional cox models with")
  metric = unlist(strsplit(metric, split = " "))
  #check for slot in summaries
  metric.exists(mxFDAobject, metric)

  mxfundata = get_data(mxFDAobject, metric, 'summaries') %>%
    filter_data(filter_cols) %>%
    dplyr::full_join(mxFDAobject@Metadata, by = mxFDAobject@sample_key)

  computed_vals = mxfundata %>%
    dplyr::group_by(dplyr::across(!!mxFDAobject@sample_key)) %>%
    dplyr::summarise(number_computable = sum(!is.na(get(value)))) %>%
    dplyr::mutate(Keep = ifelse(number_computable < 4, FALSE, TRUE),
                  Keep = factor(Keep, levels = c(TRUE, FALSE)))
  cvs = table(computed_vals$Keep) %>% data.frame()
  message(paste0(cvs[cvs$Var1 == "TRUE", "Freq"],
                 " sample have >= 4 values for FPCA; removing ",
                 cvs[cvs$Var1 == "FALSE", "Freq"], " samples"))
  #remove bad samples
  mxfundata = mxfundata %>%
    filter(get(mxFDAobject@sample_key) %in%
             computed_vals[[mxFDAobject@sample_key]][computed_vals$Keep == TRUE])

  analysis_vars = unique(c(all.vars(formula)))

  # check for missing values in the functional predictor- these must be removed
  if(smooth){
    mxfundata <- impute_fpca(mxfundata, id = mxFDAobject@sample_key, r = r, value = value,
                             analysis_vars = analysis_vars, smooth = TRUE)
    message("Functional predictor contains NA values that were imputed using FPCA")
  }
  if(anyNA(mxfundata[[value]])){
    mxfundata <- impute_fpca(mxfundata, id = mxFDAobject@sample_key , r = r, value = value, knots = knots,
                               analysis_vars = analysis_vars, smooth)
    message("Functional predictor contains NA values that were imputed using FPCA")
  }

  mxfundata <- mxfundata %>%
    select(dplyr::all_of(c(mxFDAobject@sample_key, r, value, analysis_vars))) %>%
    tidyr::pivot_wider(names_from =  dplyr::all_of(r),
                       names_prefix = "r_",
                       values_from =  dplyr::all_of(value))

  mat <- mxfundata %>%
    dplyr::select(dplyr::starts_with("r_")) %>%
    as.matrix()

  mxfundata = mxfundata %>% select(-contains("r_"))
  mxfundata$xmat = mat

  if(is.null(knots)) knots = floor(ncol(mat)/3)
  if(knots > ncol(mat) - 5) knots = floor(ncol(mat)/3)
  if(knots < 3) knots = 5

  if(family == "binomial"){
     outcome = stats::formula(formula)[[2]]
     mxfundata[[outcome]] = as.numeric(mxfundata[[outcome]])
  }

  form = deparse(stats::formula(formula))
  form =  paste0(form, ' + lf(xmat, k=', knots, ')')

  fit_sofr <- pfr(formula = stats::as.formula(form),
                  #family = family,
                  data = mxfundata)

  class(fit_sofr) <- append("sofr", class(fit_sofr))

  if(grepl("[B|b]", metric[1]) & grepl("[K|k]", metric[2])) mxFDAobject@`scalar_on_functional`$Kcross[[model_name]] = fit_sofr
  if(grepl("[B|b]", metric[1]) & grepl("[G|g]", metric[2])) mxFDAobject@`scalar_on_functional`$Gcross[[model_name]] = fit_sofr
  if(grepl("[B|b]", metric[1]) & grepl("[L|l]", metric[2])) mxFDAobject@`scalar_on_functional`$Lcross[[model_name]] = fit_sofr
  if(grepl("[U|u]", metric[1]) & grepl("[K|k]", metric[2])) mxFDAobject@`scalar_on_functional`$Kest[[model_name]] = fit_sofr
  if(grepl("[U|u]", metric[1]) & grepl("[G|g]", metric[2])) mxFDAobject@`scalar_on_functional`$Gest[[model_name]] = fit_sofr
  if(grepl("[U|u]", metric[1]) & grepl("[L|l]", metric[2])) mxFDAobject@`scalar_on_functional`$Lest[[model_name]] = fit_sofr

 return(mxFDAobject)
}

