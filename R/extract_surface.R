#' Extract Surface
#'
#' Function that transforms functional models from linear or additive functional cox models into `afcmSurface` or `lfcmSurface` objects to be plotted.
#'
#' @param mxFDAobject object of class `mxFDA` with model `model` calculated wihtin
#' @param metric spatial summary function to extract surface for
#' @param model character string for the name of the model for `metric` data
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param grid_length Length of grid on which to evaluate coefficient functions.
#' @param analysis_vars Other variables used in modeling FCM fit.
#' @param p numeric p-value used for predicting significant AFCM surface
#' @param filter_cols a named vector of factors to filter summary functions to in `c(Derived_Column = "Level_to_Filter")` format
#'
#' @return a 4 element list of either class `lfcmSurface` or `afcmSurface` depending on the class of model
#' \item{Surface}{`data.frame` for term predictions for the surface of the metric * radius area}
#' \item{Prediction}{`data.frame` for standard error of the terms for the above surface. AFCM models use the `p` to set the upper and lower standard errors of \eqn{\beta_1}}
#' \item{Metric}{character of the spatial summary function used; helps keep track if running many models}
#' \item{P-value}{a numeric value of the input p-value}
#'
#' @author Julia Wrobel \email{`r juliawrobel_email`}
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @importFrom reshape2 melt
#' @importFrom mgcv predict.gam
#' @import dplyr
#'
#' @examples
#' #load ovarian mxFDA object
#' data('ovarian_FDA')
#'
#' #run the lfcm model
#' ovarian_FDA = run_fcm(ovarian_FDA, model_name = "fit_lfcm",
#'                       formula = survival_time ~ age, event = "event",
#'                       metric = "uni g", r = "r", value = "fundiff",
#'                       analysis_vars = c("age", "survival_time"),
#'                       afcm = FALSE)
#'
#' #extract surface
#' model_surface = extract_surface(ovarian_FDA, metric = 'uni g',
#'                                 model = 'fit_lfcm',
#'                                 analysis_vars = 'age') #variables in model
#'
#' @export
extract_surface = function(mxFDAobject,
                           metric,
                           model = NULL,
                           r = "r",
                           value = "fundiff",
                           grid_length = 100,
                           analysis_vars,
                           p = 0.05,
                           filter_cols = NULL){
  #data preparation
  if(is.null(model))
    stop("Must provide a model")
  #get the right data
  if(length(metric) != 1)
    stop("Please provide a single spatial metric to extract surface for")
  metric = unlist(strsplit(metric, split = " "))

  mxfundata = get_data(mxFDAobject, metric, 'summaries') %>%
    filter_data(filter_cols) %>%
    dplyr::full_join(mxFDAobject@Metadata, by = mxFDAobject@sample_key)
  #make sure that the selected columns are present to be used for fpca
  if(!(r %in% colnames(mxfundata)))
    stop("'r' not in summary data")
  if(!(value %in% colnames(mxfundata)))
    stop("'value' not in summary data")

  #get the model
  if(grepl("[B|b]", metric[1]) & grepl("[K|k]", metric[2])) fit = mxFDAobject@functional_cox$Kcross[[model]]
  if(grepl("[B|b]", metric[1]) & grepl("[G|g]", metric[2])) fit = mxFDAobject@functional_cox$Gcross[[model]]
  if(grepl("[B|b]", metric[1]) & grepl("[L|l]", metric[2])) fit = mxFDAobject@functional_cox$Lcross[[model]]
  if(grepl("[U|u]", metric[1]) & grepl("[K|k]", metric[2])) fit = mxFDAobject@functional_cox$Kest[[model]]
  if(grepl("[U|u]", metric[1]) & grepl("[G|g]", metric[2])) fit = mxFDAobject@functional_cox$Gest[[model]]
  if(grepl("[U|u]", metric[1]) & grepl("[L|l]", metric[2])) fit = mxFDAobject@functional_cox$Lest[[model]]
  if(grepl("[M|m]", metric[1]) & grepl("[E|e]", metric[2])) fit = mxFDAobject@functional_cox$entropy[[model]]

  #model class
  mod_class = class(fit)

  #calculate
  range_s = range(mxfundata[[r]], na.rm = TRUE)
  range_x = range(mxfundata[[value]], na.rm = TRUE)
  tind_pred <- seq(range_s[1], range_s[2], length.out = grid_length)
  xind_pred <- seq(range_x[1], range_x[2], length.out = grid_length)

  grid_coef <- data.frame(t_int = rep(tind_pred, grid_length),
                          func = rep(xind_pred, each = grid_length), l_int=1)
  grid_coef <- dplyr::bind_cols(grid_coef, mxfundata[1, ] %>%
                       dplyr::select(!!analysis_vars)) #using just the first age? using the first doesnt seem right

  pred <- predict(fit, newdata = grid_coef, type = 'terms')
  if(inherits(fit,"afcm")){
    pred <- matrix(pred[,"ti(t_int,func):l_int"], ncol = sqrt(nrow(grid_coef)), nrow = sqrt(nrow(grid_coef)), byrow = FALSE)
  }else{
    pred <- matrix(pred[,"s(t_int):l_int * func"], ncol = sqrt(nrow(grid_coef)), nrow = sqrt(nrow(grid_coef)), byrow = FALSE)
  }

  rownames(pred) <- tind_pred
  colnames(pred) <- xind_pred
  out = list(Surface = reshape2::melt(pred) %>%
               dplyr::rename(r = Var1, sumfun = Var2))

  #se section
  #Defaults to FALSE, and returns surface.
  #If TRUE, returns coefficient function with standard errors for LFCM,
  #or coefficient surface for AFCM with non NA values that are statistically significant
  if(inherits(fit,"afcm")){
    pred <- predict(fit, newdata = grid_coef, type = 'terms', se.fit = TRUE)

    pred <- data.frame(beta1 = pred$fit[,"ti(t_int,func):l_int"],
                          beta1_se = pred$se.fit[,"ti(t_int,func):l_int"])%>%
      mutate(lower = beta1 - stats::qnorm(p/2, lower.tail = FALSE)*beta1_se,
             upper = beta1 + stats::qnorm(p/2, lower.tail = FALSE)*beta1_se,
             beta1_sig = ifelse(sign(lower) == sign(upper), beta1, NA))

    pred <- matrix(pred$beta1_sig, ncol = sqrt(nrow(grid_coef)), nrow = sqrt(nrow(grid_coef)), byrow = FALSE)

    rownames(pred) <- tind_pred
    colnames(pred) <- xind_pred
    out$Prediction = reshape2::melt(pred) %>%
      dplyr::rename(r = Var1, sumfun = Var2)

  }else{
    grid_coef <- data.frame(t_int = tind_pred,
                            func = 1,
                            l_int=1)
    grid_coef <- cbind(grid_coef, mxfundata[1, ] %>%
                         dplyr::select(!!analysis_vars))

    pred <- predict(fit, newdata = grid_coef, type = 'terms', se.fit = TRUE)

    pred_df <- data.frame(r = tind_pred,
                          beta1 = pred$fit[,"s(t_int):l_int * func"],
                          beta1_se = pred$se.fit[,"s(t_int):l_int * func"])
    out$Prediction = pred_df
  }

  if(mod_class[1] == "lfcm"){
    out_object = methods::new("lfcmSurface",
                     Surface = out$Surface,
                     Prediction = out$Prediction,
                     Metric = paste0(toupper(metric[2]), " value"),
                     `P-value` = p)
  } else if(mod_class[1] == "afcm"){

    out_object = methods::new("afcmSurface",
                     Surface = out$Surface,
                     Prediction = out$Prediction,
                     Metric = paste0(toupper(metric[2]), " value"),
                     `P-value` = p)
  } else {
    stop("unknown model class for surface extraction")
  }
  return(out_object)
}

methods::setClass("lfcmSurface",
                  representation(Surface = "data.frame",
                                 Prediction = "data.frame",
                                 Metric = "character",
                                 `P-value` = 'numeric'))
methods::setClass("afcmSurface",
                  representation(Surface = "data.frame",
                                 Prediction = "data.frame",
                                 Metric = "character",
                                 `P-value` = 'numeric'))


