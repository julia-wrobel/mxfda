#' extract_surface
#'
#' Internal function called by \code{TITLE: regression function} that transforms long format functional data for use in a linear or additive functional Cox model.
#'
#'
#' @param mxFDAobject object of class `mxFDA` with model `model` calculated wihtin
#' @param metric spatial summary function to extract surface for
#' @param model character string for the name of the model for `metric` data
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param grid_length Length of grid on which to evaluate coefficient functions.
#' @param analysis_vars Other variables used in modeling FCM fit.
#' @param se Defaults to FALSE, and returns surface. If TRUE, returns coefficient function with standard errors for LFCM, or coefficient surface for AFCM with non NA values that are statistically significant
#' @param filter_cols a named vector of factors to filter summary functions to in `c(Derived_Column = "Level_to_Filter")` format
#'
#' @return A \code{dataframe} with predicted surface for AFCM and LFCM fits for use in plotting
#'
#' @author Julia Wrobel \email{julia.wrobel@@emory.edu}
#' @author Alex Soupir \email{alex.soupir@@moffitt.org}
#'
#' @importFrom reshape2 melt
#' @importFrom mgcv predict.gam
#' @import dplyr
#'
#' @examples
#' # simulate data
#' set.seed(1001)
#'
#' @export
extract_surface = function(mxFDAobject,
                           metric,
                           model = NULL,
                           r = "r",
                           value = "fundiff",
                           grid_length = 100,
                           analysis_vars,
                           se = FALSE,
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
  if(se){
    if(inherits(fit,"afcm")){
      pred <- predict(fit, newdata = grid_coef, type = 'terms', se.fit = TRUE)

      pred <- data.frame(beta1 = pred$fit[,"ti(t_int,func):l_int"],
                            beta1_se = pred$se.fit[,"ti(t_int,func):l_int"])%>%
        mutate(lower = beta1 - 1.96*beta1_se,
               upper = beta1 + 1.96*beta1_se,
               beta1_sig = ifelse(sign(lower) == sign(upper), beta1, NA))

      pred <- matrix(pred$beta1_sig, ncol = sqrt(nrow(grid_coef)), nrow = sqrt(nrow(grid_coef)), byrow = FALSE)

    }else{
      grid_coef <- data.frame(t_int = tind_pred,
                              func = 1,
                              l_int=1)
      grid_coef <- cbind(grid_coef, mxfundata[1, analysis_vars])

      pred <- predict(fit, newdata = grid_coef, type = 'terms', se.fit = TRUE)

      pred_df <- data.frame(r = tind_pred,
                            beta1 = pred$fit[,"s(t_int):l_int * func"],
                            beta1_se = pred$se.fit[,"s(t_int):l_int * func"])
      return(pred_df)

    }
  }

  rownames(pred) <- tind_pred
  colnames(pred) <- xind_pred
  melt(pred) %>%
    rename(r = Var1, sumfun = Var2)
}



