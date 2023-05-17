#' extract_surface
#'
#' Internal function called by \code{TITLE: regression function} that transforms long format functional data for use in a linear or additive functional Cox model.
#'
#' @importFrom reshape2 melt
#' @import dplyr
#'
#' @return A \code{dataframe} with predicted surface for AFCM and LFCM fits for use in plotting
#'
#' @examples
#' # simulate data
#' set.seed(1001)
#'
#'
#' @param mxfundata Dataframe of spatial summary functions from multiplex imaging data, in long format.
#' @param fit AFCM or LFCM model fit object
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param grid_length Length of grid on which to evaluate coefficient functions.
#' @param analysis_vars Other variables used in modeling FCM fit.
extract_surface = function(mxfundata,
                           fit,
                       r = "r",
                       value = "fundiff",
                       grid_length = 100,
                       analysis_vars){

  range_s = range(mxfundata[[r]], na.rm = TRUE)
  range_x = range(mxfundata[[value]], na.rm = TRUE)
  tind_pred <- seq(range_s[1], range_s[2], length.out = grid_length)
  xind_pred <- seq(range_x[1], range_x[2], length.out = grid_length)

  grid_coef <- data.frame(t_int = rep(tind_pred, grid_length),
                          func = rep(xind_pred, each = grid_length), l_int=1)
  grid_coef <- cbind(grid_coef, mxfundata[1, analysis_vars])

  pred <- predict(fit, newdata = grid_coef, type = 'terms')
  if(inherits(fit,"afcm")){
    pred <- matrix(pred[,"ti(t_int,func):l_int"], ncol = sqrt(nrow(grid_coef)), nrow = sqrt(nrow(grid_coef)), byrow = FALSE)
  }else{
    pred <- matrix(pred[,"s(t_int):l_int * func"], ncol = sqrt(nrow(grid_coef)), nrow = sqrt(nrow(grid_coef)), byrow = FALSE)
  }

  rownames(pred) <- tind_pred
  colnames(pred) <- xind_pred
  melt(pred) %>%
    rename(r = Var1, sumfun = Var2)
}



