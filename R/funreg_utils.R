#' process_fcm
#'
#' Internal function called by \code{TITLE: regression function} that transforms long format functional data for use in a linear or additive functional Cox model.
#'
#' @param mxfundata Dataframe of spatial summary functions from multiplex imaging data, in long format. Can be estimated using the function \code{extract_summary_functions} or provided separately.
#' @param id Character string, the name of the variable that identifies each unique subject.
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param knots Number of knots for defining spline basis.
#' @param analysis_vars Optional list of variables to be retained for downstream analysis.
#' @param quantile_transform If TRUE, a quantile transformation is applied to the functional predictor before modeling
#'
#' @return A \code{dataframe} with matrix-valued covariates \code{l_int}, \code{t_int}, and \code{func} for use in a linear or additive functional Cox model.
#'
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
#'
#' @importFrom tidyr pivot_wider
#' @import dplyr
#'
#' @examples
#' # simulate data
#' set.seed(1001)
#'
process_fcm <- function(mxfundata,
                       id,
                       r = "r",
                       value = "fundiff",
                       analysis_vars,
                       quantile_transform){

  tind <- sort(unique(mxfundata[[r]]))

  mxfundata <- mxfundata %>%
    select(all_of(c(id, r, value, analysis_vars))) %>%
    pivot_wider(names_from = r,
                names_prefix = "r_",
                values_from = value)

  func <- mxfundata %>% select(starts_with("r_")) %>% as.matrix()
  mxfundata <- mxfundata %>% select(-starts_with("r_"))


  # set up data for use with the  AFCM
  # number of objects per subject
  I <- length(unique(mxfundata[[id]]))
  J <- length(tind)

  ### lmat: numeric integration of the functional term
  grid_spacing <- diff(tind)[floor(J/2)] # spacing between time points for Riemannian integration, assumes equal grid
  mxfundata$l_int <- matrix(grid_spacing, ncol= J, nrow= I)

  ### tmat: time indices
  mxfundata$t_int <- matrix(tind, ncol=J, nrow = I, byrow=TRUE)


  if(quantile_transform){
    # quantile transformed version of the data
    mxfundata$func <- apply(func, 2, function(y) ecdf(y)(y))
  }else{
    mxfundata$func <- func
  }

  mxfundata
}



#' impute_fpca
#'
#' Internal function called by \code{TITLE: regression function} that imputes missing data in functional predictors using FPCA.
#'
#'
#' @param mxfundata Dataframe of spatial summary functions from multiplex imaging data, in long format. Can be estimated using the function \code{extract_summary_functions} or provided separately.
#' @param id Character string, the name of the variable that identifies each unique subject.
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param knots Number of knots for defining spline basis.
#' @param analysis_vars Optional list of variables to be retained for downstream analysis.
#' @param smooth Option to smooth data using FPCA.
#'
#' @return A \code{dataframe} where the missing function values (NA) for the \code{value} variable have been replaced with estimates from FPCA.
#'
#' @author Julia Wrobel \email{julia.wrobel@@emory.edu}
#'
#' @importFrom refund fpca.face
#' @importFrom tidyr pivot_wider pivot_longer
#' @import dplyr
#'
#' @examples
#' # simulate data
#' set.seed(1001)
#'
impute_fpca = function(mxfundata,
                       id,
                       r = "r",
                       value = "fundiff",
                       knots = NULL,
                       analysis_vars,
                       smooth){

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
  mx_fpc <- fpca.face(Y = mat, knots = knots)

  # impute missing values using FPCA fitted values
  mxfundata = mxfundata %>%
    # convert to long format
    pivot_longer(contains("r_"), names_to = r, values_to = value,
                 names_prefix = "r_", names_transform = as.numeric) %>%
    mutate(Yhat = as.vector(t(mx_fpc$Yhat)),
           imputed = ifelse(is.na(get(value)) | is.nan(get(value)), Yhat, get(value)))

  if(smooth){
    mxfundata[[value]] <- mxfundata$Yhat
  }else{
    mxfundata[[value]] <- mxfundata$imputed
  }



  select(mxfundata, -imputed, -Yhat)
}

#' extract_c
#'
#' Function to calculate c-index from a an AFCM or LFCM fit
#'
#' @param fit fit AFCM or LFCM model fit object.
#' @param survival_time Vector of survival/censoring times
#' @param event Survival statust (0 = censored, 1 = event)
#'
#' @return c-index
#'
#' @author Erjia Cui
#'
#' @export
extract_c <- function(fit, survival_time, event){
  eta <- predict(fit, type = "link")
  etimes <- sort(unique(survival_time[event==1]))
  num <- denom <- 0
  for(e in seq_along(etimes)){
    ## current time
    ti  <- etimes[e]
    ## subjects who experienced an event at current time
    event_t <- which(survival_time == ti & event==1)
    ## subjects with observed times beyond event current time
    control_t <- which(survival_time > ti)

    for(i in seq_along(event_t)){
      num   <- num + sum( (eta[control_t] > eta[event_t[i]] ) ) + 0.5*sum(eta[control_t] == eta[event_t[i]])
    }
    denom <- denom + length(event_t)*length(control_t)
  }
  1-num/denom
}

