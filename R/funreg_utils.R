#' process_fcm
#'
#' Internal function called by \code{run_fcm} that transforms long format functional data for use in a linear or additive functional Cox model.
#'
#' @param mxfundata Dataframe of spatial summary functions from multiplex imaging data, in long format. Can be estimated using the function \code{extract_summary_functions} or provided separately.
#' @param id Character string, the name of the variable that identifies each unique subject.
#' @param r Character string, the name of the variable that identifies the function domain (usually a radius for spatial summary functions). Default is "r".
#' @param value Character string, the name of the variable that identifies the spatial summary function values. Default is "fundiff".
#' @param analysis_vars Optional list of variables to be retained for downstream analysis.
#' @param quantile_transform If TRUE, a quantile transformation is applied to the functional predictor before modeling
#'
#' @return A \code{dataframe} with matrix-valued covariates \code{l_int}, \code{t_int}, and \code{func} for use in a linear or additive functional Cox model.
#'
#' @keywords internal
#'
#' @author Julia Wrobel \email{`r juliawrobel_email`}
#' @author Alex Soupir \email{`r alexsoupir_email`}
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
                       quantile_transform = FALSE){

  tind <- sort(unique(mxfundata[[r]]))

  mxfundata <- mxfundata %>%
    dplyr::select(dplyr::all_of(c(id, r, value, analysis_vars))) %>%
    tidyr::pivot_wider(names_from =  dplyr::all_of(r),
                names_prefix = "r_",
                values_from =  dplyr::all_of(value))

  func <- mxfundata %>% dplyr::select(dplyr::all_of(dplyr::starts_with("r_"))) %>% as.matrix()
  mxfundata <- mxfundata %>% dplyr::select(!dplyr::starts_with("r_"))


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
#' @keywords internal
#'
#' @author Julia Wrobel \email{`r juliawrobel_email`}
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
    dplyr::select(dplyr::all_of(c(id, r, value, analysis_vars))) %>%
    tidyr::pivot_wider(names_from = r,
                names_glue = "r_{round(r, 2)}",
                values_from = value)

  mat <- mxfundata %>%
    dplyr::select(dplyr::all_of(starts_with("r_"))) %>%
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

  out = mxfundata #%>%
    #dplyr::select(dplyr::all_of(c(analysis_vars, "imputed", "Yhat")))
  return(out)
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
#' @keywords internal
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



#' pfr2
#'
#' hot-fix for the prf function from refund::
#'
#' Assumes users wont use gamm4 or lme4 packages for dependency purposes
#'
#' @details see `?refund::pfv`
#'
#' @keywords internal
#'
#' @importFrom mgcv gam gam.fit bam s te t2
#'
#' @author Jonathan Geller
#' @author Edited by Alex Soupir \email{`r alexsoupir_email`}
#'
#' @export
pfr2 =  function (formula = NULL, fitter = NA, method = "REML", ...){
  if (!is(formula, "formula")) {
    warning(paste0("The interface for pfr() has changed to using a formula ",
                   "argument, with linear functional terms specified by lf(). ",
                   "See ?pfr for details. The old interface will be deprecated ",
                   "in the next refund release."))
    call <- sys.call()
    call[[1]] <- as.symbol("pfr_old")
    fit <- eval(call)
    return(fit)
  }
  call <- match.call()
  dots <- list(...)
  if (length(dots)) {
    validDots <- if (!is.na(fitter) && fitter == "gamm4") {
      c(names(formals(gamm4)), names(formals(lmer)))
    }
    else {
      c(names(formals(gam)), names(formals(gam.fit)))
    }
    notUsed <- names(dots)[!(names(dots) %in% validDots)]
    if (length(notUsed))
      warning("Arguments <", paste(notUsed, collapse = ", "),
              "> supplied but not used.")
  }
  tf <- terms.formula(formula, specials = c("s", "te", "t2",
                                            "lf", "af", "lf.vd", "re", "peer", "fpc"))
  trmstrings <- attr(tf, "term.labels")
  terms <- sapply(trmstrings, function(trm) as.call(parse(text = trm))[[1]],
                  simplify = FALSE)
  frmlenv <- environment(formula)
  specials <- attr(tf, "specials")
  where.af <- specials$af - 1
  where.lf <- specials$lf - 1
  where.pr <- specials$peer - 1
  where.fp <- specials$fpc - 1
  where.s <- specials$s - 1
  where.te <- specials$te - 1
  where.t2 <- specials$t2 - 1
  where.re <- specials$re - 1
  where.lf.vd <- specials$lf.vd - 1
  where.all <- c(where.af, where.lf, where.s, where.te, where.t2,
                 where.re, where.lf.vd, where.pr, where.fp)
  if (length(trmstrings)) {
    where.par <- which(!(1:length(trmstrings) %in% where.all))
  }
  else where.par <- numeric(0)
  responsename <- attr(tf, "variables")[2][[1]]
  newfrml <- paste(responsename, "~", sep = "")
  newfrmlenv <- new.env()
  evalenv <- if ("data" %in% names(call))
    eval.parent(call$data)
  else NULL
  nobs <- length(eval(responsename, envir = evalenv, enclos = frmlenv))
  if (missing(fitter) || is.na(fitter)) {
    fitter <- ifelse(nobs > 1e+05, "bam", "gam")
  }
  fitter <- as.symbol(fitter)
  if (as.character(fitter) == "bam" && !("chunk.size" %in%
                                         names(call))) {
    call$chunk.size <- max(nobs/5, 10000)
  }
  if (as.character(fitter) == "gamm4")
    stopifnot(length(where.te) < 1)
  assign(x = deparse(responsename), value = as.vector(t(eval(responsename,
                                                             envir = evalenv, enclos = frmlenv))), envir = newfrmlenv)
  newtrmstrings <- attr(tf, "term.labels")
  if (!attr(tf, "intercept")) {
    newfrml <- paste(newfrml, "0 +", sep = "")
  }
  where.refund <- c(where.af, where.lf, where.lf.vd, where.pr,
                    where.fp, where.re)
  if (length(where.refund)) {
    fterms <- lapply(terms[where.refund], function(x) {
      eval(x, envir = evalenv, enclos = frmlenv)
    })
    newtrmstrings[where.refund] <- sapply(fterms, function(x) {
      refund:::safeDeparse(x$call)
    })
    lapply(fterms, function(x) {
      lapply(names(x$data), function(nm) {
        assign(x = nm, value = x$data[[nm]], envir = newfrmlenv)
        invisible(NULL)
      })
      if ("xt" %in% names(x$call)) {
        xtvars <- all.vars(x$call$xt)
        if (length(xtvars)) {
          sapply(xtvars, function(xtvar) {
            xtvarval <- eval(as.name(xtvar), envir = evalenv,
                             enclos = frmlenv)
            assign(x = xtvar, value = xtvarval, envir = parent.env(newfrmlenv))
            invisible(NULL)
          })
        }
      }
      invisible(NULL)
    })
    fterms <- lapply(fterms, function(x) x[names(x) != "data"])
  }
  else fterms <- NULL
  where.mgcv <- c(where.par, where.s, where.te, where.t2)
  if (length(where.mgcv)) {
    if ("data" %in% names(call))
      frmlenv <- list2env(eval(call$data, envir = frmlenv), frmlenv)
    lapply(terms[where.mgcv], function(x) {
      nms <- if (!is.null(names(x))) {
        all.vars(x[names(x) == ""])
      }
      else all.vars(x)
      sapply(nms, function(nm) {
        stopifnot(length(get(nm, envir = frmlenv)) ==
                    nobs)
        assign(x = nm, value = get(nm, envir = frmlenv),
               envir = newfrmlenv)
        invisible(NULL)
      })
      invisible(NULL)
    })
  }
  newfrml <- formula(paste(newfrml, paste(newtrmstrings, collapse = "+")))
  environment(newfrml) <- newfrmlenv
  pfrdata <- refund:::list2df(as.list(newfrmlenv))
  datameans <- sapply(as.list(newfrmlenv), function(x) {
    if (is.numeric(x) | is.logical(x)) {
      mean(x)
    }
    else if (is.factor(x)) {
      names(which.max(table(x)))
    }
    else NA
  }, simplify = FALSE)
  newcall <- refund:::expand.call(pfr, call)
  newcall$fitter <- newcall$bs.int <- newcall$bs.yindex <- NULL
  newcall$formula <- newfrml
  newcall$data <- quote(pfrdata)
  newcall$method <- method
  newcall[[1]] <- fitter
  res <- eval(newcall)
  res.smooth <- if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$smooth
  }
  else res$smooth
  names(res.smooth) <- sapply(res.smooth, function(x) x$label)
  if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$smooth <- res.smooth
  }
  else {
    res$smooth <- res.smooth
  }
  termtype <- rep("par", length(terms))
  for (i in 1:length(specials)) termtype[specials[[i]] - 1] <- names(specials)[i]
  ret <- list(formula = formula, responsename = responsename,
              nobs = nobs, termnames = names(terms), termtype = termtype,
              datameans = datameans, ft = fterms)
  if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$pfr <- ret
    class(res$gam) <- c("pfr", class(res$gam))
  }
  else {
    res$pfr <- ret
    class(res) <- c("pfr", class(res))
  }
  return(res)
}

