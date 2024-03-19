#' Plot mxFDA object
#'
#' @param x object of class `mxFDA` to be plotted
#' @param ... additional paramters including `y`, `what`, and `sampleID` to inform whats to be plotted
#' @param filter_cols column key to filter
#'
#' @details `r lifecycle::badge('stable')`
#'
#' If there are multiple metrics that are included in the derived table, an extra parameter `filter_cols`
#' in the format of `c(Derived_Column = "Level_to_Filter")` will return curves from the `Derived_Column`
#' with the level `Level_to_Filter`
#'
#' When plotting mFPCA objects, additional arguments `level1` and `level2` help indicate which FPCA from level 1 and level 2 to plot
#'
#' @return object of class `ggplot` compatible the `ggplot2` aesthetics
#' @export
#'
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @examples
#' #set seed
#' set.seed(333)
#' #plotting summary
#' data("ovarian_FDA")
#' plot(ovarian_FDA, y = 'fundiff', what = 'uni g')
#' #running fpca
#' ovarian_FDA = run_fpca(ovarian_FDA, metric = "uni g", r = "r", value = "fundiff",
#'                        lightweight = TRUE,
#'                        pve = .99)
#' #plot fpca
#' plot(ovarian_FDA, what = 'uni g fpca', pc_choice = 1)
#'
plot.mxFDA = function(x, filter_cols = NULL, ...){
  params = as.list(substitute(list(...)))
  if(is.null(params$what)) stop("need to provide what to plot")
  params$what = unlist(strsplit(params$what, split = " "))
  #return(params)

  if(length(params$what) == 2){
    if(is.null(params$y)) stop("for summary functions, need to provide a column name")
    dat = get_data(x, params$what, type = 'summaries') %>%
      filter_data(filter_cols)

    pl = ggplot2::ggplot(data = dat,
                    ggplot2::aes(x = r, y = get(params$y))) +
      ggplot2::geom_line(alpha = 0.2, ggplot2::aes(group = get(x@sample_key))) +
      ggplot2::labs(title = params$what[1], y = params$what[2])

    return(pl)
  }
  #non-mixed fpca
  if(grepl("^fpca", params$what[3], ignore.case = TRUE)){
    dat = get_data(x, params$what, type = "fpca")
    if(is.null(params$pc_choice)) stop("need to provide a numeric principal component to plot")

    plus = minus = mu = index = NULL
    obj = dat$fpc_object
    efunctions = matrix(obj$efunctions, ncol = obj$npc)
    sqrt.evalues = diag(sqrt(obj$evalues), obj$npc, obj$npc)
    scaled_efunctions = efunctions %*% sqrt.evalues
    scaled_efuncs = scaled_efunctions[,params$pc_choice]

    df = data.frame(id = 1,
                    plus = obj$mu + 2 * scaled_efuncs,
                    minus = obj$mu - 2 * scaled_efuncs,
                    mu  = obj$mu,
                    index = seq(obj$index_range[1], obj$index_range[2], length.out = length(obj$mu)))

    pl = ggplot2::ggplot(df, ggplot2::aes(index, mu)) +
      ggplot2::geom_line(lwd = 1) +
      ggplot2::geom_point(ggplot2::aes(y = plus), color = "blue", size = 2, shape = "+") +
      ggplot2::geom_point(ggplot2::aes(y = minus), color = "darkred", size = 2, shape = "-") +
      ggplot2::ggtitle(bquote(psi[.(params$pc_choice)]~(t) ~ ","
                     ~.(100*round(obj$evalues[params$pc_choice]/sum(obj$evalues),3)) ~ "% Variance")) +
      ggplot2::labs(x = "r", y = paste("fpc", params$pc_choice)) +
      ggplot2::theme_minimal()

    return(pl)
  }
  #mixed fpca (mfpca)
  if(grepl("^mfpca", params$what[3], ignore.case = TRUE)){
    dat = get_data(x, params$what, type = "mfpca")

    #extract object
    obj = dat$mfpc_object
    efunctions = obj$efunctions
    sqrt.evalues = lapply(1:2, function(i) diag(sqrt(obj$evalues[[i]]), nrow = obj$npc[[i]]))
    scaled_efunctions = lapply(1:2, function(i) efunctions[[i]] %*% sqrt.evalues[[i]])
    names(scaled_efunctions) <- names(sqrt.evalues) <- names(efunctions) <- c("level1", "level2")

    PCchoice = list(as.numeric(params$level1), as.numeric(params$level2))
    names(PCchoice) <- c("level1", "level2")
    scaled_efuncs = lapply(1:2, function(i) scaled_efunctions[[i]][,PCchoice[[i]]])
    mu = data.frame(grid = 1:length(obj$mu), value = obj$mu)

    pl <- lapply(1:2, function(i){
      ggplot(mu, aes(x = grid, y = value)) +
        geom_line(lwd=1) +
        geom_point(data = data.frame(grid =1:length(obj$mu),value =  obj$mu + 2*scaled_efuncs[[i]]), color = "blue", size = 2, shape = '+') +
        geom_point(data = data.frame(grid =1:length(obj$mu), value = obj$mu - 2*scaled_efuncs[[i]]), color = "darkred", size = 2, shape = "-") +
        ggtitle(paste("level ", i, "fpc", PCchoice[[i]])) +
        labs(x = "r", y = paste("level ", i, "fpc", PCchoice[[i]])) +
        theme_minimal()
    })

    # return list of plots (one for each level)
    return(pl)
  }

}

#' Plot lfcm surface
#'
#' @param x object of class `lfcmSurface` to be plotted
#' @param ... currently ignored
#'
#' @return object compatable with ggplot2
#' @export
#'
#' @author Julia Wrobel \email{`r juliawrobel_email`}
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
plot.lfcmSurface = function(x, ...){
  x@Prediction %>%
    ggplot(aes(r, exp(beta1))) +
    geom_line()+
    geom_line(aes(x=r, y = exp(beta1-2*beta1_se)), linetype = 'longdash')+
    geom_line(aes(x=r, y = exp(beta1+2*beta1_se)), linetype = 'longdash')+
    geom_hline(yintercept = 1, color = "red", linetype = 3) +
    labs(y = expression(e^hat(beta)(r)), title = "LFCM Hazard Ratio") +
    theme(axis.title = element_text(size = 15))
}

#' Plot afcm object
#'
#' @param x object of class `afcmSurface` to be plotted
#' @param ... currently ignored
#'
#' @return object compatable with ggplot2
#' @export
#'
#' @author Julia Wrobel \email{`r juliawrobel_email`}
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
plot.afcmSurface = function(x, ...){
  bind_rows(
    x@Surface %>% mutate(type = "AFCM surface"),
    x@Prediction %>% mutate(type = paste0("p < ", x@`P-value`))
  ) %>%
    mutate(hr = exp(value)) %>%
    ggplot(aes(r, sumfun, fill = hr)) +
    geom_tile() +
    labs(y = x@Metric) +
    scale_fill_distiller(name = expression(e^hat(F)(.,.)), na.value = "transparent",
                         type = "div") +
    facet_wrap(~type, ncol = 2)
}


#' Plot sofr object
#'
#' @param x object of class `sofr` to be plotted
#' @param ... currently ignored
#'
#' @return object compatable with ggplot2
#' @export
#'
#' @author Julia Wrobel \email{`r juliawrobel_email`}
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
plot.sofr = function(x, ...){
  class(x) = class(x)[-1]
  class(x) = class(x)[-1]
  #class(x) <- append("get_dat", class(x))
  v = get_plot_dat(x, ...)
  #dev.off() #because apparently you cant just turn off the default printing for some reason....terrible design for R
  df = data.frame(x = v$x, fit = v$fit) %>%
    dplyr::mutate(`fit+se` = fit + v$se,
                  `fit-se` = fit - v$se)
  pl = df %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(aes(x = x, y = fit)) +
    ggplot2::geom_line(aes(x = x, y = `fit+se`), linetype = "dashed") +
    ggplot2::geom_line(aes(x = x, y = `fit-se`), linetype = 'dashed') +
    ggplot2::labs(x = v$xlab, y = v$ylab)
  return(pl)
}



#' @keywords internal
#'
#' @author Refunders
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
get_plot_dat = function (x, residuals = FALSE, rug = NULL, se = TRUE, pages = 0,
                     select = NULL, scale = -1, n = 100, n2 = 40, n3 = 3, theta = 30,
                     phi = 30, jit = FALSE, xlab = NULL, ylab = NULL, main = NULL,
                     ylim = NULL, xlim = NULL, too.far = 0.1, all.terms = FALSE,
                     shade = FALSE, shade.col = "gray80", shift = 0, trans = I,
                     seWithMean = FALSE, unconditional = FALSE, by.resids = FALSE,
                     scheme = 0, ...)
{
  sub.edf <- function(lab, edf) {
    pos <- regexpr(":", lab)[1]
    if (pos < 0) {
      pos <- nchar(lab) - 1
      lab <- paste(substr(lab, start = 1, stop = pos),
                   ",", round(edf, digits = 2), ")", sep = "")
    }
    else {
      lab1 <- substr(lab, start = 1, stop = pos - 2)
      lab2 <- substr(lab, start = pos - 1, stop = nchar(lab))
      lab <- paste(lab1, ",", round(edf, digits = 2),
                   lab2, sep = "")
    }
    lab
  }
  if (is.null(rug))
    rug <- if (nrow(x$model) > 10000){
      FALSE
    } else TRUE
  if (unconditional) {
    if (is.null(x$Vc)){
      warning("Smoothness uncertainty corrected covariance not available")
    } else x$Vp <- x$Vc
  }
  w.resid <- NULL
  if (length(residuals) > 1) {
    if (length(residuals) == length(x$residuals)){
      w.resid <- residuals
    } else warning("residuals argument to plot.gam is wrong length: ignored")
    partial.resids <- TRUE
  } else partial.resids <- residuals
  m <- length(x$smooth)
  if (length(scheme) == 1)
    scheme <- rep(scheme, m)
  if (length(scheme) != m) {
    warn <- paste("scheme should be a single number, or a vector with",
                  m, "elements")
    warning(warn)
    scheme <- rep(scheme[1], m)
  }
  order <- if (is.list(x$pterms)){
    unlist(lapply(x$pterms, attr, "order"))
  } else attr(x$pterms, "order")
  if (all.terms){
    n.para <- sum(order == 1)
  } else n.para <- 0
  if (se) {
    if (is.numeric(se))
      se2.mult <- se1.mult <- se
    else {
      se1.mult <- 2
      se2.mult <- 1
    }
    if (se1.mult < 0)
      se1.mult <- 0
    if (se2.mult < 0)
      se2.mult <- 0
  } else se1.mult <- se2.mult <- 1
  if (se && x$Vp[1, 1] < 0) {
    se <- FALSE
    warning("No variance estimates available")
  }
  if (partial.resids) {
    if (is.null(w.resid)) {
      if (is.null(x$residuals) || is.null(x$weights))
        partial.resids <- FALSE
      else {
        wr <- sqrt(x$weights)
        w.resid <- x$residuals * wr
      }
    }
    if (partial.resids)
      fv.terms <- predict(x, type = "terms")
  }
  pd <- list()
  i <- 1
  if (m > 0)
    for (i in 1:m) {
      first <- x$smooth[[i]]$first.para
      last <- x$smooth[[i]]$last.para
      edf <- sum(x$edf[first:last])
      term.lab <- sub.edf(x$smooth[[i]]$label, edf)
      attr(x$smooth[[i]], "coefficients") <- x$coefficients[first:last]
      P <- refund:::plot.mgcv.smooth(x$smooth[[i]], P = NULL, data = x$model,
                partial.resids = partial.resids, rug = rug,
                se = se, scale = scale, n = n, n2 = n2, n3 = n3,
                theta = theta, phi = phi, jit = jit, xlab = xlab,
                ylab = ylab, main = main, label = term.lab,
                ylim = ylim, xlim = xlim, too.far = too.far,
                shade = shade, shade.col = shade.col, se1.mult = se1.mult,
                se2.mult = se2.mult, shift = shift, trans = trans,
                by.resids = by.resids, scheme = scheme[i])
      if (is.null(P))
        pd[[i]] <- list(plot.me = FALSE)
      else if (is.null(P$fit)) {
        p <- x$coefficients[first:last]
        offset <- attr(P$X, "offset")
        if (is.null(offset)){
          P$fit <- P$X %*% p
        } else P$fit <- P$X %*% p + offset
        if (!is.null(P$exclude))
          P$fit[P$exclude] <- NA
        if (se && P$se) {
          if (seWithMean && attr(x$smooth[[i]], "nCons") > 0) {
            if (length(x$cmX) < ncol(x$Vp))
              x$cmX <- c(x$cmX, rep(0, ncol(x$Vp) -
                                      length(x$cmX)))
            if (seWithMean == 2)
              x$cmX[-(1:x$nsdf)] <- 0
            X1 <- matrix(x$cmX, nrow(P$X), ncol(x$Vp),
                         byrow = TRUE)
            meanL1 <- x$smooth[[i]]$meanL1
            if (!is.null(meanL1))
              X1 <- X1/meanL1
            X1[, first:last] <- P$X
            lpi <- attr(x$formula, "lpi")
            if (is.null(lpi))
              se.fit <- sqrt(pmax(0, rowSums(methods::as(X1 %*%
                                                  x$Vp, "matrix") * X1)))
            else {
              ii <- rep(0, 0)
              for (q in 1:length(lpi)) if (any(first:last %in%
                                               lpi[[q]]))
                ii <- c(ii, lpi[[q]])
              se.fit <- sqrt(pmax(0, rowSums(methods::as(X1[,
                                                   ii] %*% x$Vp[ii, ii], "matrix") * X1[,
                                                                                        ii])))
            }
          } else se.fit <- sqrt(pmax(0, rowSums(methods::as(P$X %*%
                                                   x$Vp[first:last, first:last, drop = FALSE],
                                                 "matrix") * P$X)))
          if (!is.null(P$exclude))
            se.fit[P$exclude] <- NA
        }
        if (partial.resids) {
          P$p.resid <- fv.terms[, length(order) + i] +
            w.resid
        }
        if (se && P$se)
          P$se <- se.fit * P$se.mult
        P$X <- NULL
        P$plot.me <- TRUE
        pd[[i]] <- P
        rm(P)
      }
      else {
        if (partial.resids) {
          P$p.resid <- fv.terms[, length(order) + i] +
            w.resid
        }
        P$plot.me <- TRUE
        pd[[i]] <- P
        rm(P)
      }
    }

  return(pd[[1]])
}

