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
    ylim(0, 15) +
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
