#' Plot mxFDA object
#'
#' @param x object of class `mxFDA` to be plotted
#' @param ... additional paramters including `y`, `what`, and `sampleID` to inform whats to be plotted
#'
#' @return object compatable with ggplot2
#' @export
#'
plot.mxFDA = function(x,...){
  params = as.list(substitute(list(...)))
  params$what = unlist(strsplit(params$what, split = " "))

  if(length(params$what) == 2){
    dat = get_data(x, params$what, type = 'summaries')
    # if(grepl("[B|b]", params$what[1]) & grepl("[K|k]", params$what[2])) dat = x@`Bivariate Summaries`$Kcross
    # if(grepl("[B|b]", params$what[1]) & grepl("[G|g]", params$what[2])) dat = x@`Bivariate Summaries`$Gcross
    # if(grepl("[B|b]", params$what[1]) & grepl("[L|l]", params$what[2])) dat = x@`Bivariate Summaries`$Lcross
    # if(grepl("[U|u]", params$what[1]) & grepl("[K|k]", params$what[2])) dat = x@`Univariate Summaries`$Kest
    # if(grepl("[U|u]", params$what[1]) & grepl("[G|g]", params$what[2])) dat = x@`Univariate Summaries`$Gest
    # if(grepl("[U|u]", params$what[1]) & grepl("[L|l]", params$what[2])) dat = x@`Univariate Summaries`$Lest

    pl = ggplot2::ggplot(data = dat,
                    ggplot2::aes(x = r, y = get(params$y))) +
      ggplot2::geom_line(alpha = 0.2, ggplot2::aes(group = get(params$sampleID))) +
      ggplot2::labs(title = params$what[1], y = params$what[2])

    return(pl)
  }

  if(grepl("fpca", params$what[3], ignore.case = TRUE)){
    dat = get_data(x, params$what, type = "fpca")
    # if(grepl("[B|b]", params$what[1]) & grepl("[K|k]", params$what[2])) dat = x@`Functional Data`$Kcross
    # if(grepl("[B|b]", params$what[1]) & grepl("[G|g]", params$what[2])) dat = x@`Functional Data`$Gcross
    # if(grepl("[B|b]", params$what[1]) & grepl("[L|l]", params$what[2])) dat = x@`Functional Data`$Lcross
    # if(grepl("[U|u]", params$what[1]) & grepl("[K|k]", params$what[2])) dat = x@`Functional Data`$Kest
    # if(grepl("[U|u]", params$what[1]) & grepl("[G|g]", params$what[2])) dat = x@`Functional Data`$Gest
    # if(grepl("[U|u]", params$what[1]) & grepl("[L|l]", params$what[2])) dat = x@`Functional Data`$Lest

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

}

