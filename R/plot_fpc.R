

#' Create plot of mean +/- scaled eigenfunctions from FPCA
#'
#' Produces a ggplot with mean plus or minus two standard deviations of a selected FPC.
#'
#' @param obj fpca object to be plotted.
#' @param pc_choice FPC to be plotted.
#'
#' @details `r lifecycle::badge('superseded')`
#'
#' @return object of class ggplot
#'
#' @author Julia Wrobel \email{`r juliawrobel_email`}
#'
#' @import ggplot2
#'
#' @export
plot_fpc = function(obj, pc_choice){

  ## NULLify global values called in ggplot
  plus = minus = mu = index = NULL

  efunctions = matrix(obj$efunctions, ncol = obj$npc)
  sqrt.evalues = diag(sqrt(obj$evalues), obj$npc, obj$npc)
  scaled_efunctions = efunctions %*% sqrt.evalues
  scaled_efuncs = scaled_efunctions[,pc_choice]

  df = data.frame(id = 1,
                  plus = obj$mu + 2 * scaled_efuncs,
                  minus = obj$mu - 2 * scaled_efuncs,
                  mu  = obj$mu,
                  index = seq(obj$index_range[1], obj$index_range[2], length.out = length(obj$mu)))

  ggplot(df, aes(index, mu)) +
    geom_line(lwd = 1) +
    geom_point(aes(y = plus), color = "blue", size = 2, shape = "+") +
    geom_point(aes(y = minus), color = "darkred", size = 2, shape = "-") +
    ggtitle(bquote(psi[.(pc_choice)]~(t) ~ ","
                   ~.(100*round(obj$evalues[pc_choice]/sum(obj$evalues),3)) ~ "% Variance")) +
    labs(x = "r", y = paste("fpc", pc_choice)) +
    theme_minimal()
}
