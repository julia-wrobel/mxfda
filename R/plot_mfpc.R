#' Create plot of mean +/- scaled eigenfunctions from FPCA
#'
#' Produces a ggplot with mean plus or minus two standard deviations of a selected FPC.
#'
#' @param obj fpca object to be plotted.
#' @param pc_choice_level1,pc_choice_level2 FPC to be plotted.
#'
#' @details `r lifecycle::badge('superseded')`
#'
#' @return list of objects of class ggplot
#'
#' @author Julia Wrobel \email{julia.wrobel@@emory.edu}
#'
#' @import ggplot2
#'
#' @export
#'
plot_mfpc = function(obj, pc_choice_level1, pc_choice_level2){

  efunctions = obj$efunctions
  sqrt.evalues = lapply(1:2, function(i) diag(sqrt(obj$evalues[[i]]), nrow = obj$npc[[i]]))
  scaled_efunctions = lapply(1:2, function(i) efunctions[[i]] %*% sqrt.evalues[[i]])
  names(scaled_efunctions) <- names(sqrt.evalues) <- names(efunctions) <- c("level1", "level2")

  PCchoice = list(as.numeric(pc_choice_level1), as.numeric(pc_choice_level2))
  names(PCchoice) <- c("level1", "level2")
  scaled_efuncs = lapply(1:2, function(i) scaled_efunctions[[i]][,PCchoice[[i]]])
  mu = data.frame(grid = 1:length(obj$mu), value = obj$mu)

  p <- lapply(1:2, function(i){
    ggplot(mu, aes(x = grid, y = value)) +
      geom_line(lwd=1) +
      geom_point(data = data.frame(grid =1:length(obj$mu),value =  obj$mu + 2*scaled_efuncs[[i]]), color = "blue", size = 2, shape = '+') +
      geom_point(data = data.frame(grid =1:length(obj$mu), value = obj$mu - 2*scaled_efuncs[[i]]), color = "darkred", size = 2, shape = "-") +
      ggtitle(paste("level ", i, "fpc", PCchoice[[i]])) +
      labs(x = "r", y = paste("level ", i, "fpc", PCchoice[[i]])) +
      theme_minimal()
  })

  # return list of plots (one for each level)
  p

# end this
}
