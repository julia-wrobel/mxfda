#' Entropy
#'
#' @param X object of class `ppp` from spatstat with 2 marks
#' @param i ignored
#' @param j ignored
#' @param r_vec vector of length wanted for breaks (will be rescaled) with max value at max for measuring entropy
#' @param correction ignored
#'
#' @details `r lifecycle::badge('experimental')`
#'
#' @return data frame with entropy calculated for `length(r_vec)` bins within 0 to `max(r_vec)`
#'
#' @author Thao Vu \email{`r thaovu_email`}
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @references Vu, T., Seal, S., Ghosh, T., Ahmadian, M., Wrobel, J., & Ghosh, D. (2023).
#' FunSpace: A functional and spatial analytic approach to cell imaging data using entropy measures.
#' \emph{PLOS Computational Biology}, 19(9), e1011490.
#' @references Altieri, L., Cocchi, D., & Roli, G. (2018).
#' A new approach to spatial entropy measures.
#' \emph{Environmental and ecological statistics}, 25, 95-110.
#'
#' @export
entropy = function(X,
                   i,
                   j,
                   r_vec,
                   correction
                   #n_break = 50,
                   #rmax = 400
                             ){
  breaks = length(r_vec)
  # Compute spatial entropy
  spa_entropy = SimDesign::quiet(
    SpatEntropy::altieri(X, distbreak = r_vec[-breaks],
                         verbose = FALSE, plotout = FALSE)
  )

  if (length(spa_entropy$SPI.terms) < breaks){
    breaks = length(spa_entropy$SPI.terms)
    r = unique(c(spa_entropy$distance.breaks[,1],spa_entropy$distance.breaks[,2]))[-c(1)]
    r[breaks] = max(r_vec)
    r_vec = r
    }

  df = data.frame(r = r_vec, spatial_entropy = spa_entropy$SPI.terms,
                  residual_entropy  =  spa_entropy$RES.terms)

  df

}






