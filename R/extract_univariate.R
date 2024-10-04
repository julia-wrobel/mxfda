#' univariate
#'
#' Internal function called by [extract_summary_functions()] to calculate a univariate spatial summary function for a single image.
#'
#' @param mximg Dataframe of cell-level multiplex imaging data for a single image.
#' Should have variables \code{x} and \code{y} to denote x and y spatial locations of each cell.
#' @param markvar The name of the variable that denotes cell type(s) of interest. Character.
#' @param mark1 dummy filler, unused
#' @param mark2 dummy filler, unused
#' @param r_vec Numeric vector of radii over which to evaluate spatial summary functions. Must begin at 0.
#' @param func Spatial summary function to calculate. Options are c(Kest, Lest, Gest) which denote Ripley's K, Besag's L, and nearest neighbor G function, respectively.
#' @param edge_correction Character string that denotes the edge correction method for spatial summary function. For Kest and Lest choose one of c("border", "isotropic", "Ripley", "translate", "none"). For Gest choose one of c("rs", "km", "han")
#' @param empirical_CSR logical to indicate whether to use the permutations to identify the sample-specific complete spatial randomness (CSR) estimation.
#' @param permutations integer for the number of permtuations to use if empirical_CSR is `TRUE` and exact CSR not calculable
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @references
#'
#' Creed, J. H., Wilson, C. M., Soupir, A. C., Colin-Leitzinger, C. M., Kimmel, G. J.,
#' Ospina, O. E., Chakiryan, N. H., Markowitz, J., Peres, L. C., Coghill, A., & Fridley, B. L. (2021).
#' spatialTIME and iTIME: R package and Shiny application for visualization and analysis of
#' immunofluorescence data. \emph{Bioinformatics} (Oxford, England), 37(23), 4584–4586.
#' https://doi.org/10.1093/bioinformatics/btab757
#'
#' @return A \code{data.frame} containing:
#' \item{r}{the radius of values over which the spatial summary function is evaluated}
#' \item{sumfun}{the values of the spatial summary function}
#' \item{csr}{the values of the spatial summary function under complete spatial randomness}
#' \item{fundiff}{sumfun - csr, positive values indicate clustering and negative values repulsion}
#'
#' @author Julia Wrobel \email{`r juliawrobel_email`}
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @importFrom spatstat.geom ppp convexhull.xy
#' @import dplyr
#'
#' @export
univariate = function(mximg,
                      markvar,
                      mark1,
                      mark2,
                      r_vec,
                      func = c(Kest, Lest, Gest),
                      edge_correction,
                      empirical_CSR = FALSE,
                      permutations = 1000){
  #### note to switch edge correction based on choice of func, this should be automated

  if(nrow(mximg) < 3) return(NA)

  # create a convex hull as the observation window
  w = convexhull.xy(mximg[["x"]], mximg[["y"]])

  # create ppp object
  pp_obj = ppp(mximg[["x"]], mximg[["y"]], window = w,
               marks = mximg[[markvar]], checkdup = FALSE)

  # estimate spatial summary function
  #since removed filtering in the extract_summary_function function, must apply filtering to
  sumfun = func(subset(pp_obj, marks == mark1), #subset ppp to mark1
                r = r_vec,
                correction = edge_correction)

  #if the user wants permutations to estimate CSR then will have to run them
  if(empirical_CSR == TRUE){
    if(identical(Kest, func)){
      #message("Using Exact Complete Spatial Randomness for Ripley's K")
      sumfun[[2]] = func(pp_obj, #exact CSR is across all points
                       r = r_vec,
                       correction = edge_correction)[[3]]
    } else {
      sumfun[[2]] = sapply(1:permutations, function(x){
        func(pp_obj[sample(1:pp_obj$n, sum(pp_obj$marks == mark1))], #sample ppp to length of the correct marks
             r = r_vec,
             correction = edge_correction)[[3]]
      }) %>%
        rowMeans(na.rm = TRUE) #take average of all permutations
    }

  }

  if(edge_correction == "none") colnames(sumfun)[3] = "none"

  df = data.frame(sumfun) %>%
    select(r, sumfun = all_of(edge_correction), csr = theo) %>%
    mutate(fundiff = sumfun - csr) %>%
    select(r, sumfun, csr, fundiff)

  df

}

