#' extract_univariate
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
#' @param breaks ignored
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @return A \code{data.frame} containing:
#' \item{r}{the radius of values over which the spatial summary function is evaluated}
#' \item{sumfun}{the values of the spatial summary function}
#' \item{csr}{the values of the spatial summary function under complete spatial randomness}
#' \item{fundiff}{sumfun - csr, positive values indicate clustering and negative values repulsion}
#'
#' @author Julia Wrobel \email{`r juliawrobel_email`}
#'
#' @importFrom spatstat.geom ppp convexhull.xy
#' @import dplyr
#'
#' @export
extract_univariate = function(mximg,
                              markvar,
                              mark1,
                              mark2,
                              r_vec,
                              func = c(Kest, Lest, Gest),
                              edge_correction,
                              breaks = NULL){
  #### note to switch edge correction based on choice of func, this should be automated

  if(nrow(mximg) < 3) return(NA)

  # create a convex hull as the observation window
  w = convexhull.xy(mximg[["x"]], mximg[["y"]])

  # create ppp object
  pp_obj = ppp(mximg[["x"]], mximg[["y"]], window = w,
               marks = mximg[[markvar]], checkdup = FALSE)

  # estimate spatial summary function
  sumfun = func(pp_obj,
                r = r_vec,
                correction = edge_correction)

  if(edge_correction == "none") colnames(sumfun)[3] = "none"

  df = data.frame(sumfun) %>%
    select(r, sumfun = all_of(edge_correction), csr = theo) %>%
    mutate(fundiff = sumfun - csr) %>%
    select(r, sumfun, csr, fundiff)

  df

}

