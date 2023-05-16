#' extract_univariate
#'
#' Internal function called by \code{extract_summary_functions} to calculate a univariate spatial summary function for a single image.
#'
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
#' @importFrom spatstat.explore Kest Lest Gest
#' @importFrom spatstat.geom ppp convexhull.xy
#' @importFrom tibble as_tibble
#' @import dplyr
#'
#' @return A \code{data.frame} containing:
#' \item{r}{the radius of values over which the spatial summary function is evaluated}
#' \item{sumfun}{the values of the spatial summary function}
#' \item{csr}{the values of the spatial summary function under complete spatial randomness}
#' \item{fundiff}{sumfun - csr, positive values indicate clustering and negative values repulsion}
#'
#' @references Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
#' Fast covariance estimation for high-dimensional functional data.
#' \emph{Statistics and Computing}, 26, 409-421.
#' DOI: 10.1007/s11222-014-9485-x.

#' @examples
#' # simulate data
#' set.seed(1001)
#'
#'
#' @param mximg Dataframe of cell-level multiplex imaging data for a single image.
#' Should have variables \code{x} and \code{y} to denote x and y spatial locations of each cell.
#' @param markvar The name of the variable that denotes cell type(s) of interest. Character.
#' @param mark1 Character string that denotes cell type of interest.
#' @param r_vec Numeric vector of radii over which to evaluate spatial summary functions. Must begin at 0.
#' @param func Spatial summary function to calculate. Options are c(Kest, Lest, Gest) which denote Ripley's K, Besag's L, and nearest neighbor G function, respectively.
#' @param edge_correction Character string that denotes the edge correction method for spatial summary function. For Kest and Lest choose one of c("border", "isotropic", "Ripley", "translate", "none"). For Gest choose one of c("rs", "km", "han")
extract_univariate = function(mximg,
                              markvar,
                              mark1,
                              mark2,
                              r_vec,
                              func = c(Kest, Lest, Gest),
                              edge_correction){
  #### note to switch edge correction based on choice of func, this should be automated


  # subset data to only cells if interest
  # mximg = filter(mximg, markvar == mark1)

  if(nrow(mximg) < 3) return(NA)

  # create a convex hull as the observation window
  w = convexhull.xy(mximg[["x"]], mximg[["y"]])

  # create ppp object
  pp_obj = ppp(mximg[["x"]], mximg[["y"]], window = w, marks = mximg[[markvar]])

  # estimate L using spatstat
  sumfun = func(pp_obj,
                r = r_vec,
                correction = edge_correction)

  df = as_tibble(sumfun) %>%
    select(r, sumfun = all_of(edge_correction), csr = theo) %>%
    mutate(fundiff = sumfun - csr) %>%
    select(r, sumfun, csr, fundiff)

  df

}

