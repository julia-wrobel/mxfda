#' extract_bivariate
#'
#' Internal function called by \code{extract_summary_functions} to calculate a bivariate spatial summary function for a single image.
#'
#'
#' @param mximg Dataframe of cell-level multiplex imaging data for a single image.
#' Should have variables \code{x} and \code{y} to denote x and y spatial locations of each cell.
#' @param markvar The name of the variable that denotes cell type(s) of interest. Character.
#' @param mark1 Character string that denotes first cell type of interest.
#' @param mark2 Character string that denotes second cell type of interest.
#' @param r_vec Numeric vector of radii over which to evaluate spatial summary functions. Must begin at 0.
#' @param func Spatial summary function to calculate. Options are c(Kcross, Lcross, Gcross) which denote Ripley's K, Besag's L, and nearest neighbor G function, respectively.
#' @param edge_correction Character string that denotes the edge correction method for spatial summary function. For Kcross and Lcross choose one of c("border", "isotropic", "Ripley", "translate", "none"). For Gcross choose one of c("rs", "km", "han")
#'
#' @return A \code{data.frame} containing:
#' \item{r}{the radius of values over which the spatial summary function is evaluated}
#' \item{sumfun}{the values of the spatial summary function}
#' \item{csr}{the values of the spatial summary function under complete spatial randomness}
#' \item{fundiff}{sumfun - csr, positive values indicate clustering and negative values repulsion}
#'
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
#'
#' @references Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
#' Fast covariance estimation for high-dimensional functional data.
#' \emph{Statistics and Computing}, 26, 409-421.
#' DOI: 10.1007/s11222-014-9485-x.
#'
#' @importFrom spatstat.geom ppp convexhull.xy
#' @importFrom tibble as_tibble
#' @import dplyr
#'
#' @examples
#' # simulate data
#' set.seed(1001)
#'
#' @export
extract_bivariate = function(mximg,
                              markvar,
                              mark1,
                              mark2,
                              r_vec,
                              func = c(Kcross, Lcross, Gcross),
                              edge_correction){
  #### note to switch edge correction based on choice of func, this should be automated

  # subset data to only cells if interest
  #mximg = filter(mximg, get(markvar) %in% c(mark1, mark2))
  # mximg = filter(mximg, markvar == mark1)
  n = mximg %>% count(get(markvar)) %>% pull(n)
  if(any(n < 2)) return(NA)
  if(length(n) < 2) return(NA)

  # create a convex hull as the observation window
  w = convexhull.xy(mximg[["x"]], mximg[["y"]])

  # create ppp object
  pp_obj = ppp(mximg[["x"]], mximg[["y"]], window = w,
               marks = as.factor(mximg[[markvar]]), checkdup = FALSE)

  # estimate L using spatstat
  sumfun = func(pp_obj, mark1, mark2,
                r = r_vec,
                correction = edge_correction)

  df = as_tibble(sumfun) %>%
    select(r, sumfun = all_of(edge_correction), csr = theo) %>%
    mutate(fundiff = sumfun - csr) %>%
    select(r, sumfun, csr, fundiff)

  df

}




