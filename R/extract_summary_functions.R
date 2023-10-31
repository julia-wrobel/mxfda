#' extract_summary_functions
#'
#' This is the main function for the fastGFPCA package.
#' The function requires input data \code{Y} to be a dataframe in long format with variables
#' \code{id}, \code{index}, and \code{value} to indicate subject IDs,
#' observation times on the domain, and observations, respectively.
#' The \code{index} must contain the same, equally spaced grid points for each subject.

#' The number of functional principal components (FPCs) can either be specified
#' directly (argument \code{npc}) or chosen based on the explained share of
#' variance (\code{pve}). Families and link functions can be specified using the same
#' syntax as the \code{family} argument from the \code{lme4::glmer()} function.
#'
#'
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
#' @import dplyr
#' @importFrom tidyr nest unnest
#' @importFrom purrr map
#'
#' @return A \code{data.frame} containing:
#' \item{image_patient_id}{the unique image id}
#' \item{r}{the radius of values over which the spatial summary function is evaluated}
#' \item{sumfun}{the values of the spatial summary function}
#' \item{csr}{the values of the spatial summary function under complete spatial randomness}
#' \item{fundiff}{sumfun - csr, positive values indicate clustering and negative values repulsion}
#' @references Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
#' Fast covariance estimation for high-dimensional functional data.
#' \emph{Statistics and Computing}, 26, 409-421.
#' DOI: 10.1007/s11222-014-9485-x.

#' @examples
#' # simulate data
#' set.seed(1001)
#'
#'
#' @param mxdata Dataframe of cell-level multiplex imaging data. Should have variables x and y to denote x and y spatial locations of each cell.
#' @param image_patient_id The name of the variable that identifies each unique image.
#' @param r_vec Numeric vector of radii over which to evaluate spatial summary functions. Must begin at 0.
#' @param extract_func Defaults to extract_univariate, which calculates univariate spatial summary functions. Choose extract_bivariate for bivariate spatial summary functions.
#' @param summary_func Spatial summary function to calculate. Options are c(Kest, Lest, Gest) which denote Ripley's K, Besag's L, and nearest neighbor G function, respectively.
#' @param markvar The name of the variable that denotes cell type(s) of interest. Character.
#' @param mark1 Character string that denotes first cell type of interest.
#' @param mark2 Character string that denotes second cell type of interest for calculating bivariate summary statistics. Not used when calculating univariate statistics.
#' @param edge_correction Character string that denotes the edge correction method for spatial summary function. For Kest and Lest choose one of c("border", "isotropic", "Ripley", "translate", "none"). For Gest choose one of c("rs", "km", "han")
#' @param analysis_vars Optional list of variables to be retained for downstream analysis.
#' @export
extract_summary_functions <- function(mxFDAobject, image_patient_id, r_vec = seq(0, 100, by = 10),
                                      extract_func = c(extract_univariate, extract_bivariate),
                                      summary_func = c(Kest, Lest, Gest),
                                      markvar,
                                      mark1,
                                      mark2 = NULL,
                                      edge_correction
                                      ){



  df_nest = mxFDAobject@Spatial %>%
    select(all_of(image_patient_id), x, y, all_of(markvar)) %>%
    filter(get(markvar) %in% c(mark1, mark2)) %>%
    nest(data = c(x, y, all_of(markvar)))

   ndat = df_nest %>% mutate(sumfuns = map(df_nest$data, extract_func,
                                    markvar = markvar,
                                    mark1 = mark1,
                                    mark2 = mark2,
                                    r_vec = r_vec,
                                    func = summary_func,
                                    edge_correction = edge_correction)) %>%
     select(-data) %>%
     unnest(sumfuns)

   if(deparse(substitute(extract_func)) == "extract_univariate"){
     if(deparse(substitute(summary_func)) == "Kest") mxFDAobject@`Univariate Summaries`$Kest = ndat
     if(deparse(substitute(summary_func)) == "Lest") mxFDAobject@`Univariate Summaries`$Lest = ndat
     if(deparse(substitute(summary_func)) == "Gest") mxFDAobject@`Univariate Summaries`$Gest = ndat
   } else {
     if(deparse(substitute(summary_func)) == "Kcross") mxFDAobject@`Bivariate Summaries`$Kcross = ndat
     if(deparse(substitute(summary_func)) == "Lcross") mxFDAobject@`Bivariate Summaries`$Lcross = ndat
     if(deparse(substitute(summary_func)) == "Gcross") mxFDAobject@`Bivariate Summaries`$Gcross = ndat
   }
   return(mxFDAobject)
}


