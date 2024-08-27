#' bivariate
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
#' @param func Spatial summary function to calculate. Options are c(Kcross, Lcross, Gcross) which denote Ripley's K, Besag's L, and nearest neighbor G function, respectively, or entropy from Vu et al, 2023.
#' @param edge_correction Character string that denotes the edge correction method for spatial summary function. For Kcross and Lcross choose one of c("border", "isotropic", "Ripley", "translate", "none"). For Gcross choose one of c("rs", "km", "han")
#' @param breaks an integer for the number of breaks used for entropy
#' @param permute_CSR logical to indicate whether to use the permutations to identify the sample-specific complete spatial randomness (CSR) estimation.
#' @param permutations integer for the number of permtuations to use if permute_CSR is `TRUE` and exact CSR not calculable
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
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @references
#'
#' Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
#' Fast covariance estimation for high-dimensional functional data.
#' \emph{Statistics and Computing}, 26, 409-421.
#' DOI: 10.1007/s11222-014-9485-x.
#'
#' Vu, T., Seal, S., Ghosh, T., Ahmadian, M., Wrobel, J., & Ghosh, D. (2023).
#' FunSpace: A functional and spatial analytic approach to cell imaging data using entropy measures.
#' \emph{PLOS Computational Biology}, 19(9), e1011490.
#'
#' Creed, J. H., Wilson, C. M., Soupir, A. C., Colin-Leitzinger, C. M., Kimmel, G. J.,
#' Ospina, O. E., Chakiryan, N. H., Markowitz, J., Peres, L. C., Coghill, A., & Fridley, B. L. (2021).
#' spatialTIME and iTIME: R package and Shiny application for visualization and analysis of
#' immunofluorescence data. \emph{Bioinformatics} (Oxford, England), 37(23), 4584–4586.
#' https://doi.org/10.1093/bioinformatics/btab757
#'
#' @importFrom spatstat.geom ppp convexhull.xy
#' @import dplyr
#'
#' @export
bivariate = function(mximg,
                             markvar,
                             mark1,
                             mark2,
                             r_vec,
                             func = c(Kcross, Lcross, Gcross, entropy),
                             edge_correction,
                             breaks = NULL,
                     permute_CSR = FALSE,
                     permutations = 1000){
  #### note to switch edge correction based on choice of func, this should be automated

  # subset data to only cells if interest
  #some of the samples will have other cell types, as well as missing cell type of interest
  n = mximg %>% count(get(markvar)) %>% pull(n, name = 1)
  if(any(n < 2)) return(NA)
  if(sum(names(n) %in% c(mark1, mark2)) < 2) return(NA) #make sure the marks wanted are in n vector

  # create a convex hull as the observation window
  w = convexhull.xy(mximg[["x"]], mximg[["y"]])

  # create ppp object
  pp_obj = ppp(mximg[["x"]], mximg[["y"]], window = w,
               marks = as.factor(mximg[[markvar]]), checkdup = FALSE)
  if(!is.null(breaks)){
    r_vec = exp(seq(log(0.05*max(r_vec)), log(max(r_vec)), length.out = breaks))
  }
  # estimate L using spatstat
  sumfun = func(pp_obj,
                mark1,
                mark2,
                r = r_vec,
                correction = edge_correction)

  #if the user wants permutations to estimate CSR then will have to run them
  if(permute_CSR == TRUE){
    if(identical(Kcross, func)){
      #message("Using Exact Complete Spatial Randomness for Ripley's K")
      sumfun[[2]] = func(pp_obj, #exact CSR is across all points
                         r = r_vec,
                         correction = edge_correction)[[3]]
    } else {
      sumfun[[2]] = sapply(1:permutations, function(x){
        pp_obj_tmp = pp_obj
        pp_obj_tmp$marks = sample(pp_obj_tmp$marks)
        func(pp_obj_tmp, #sample ppp to length of the correct marks
             mark1,
             mark2,
             r = r_vec,
             correction = edge_correction)[[3]]
      }) %>%
        rowMeans(na.rm = TRUE) #take average of all permutations
    }

  }

  if(is.null(breaks)){
    df = data.frame(sumfun) %>%
      select(r, sumfun = all_of(edge_correction), csr = theo) %>%
      mutate(fundiff = sumfun - csr) %>%
      select(r, sumfun, csr, fundiff)
    return(df)
  } else {
    return(data.frame(sumfun))
  }



}




