#' Extract Summary Functions
#'
#' Function to extract spatial summary functions from the `Spatial` slot of an `mxFDA` object
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param r_vec Numeric vector of radii over which to evaluate spatial summary functions. Must begin at 0.
#' @param extract_func Defaults to extract_univariate, which calculates univariate spatial summary functions. Choose extract_bivariate for bivariate spatial summary functions.
#' @param summary_func Spatial summary function to calculate. Options are c(Kest, Lest, Gest) which denote Ripley's K, Besag's L, and nearest neighbor G function, respectively.
#' @param markvar The name of the variable that denotes cell type(s) of interest. Character.
#' @param mark1 Character string that denotes first cell type of interest.
#' @param mark2 Character string that denotes second cell type of interest for calculating bivariate summary statistics. Not used when calculating univariate statistics.
#' @param edge_correction Character string that denotes the edge correction method for spatial summary function. For Kest and Lest choose one of c("border", "isotropic", "Ripley", "translate", "none"). For Gest choose one of c("rs", "km", "han")
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @return an object of class `mxFDA` containing the corresponding spatial summary function slot filled. See [make_mxfda()] for object structure details.
#'
#' @author Julia Wrobel \email{`r juliawrobel_email`}
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @references Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
#' Fast covariance estimation for high-dimensional functional data.
#' \emph{Statistics and Computing}, 26, 409-421.
#' DOI: 10.1007/s11222-014-9485-x.
#'
#' [spatstat.explore::Kest()]
#'
#' [spatstat.explore::Gest()]
#'
#' [spatstat.explore::Lest()]
#'
#' [spatstat.explore::Kcross()]
#'
#' [spatstat.explore::Gcross()]
#'
#' [spatstat.explore::Lcross()]
#'
#' @import dplyr
#' @importFrom tidyr nest unnest
#' @importFrom purrr map
#'
#' @examples
#' #load ovarian FDA object
#' data('ovarian_FDA')
#'
#' #run function
#' ovarian_FDA = extract_summary_functions(ovarian_FDA, r_vec = 0:100,
#'                                         extract_func = extract_univariate,
#'                                         summary_func = Kest,
#'                                         markvar = "immune",
#'                                         mark1 = "immune",
#'                                         edge_correction = "trans")
#'
#' @export
extract_summary_functions <- function(mxFDAobject, r_vec = seq(0, 100, by = 10),
                                      extract_func = c(extract_univariate, extract_bivariate),
                                      summary_func = c(Kest, Lest, Gest),
                                      markvar,
                                      mark1,
                                      mark2 = NULL,
                                      edge_correction
                                      ){
  if(!inherits(mxFDAobject, "mxFDA"))
    stop("Object must be of class `mxFDA`.")
  #need spatial data to calculate spatial summary functions
  if(nrow(mxFDAobject@Spatial) == 0)
    stop("No summary function to be calculated - missing spatial")
  #check correction methods
  k_l_correction = c("border", "isotropic", "Ripley", "translate")
  gest_correction = c("rs", "km", "han")
  #check function
  if(identical(all.equal(Kest, summary_func), TRUE) |
     identical(all.equal(Kcross, summary_func), TRUE) |
     identical(all.equal(Lest, summary_func), TRUE) |
     identical(all.equal(Lcross, summary_func), TRUE)){
    if(!any(grepl(edge_correction, k_l_correction)))#will capture short names like iso and trans
      stop("edge correction must match summary function")
  }
  if(identical(all.equal(Gest, summary_func), TRUE) |
     identical(all.equal(Gcross, summary_func), TRUE) ){
    if(!any(grepl(edge_correction, gest_correction)))
      stop("edge correction must match summary function")
  }

  df_nest = mxFDAobject@Spatial %>%
    select(all_of(mxFDAobject@sample_key), x, y, all_of(markvar)) %>%
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

   cell_counts = mxFDAobject@Spatial %>%
     dplyr::filter(get(markvar) %in% !!c(mark1, mark2)) %>%
     dplyr::group_by(across(!!c(mxFDAobject@sample_key, markvar))) %>%
     dplyr::summarise(counts = dplyr::n()) %>%
     #summarise number of cells for each class
     dplyr::mutate(!!markvar := paste0(get(markvar), " cells")) %>%
     tidyr::spread(key = markvar, value = 'counts')

   ndat = full_join(ndat, cell_counts, by = mxFDAobject@sample_key)

   if(deparse(substitute(extract_func)) == "extract_univariate"){
     if(deparse(substitute(summary_func)) == "Kest") mxFDAobject@`univariate_summaries`$Kest = ndat
     if(deparse(substitute(summary_func)) == "Lest") mxFDAobject@`univariate_summaries`$Lest = ndat
     if(deparse(substitute(summary_func)) == "Gest") mxFDAobject@`univariate_summaries`$Gest = ndat
   } else {
     if(deparse(substitute(summary_func)) == "Kcross") mxFDAobject@`bivariate_summaries`$Kcross = ndat
     if(deparse(substitute(summary_func)) == "Lcross") mxFDAobject@`bivariate_summaries`$Lcross = ndat
     if(deparse(substitute(summary_func)) == "Gcross") mxFDAobject@`bivariate_summaries`$Gcross = ndat
   }
   return(mxFDAobject)
}


