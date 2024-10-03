#' Extract Summary Functions
#'
#' Function to extract spatial summary functions from the `Spatial` slot of an `mxFDA` object
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param r_vec Numeric vector of radii over which to evaluate spatial summary functions. Must begin at 0.
#' @param extract_func Defaults to univariate, which calculates univariate spatial summary functions. Choose bivariate for bivariate spatial summary functions.
#' @param summary_func Spatial summary function to calculate. Options are c(Kest, Lest, Gest) which denote Ripley's K, Besag's L, and nearest neighbor G function, respectively.
#' @param markvar The name of the variable that denotes cell type(s) of interest. Character.
#' @param mark1 Character string that denotes first cell type of interest.
#' @param mark2 Character string that denotes second cell type of interest for calculating bivariate summary statistics. Not used when calculating univariate statistics.
#' @param edge_correction Character string that denotes the edge correction method for spatial summary function. For Kest and Lest choose one of c("border", "isotropic", "Ripley", "translate", "none"). For Gest choose one of c("rs", "km", "han")
#' @param emperical_CSR logical to indicate whether to use the permutations to identify the sample-specific complete spatial randomness (CSR) estimation. If there are not enough levels present in `markvar` column for permutations, the theoretical will be used.
#' @param permutations integer for the number of permtuations to use if emperical_CSR is `TRUE` and exact CSR not calculable
#'
#' @details `r lifecycle::badge('stable')`
#'
#' Complete spatial randomness (CSR) is the estimation or measure of a spatial summary function when the points or cells in a sample are randomly distributed,
#' following no clustering or dispersion pattern. Some samples do have artifacts that may influence what CSR is under the distribution of points as they
#' are found in the sample such as large regions of missing points or possibly in the case of tissue sections, necrotic tissue where cells are dead. Theoretical
#' CSR requires points have an equal chance of occurring anywhere in the sample that these artifacts violate, necessitating the need to estimate or
#' calculate what this CSR would be for each sample independently. Previously Wilson et al. had demonstrated cases in which sample-specific
#' CSR was important over the use of the theoretical in calculating how much the observed deviates from expected.
#'
#' @return an object of class `mxFDA` containing the corresponding spatial summary function slot filled. See [make_mxfda()] for object structure details.
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
#' Wilson, C., Soupir, A. C., Thapa, R., Creed, J., Nguyen, J., Segura, C. M.,
#' Gerke, T., Schildkraut, J. M., Peres, L. C., & Fridley, B. L. (2022).
#' Tumor immune cell clustering and its association with survival in African
#' American women with ovarian cancer. PLoS computational biology, 18(3),
#' e1009900. https://doi.org/10.1371/journal.pcbi.1009900
#'
#' Creed, J. H., Wilson, C. M., Soupir, A. C., Colin-Leitzinger, C. M., Kimmel, G. J.,
#' Ospina, O. E., Chakiryan, N. H., Markowitz, J., Peres, L. C., Coghill, A., & Fridley, B. L. (2021).
#' spatialTIME and iTIME: R package and Shiny application for visualization and analysis of
#' immunofluorescence data. \emph{Bioinformatics} (Oxford, England), 37(23), 4584–4586.
#' https://doi.org/10.1093/bioinformatics/btab757
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
#'                                         extract_func = univariate,
#'                                         summary_func = Kest,
#'                                         markvar = "immune",
#'                                         mark1 = "immune",
#'                                         edge_correction = "trans")
#'
#' @export
extract_summary_functions <- function(mxFDAobject, r_vec = seq(0, 100, by = 10),
                                      extract_func = c(univariate, bivariate),
                                      summary_func = c(Kest, Lest, Gest),
                                      markvar,
                                      mark1,
                                      mark2 = NULL,
                                      edge_correction,
                                      emperical_CSR = FALSE,
                                      permutations = 1000
                                      ){
  if(!inherits(mxFDAobject, "mxFDA"))
    stop("Object must be of class `mxFDA`.")
  #need spatial data to calculate spatial summary functions
  if(nrow(mxFDAobject@Spatial) == 0)
    stop("No summary function to be calculated - missing spatial")
  #if running entropy, break out and do different
  if(identical(summary_func, entropy))
    stop("To calculate entropy, please use `extract_entropy()`")
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

  markvar_levels = length(unique(mxFDAobject@Spatial[[markvar]]))

  #check if enough marks to permute with if needed
  emperical_CSR = can_permute(extract_func, summary_func, emperical_CSR, markvar_levels)

  df_nest = mxFDAobject@Spatial %>%
    select(all_of(mxFDAobject@sample_key), x, y, all_of(markvar)) %>%
    #need all points to be sent to sub functions to calculate sample specific CSR when emperical_CSR == TRUE
    #filter(get(markvar) %in% c(mark1, mark2)) %>%
    nest(data = c(x, y, all_of(markvar)))

  ndat = df_nest %>% mutate(sumfuns = map(df_nest$data, extract_func,
                                          markvar = markvar,
                                          mark1 = mark1,
                                          mark2 = mark2,
                                          r_vec = r_vec,
                                          func = summary_func,
                                          edge_correction = edge_correction,
                                          emperical_CSR = emperical_CSR,
                                          permutations = permutations,
                                          .progress = TRUE)) %>%
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

   if(deparse(substitute(extract_func)) == "univariate"){
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


