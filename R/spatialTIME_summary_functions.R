#' `spatialTIME` summary function extraction
#'
#' @param mxFDAobject object of class `mxFDA` for which to calculate summary functions
#' @param metric string vector with combination of uni/bi k/g
#' @param markers column names of spatial
#' @param r_range vector of radii to calculate metrics
#' @param num_permutations integer for the number of permutations to perform
#' @param edge_correction character vector of length 1 or 2 depending on the metrics for edge corrections to perform
#' @param permute logical whether to permute CSR or not. Required `TRUE` for Nearest Neighbor G
#' @param workers integer for number of cores to use to calculate derived variables
#' @param xloc,yloc columns in the spatial data which to use for point locations
#'
#' @details `r lifecycle::badge('experimental')`
#'
#' @returns `mxFDA` object with the summary function slots filled according to input metric
#'
#' @keywords internal
#'
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @examples
#' # example code
#'
#' @export
#'
spatialTIME_summary_functions = function(mxFDAobject,
                                         metric = 'uni k',
                                         markers = NULL,
                                         r_range = 0:100,
                                         num_permutations = 50,
                                         edge_correction = "trans", #parse out options for G and K if a vector and check wanted metric matches edge correction
                                         permute = TRUE,
                                         workers = 1,
                                         xloc = "x",
                                         yloc = "y"){
  if(!inherits(mxFDAobject, "mxFDA"))
    stop("mxFDAobject must be created witih make_mxfda()")
  if(is.null(mxFDAobject@Spatial))
    stop("To calculate summary functions, spatial cannot be empty")

  #split the metrics
  metrics = lapply(metric, function(x){
    m = unlist(strsplit(x, " "))
    if(!(grepl("[B|b]", m[1]) | grepl("[U|u]", m[1])))
      stop("metrics must be either univariate or bivariate")
    if(!(grepl("[K|k]", m[2]) | grepl("[G|g]", m[2])))
      stop("metrics must be either K or G")
    if(grepl("[G|g]", m[2]) & !permute)
      stop("Nearest Neighbor G requires permutation")
    #edge correction check
    if(grepl("[G|g]", m[2]))
      if(length(intersect(c('rs', 'han', 'none'), edge_correction)) == 0)
        stop("Need edge correction appropriate for Nearest Neighbor G(r)")
    #edge correction check
    if(grepl("[K|k]", m[2]))
      if(length(intersect(c('trans', 'none'), edge_correction)) == 0)
        stop("Need edge correction appropriate for Ripley's K(r)")
    #if passes checks, return the actual vector
    return(m)
  })

  #markers
  #checks if markers is null
  if(is.null(markers)){
    #tries to identify columns where there might be marker information
    markers = colnames(mxFDAobject@Spatial) %>%
      grep("pheno|cd|fox|\\+|\\-", ., value = TRUE)
    keep_col_vec = sapply(markers, function(x){
      get_marker_data(mxFDAobject@Spatial[[x]])
    }) %>% data.frame(check.names = FALSE)

  } else {
    #make sure that that there is some intersection of input and available
    markers = intersect(markers, colnames(mxFDAobject@Spatial))
    keep_col_vec = sapply(markers, function(x){
      get_marker_data(mxFDAobject@Spatial[[x]])
    }) %>% data.frame(check.names = FALSE)
  }
  if(is.null(markers))
    stop("Could not find markers in the spatial data")

  #convert spatial to list with sample name in first column
  #remove the previous marker columns and add in the data frame with 0/1 for -/+
  spat_df = mxFDAobject@Spatial %>%
    dplyr::relocate(!!mxFDAobject@sample_key, .before = 1) %>%
    dplyr::select(!dplyr::all_of(markers)) %>%
    dplyr::bind_cols(keep_col_vec)
  spat_list = split(spat_df, spat_df[[mxFDAobject@sample_key]])

  #tmp summary
  spatial_summary = extract_spatial_summary(mxFDAobject, markers) %>%
    dplyr::full_join(mxFDAobject@Metadata %>%
                       dplyr::select(dplyr::all_of(c(mxFDAobject@subject_key,
                                       mxFDAobject@sample_key))))

  mif = spatialTIME::create_mif(clinical_data = mxFDAobject@Metadata,
                                sample_data = spatial_summary,
                                spatial_list = spat_list,
                                patient_id = mxFDAobject@subject_key,
                                sample_id = mxFDAobject@sample_key)

  for(i in seq(metrics)){
    derived_var = metrics[[i]]
    #Univariate Ripley's K
    if(grepl("[U|u]", derived_var[1]) & grepl("[K|k]", derived_var[2])) {
       dat = spatialTIME::ripleys_k(mif = mif,
                                    mnames = markers,
                                    r_range = r_range,
                                    num_permutations = num_permutations,
                                    edge_correction =
                                      intersect(c('trans', 'none'), edge_correction),
                                    permute = permute,
                                    keep_permutation_distribution = FALSE,
                                    workers = workers,
                                    overwrite = TRUE,
                                    xloc = xloc,
                                    yloc = yloc)$derived$univariate_Count
      mxFDAobject@`univariate_summaries`$Kest = dat
    }
    #Univariate Nearest Neighbor G
    if(grepl("[U|u]", derived_var[1]) & grepl("[G|g]", derived_var[2])) {
      dat = spatialTIME::NN_G(mif = mif,
                              mnames = markers,
                              r_range = 0:100,
                              num_permutations = num_permutations,
                              edge_correction = intersect(c('rs', 'han', 'none'), edge_correction),
                              keep_perm_dis = FALSE,
                              workers = workers,
                              overwrite = TRUE,
                              xloc = xloc,
                              yloc = yloc)$derived$univariate_NN

      mxFDAobject@`univariate_summaries`$Gest = dat
    }
    #Bivariate Ripley's K
    if(grepl("[B|b]", derived_var[1]) & grepl("[K|k]", derived_var[2])){
      dat = spatialTIME::bi_ripleys_k(mif = mif,
                                      mnames = markers,
                                      r_range = r_range,
                                      num_permutations = num_permutations,
                                      edge_correction =
                                        intersect(c('trans', 'none'), edge_correction),
                                      permute = permute,
                                      keep_permutation_distribution = FALSE,
                                      workers = workers,
                                      overwrite = TRUE,
                                      xloc = xloc,
                                      yloc = yloc)$derived$bivariate_Count
      mxFDAobject@`univariate_summaries`$Kcross = dat
    }
    #Bivariate Nearest Neighbor G
    if(grepl("[B|b]", derived_var[1]) & grepl("[G|g]", derived_var[2])) {
      dat = spatialTIME::bi_NN_G(mif = mif,
                                 mnames = markers,
                                 r_range = 0:100,
                                 num_permutations = num_permutations,
                                 edge_correction = intersect(c('rs', 'han', 'none'), edge_correction),
                                 keep_perm_dis = FALSE,
                                 workers = workers,
                                 overwrite = TRUE,
                                 xloc = xloc,
                                 yloc = yloc)$derived$bivariate_NN
      mxFDAobject@`bivariate_summaries`$Gcross = dat
    }
  }

  return(mxFDAobject)

}

get_marker_data = function(column_data){
  if(one_zero(column_data))
    return(column_data)
  if(TRUE %in% grepl("\\+", column_data) & TRUE %in% grepl("\\-", column_data)) #if not 0/1, look for -/+
    return(ifelse(grepl("\\+", column_data), 1, 0))
  return(rep(0, length(column_data)))
}
