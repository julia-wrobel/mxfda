#' Add Summary Function
#'
#' Sometimes other ways of calculating summary functions is wanted and is done in other packages,
#' in this instance the data can be loaded into the `mxFDA` object.
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param summary_function_data data frame with `summary_key` from `mxFDA` object as key column for summary function
#' @param metric character vector with either 'uni' or 'bi' and 'k', 'l', or 'g'; e.g. 'uni g'
#'
#' @return an updated `mxFDA` object with a derived value added
#' @export
#'
add_summary_function = function(mxFDAobject,
                                summary_function_data,
                                metric){
  #get the right data
  if(length(metric) != 1)
    stop("Please provide a single spatial metric to calculate functional PCA with")
  metric = unlist(strsplit(metric, split = " "))

  if(!(mxFDAobject@sample_key %in% colnames(summary_function_data)))
    stop("summary_function_data must have sample_key column name")

  #fill slot with new data
  if(grepl("[B|b]", metric[1]) & grepl("[K|k]", metric[2])) mxFDAobject@bivariate_summaries$Kcross = summary_function_data
  if(grepl("[B|b]", metric[1]) & grepl("[G|g]", metric[2])) mxFDAobject@`bivariate_summaries`$Gcross = summary_function_data
  if(grepl("[B|b]", metric[1]) & grepl("[L|l]", metric[2])) mxFDAobject@`bivariate_summaries`$Lcross = summary_function_data
  if(grepl("[U|u]", metric[1]) & grepl("[K|k]", metric[2])) mxFDAobject@`univariate_summaries`$Kest = summary_function_data
  if(grepl("[U|u]", metric[1]) & grepl("[G|g]", metric[2])) mxFDAobject@`univariate_summaries`$Gest = summary_function_data
  if(grepl("[U|u]", metric[1]) & grepl("[L|l]", metric[2])) mxFDAobject@`univariate_summaries`$Lest = summary_function_data

  #return updated object
  return(mxFDAobject)
}
