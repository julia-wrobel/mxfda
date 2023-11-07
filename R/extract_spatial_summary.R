#' Summarise spatial data in mxFDA object
#'
#' @param mxFDAobject object of class `mxFDA`
#' @param columns character vector for column heading for cells to summarise
#' @param grouping_columns character vector of other columns to use as grouping, such as region classification column
#'
#' @return data frame with percent of total points per spatial sample `column` are positive for
#' @export
#'
extract_spatial_summary = function(mxFDAobject, columns, grouping_columns = NULL){
  #set function specific dplyr option to mute summarise
  options(dplyr.summarise.inform = FALSE)
  #check classes
  if(!inherits(mxFDAobject, "mxFDA"))
    stop('Need to use object of `mxFDA` class')
  #extract data
  dat = mxFDAobject@Spatial
  if(FALSE %in% (columns %in% colnames(dat)))
    stop('Need to provide at least one cell type column to summarise')
  if(FALSE %in% (grouping_columns %in% colnames(dat)) & !is.null(grouping_columns))
    stop('If providing other grouping columns, need to be in data')
  #check class of columns to see if user did 1/0s or names
  col_classes = apply(dat[,grouping_columns], 2, class)
  if(!("character" %in% col_classes))
    stop("Not yet working with factor 1/0 for cell classes")

  dat2 = dat %>%
    dplyr::group_by(dplyr::across(mxFDAobject@sample_key)) %>%
    dplyr::mutate(`Total Cells` = dplyr::n())

  summ_1 = lapply(columns, function(col_of_interest){
    dat2 %>%
      dplyr::select(dplyr::any_of(c(mxFDAobject@sample_key, grouping_columns, "Total Cells", col_of_interest))) %>%
      dplyr::group_by(dplyr::across(!!c(mxFDAobject@sample_key, grouping_columns, "Total Cells")), dplyr::across(!!col_of_interest)) %>%
      dplyr::summarise(Cells = dplyr::n()) %>%
      tidyr::spread(col_of_interest, "Cells", fill = 0)
  }) %>%
    purrr::reduce(dplyr::full_join, by = c(mxFDAobject@sample_key, grouping_columns, "Total Cells")) %>%
    dplyr::ungroup()

  summ_columns = setdiff(colnames(summ_1), c("Total Cells", grouping_columns, mxFDAobject@sample_key))

  spatial_sum = summ_1 %>%
    replace(is.na(.), 0) %>%
    dplyr::mutate(across(!!summ_columns, ~ .x / `Total Cells` * 100, .names = "{.col} Percent"))
  return(spatial_sum)
}
