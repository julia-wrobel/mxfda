#' Summary method for object of class `mxFDA`
#'
#' @param x object of class `mxFDA`
#' @param ... unused currently
#'
#' @return summary of object
#' @export
#'
summary.mxFDA = function(x, ...){
  #basic object
  cat("mxFDA Object:\n")
  cat("\tSubjects: ", length(unique(x@Metadata[[x@subject_key]])), "\n", sep = "")
  cat("\tSamples: ", length(unique(x@Metadata[[x@sample_key]])), "\n", sep = "")
  #any FPCs calculated
  if(length(x@functional_pca) == 0){
    cat("FPCs not yet calculated\n", sep = "")
  } else {
    fpc_slots = names(x@functional_pca)
    cat("FPCs Calculated:\n", sep = "")
    for(f in fpc_slots){
      cat("\t", f, ": ", ncol(x@functional_pca[[f]]$score_df), " FPCs\n", sep = "")
    }
  }
  #any models run
  if(length(x@functional_cox) == 0){
    cat("FPCs not yet calculated\n", sep = "")
  } else {
    f_cox_slots = names(x@functional_cox)
    cat("Models Fit:\n", sep = "")
    for(f in f_cox_slots){
      cat("\t", f, ": ", paste0(sapply(x@functional_cox[[f]], function(i){ class(i)[1]}) %>% toupper(),
                                collapse = ", "),
          " models\n", sep = "")
    }
  }
}

setMethod(f = 'show',
          signature = 'mxFDA',
          definition = function(object){
            summary(object)
          })
