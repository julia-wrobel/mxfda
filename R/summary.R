#' Summary method for object of class `mxFDA`
#'
#' @param x object of class `mxFDA`
#' @param ... unused currently
#'
#' @return summary of object
#' @export
#'
#' @author Alex Soupir \email{alex.soupir@@moffitt.org}
#'
#' @examples
#' #set seed
#' set.seed(333)
#'
summary.mxFDA = function(x, ...){
  #basic object
  cat("mxFDA Object:\n")
  cat("\tSubjects: ", length(unique(x@Metadata[[x@subject_key]])), "\n", sep = "")
  cat("\tSamples: ", length(unique(x@Metadata[[x@sample_key]])), "\n", sep = "")
  #alert about spatial slot
  if(nrow(x@Spatial) == 0){
    cat("\tSpatial slot empty\n")
  } else {
    cat("\tHas spatial data\n")
  }
  #spatial summaries
  uni_sums = names(x@univariate_summaries)
  bi_sums = names(x@bivariate_summaries)
  cat("\tUnivariate Summaries: ", ifelse(length(uni_sums) == 0, "None", paste0(uni_sums, collapse =  ", ")), "\n", sep = "")
  cat("\tBivariate Summaries: ", ifelse(length(bi_sums) == 0, "None", paste0(bi_sums, collapse =  ", ")), "\n", sep = "")
  #any FPCs calculated
  if(length(x@functional_pca) == 0){
    cat("FPCs not yet calculated\n", sep = "")
  } else {
    fpc_slots = names(x@functional_pca)
    cat("FPCs Calculated:\n", sep = "")
    for(f in fpc_slots){
      cat("\t", f, ": ", ncol(x@functional_pca[[f]]$score_df), " FPCs describe ",
          round((ovarian_FDA@functional_pca$Gest$fpc_object$pve * 100), digits = 1), "% variance\n", sep = "")
    }
  }
  #mixed pricipal components
  if(length(x@functional_mpca) == 0){
    cat("MFPCs not yet calculated\n", sep = "")
  } else {
    fpc_slots = names(x@functional_mpca)
    cat("MFPCs Calculated:\n", sep = "")
    for(f in fpc_slots){
      cat("\t", f, ": ",
          ncol(x@functional_mpca[[f]]$score_df), " Level1 MFPCs and ",
          ncol(x@functional_mpca[[f]]$score_df), " Level2 MFPCs", "\n", sep = "") #need to play with the output to determine how to report
    }
  }
  #any models run
  if(length(x@functional_cox) == 0){
    cat("FCMs not yet calculated\n", sep = "")
  } else {
    f_cox_slots = names(x@functional_cox)
    cat("Models Fit:\n", sep = "")
    for(f in f_cox_slots){
      cat("\t", f, ": ", paste0(sapply(x@functional_cox[[f]], function(i){ class(i)[1]}) %>% toupper(),
                                collapse = ", "),
          " models\n", sep = "")
    }
  }
  #any models run
  if(length(x@functional_mcox) == 0){
    cat("MFCMs not yet calculated\n", sep = "")
  } else {
    f_cox_slots = names(x@functional_mcox)
    cat("Models Fit:\n", sep = "")
    for(f in f_cox_slots){
      cat("\t", f, ": ", paste0(sapply(x@functional_mcox[[f]], function(i){ class(i)[1]}) %>% toupper(),#need to play with the output to determine how to report
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
