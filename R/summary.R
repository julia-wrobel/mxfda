#' Summary method for object of class `mxFDA`
#'
#' @param object object of class `mxFDA`
#' @param ... unused currently
#'
#' @details `r lifecycle::badge('stable')`
#'
#' @return summary of object to the R console
#'
#' @author Alex Soupir \email{`r alexsoupir_email`}
#'
#' @export
summary.mxFDA = function(object, ...){
  #basic object
  cat("mxFDA Object:\n")
  cat("\tSubjects: ", length(unique(object@Metadata[[object@subject_key]])), "\n", sep = "")
  cat("\tSamples: ", length(unique(object@Metadata[[object@sample_key]])), "\n", sep = "")
  #alert about spatial slot
  if(nrow(object@Spatial) == 0){
    cat("\tSpatial slot empty\n")
  } else {
    cat("\tHas spatial data\n")
  }
  #spatial summaries
  uni_sums = names(object@univariate_summaries)
  bi_sums = names(object@bivariate_summaries)
  cat("\tUnivariate Summaries: ", ifelse(length(uni_sums) == 0, "None", paste0(uni_sums, collapse =  ", ")), "\n", sep = "")
  cat("\tBivariate Summaries: ", ifelse(length(bi_sums) == 0, "None", paste0(bi_sums, collapse =  ", ")), "\n", sep = "")
  #any FPCs calculated
  if(length(object@functional_pca) == 0){
    cat("FPCs not yet calculated\n", sep = "")
  } else {
    fpc_slots = names(object@functional_pca)
    cat("FPCs Calculated:\n", sep = "")
    for(f in fpc_slots){
      cat("\t", f, ": ", ncol(object@functional_pca[[f]]$score_df), " FPCs describe ",
          round((object@functional_pca[[f]]$fpc_object$pve * 100), digits = 1), "% variance\n", sep = "")
    }
  }
  #mixed pricipal components
  if(length(object@functional_mpca) == 0){
    cat("MFPCs not yet calculated\n", sep = "")
  } else {
    fpc_slots = names(object@functional_mpca)
    cat("MFPCs Calculated:\n", sep = "")
    for(f in fpc_slots){
      cat("\t", f, ": ",
          ncol(object@functional_mpca[[f]]$score_df), " Level1 MFPCs and ",
          ncol(object@functional_mpca[[f]]$scores_level2), " Level2 MFPCs explaining ",
          round((object@functional_mpca[[f]]$mfpc_object$pve * 100), digits = 1),
          "% variance\n", sep = "") #need to play with the output to determine how to report
    }
  }
  #any models run
  if(length(object@functional_cox) == 0){
    cat("FCMs not yet calculated\n", sep = "")
  } else {
    f_cox_slots = names(object@functional_cox)
    cat("Models Fit:\n", sep = "")
    for(f in f_cox_slots){
      cat("\t", f, ": ", paste0(sapply(object@functional_cox[[f]], function(i){ class(i)[1]}) %>% toupper(),
                                collapse = ", "),
          " models\n", sep = "")
    }
  }
  #any models run
  if(length(object@functional_mcox) == 0){
    cat("MFCMs not yet calculated\n", sep = "")
  } else {
    f_cox_slots = names(object@functional_mcox)
    cat("Models Fit:\n", sep = "")
    for(f in f_cox_slots){
      cat("\t", f, ": ", paste0(sapply(object@functional_mcox[[f]], function(i){ class(i)[1]}) %>% toupper(),#need to play with the output to determine how to report
                                collapse = ", "),
          " models\n", sep = "")
    }
  }
  if(length(object@scalar_on_functional) == 0){
    cat("Scalar on Functional Regression not calculated\n", sep = "")
  } else {
    sofr_slots = names(object@scalar_on_functional)
    cat("Models Fit:\n", sep = "")
    for(f in sofr_slots){
      cat("\t", f, ": ", paste0(sapply(object@scalar_on_functional[[f]], function(i){ class(i)[1]}) %>% toupper(),#need to play with the output to determine how to report
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
