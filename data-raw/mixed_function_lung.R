# library(dplyr)
# #making lung FDA with uni g for mFPCA example
# data(lung_df)
#
# #get the clinical information
# clinical = lung_df %>%
#   dplyr::select(image_id, patient_id, patientImage_id, gender,
#                 age, survival_days, survival_status, stage) %>%
#   dplyr::distinct()
#
# #make small, just need to make sure it runs
# spatial = lung_df %>%
#   dplyr::select(-image_id, -gender, -age, -survival_days, -survival_status, -stage)
#
# #create `mxFDA` object
# mxFDAobject = make_mxfda(metadata = clinical,
#                          spatial = spatial,
#                          subject_key = "patient_id",
#                          sample_key = "patientImage_id")
#
# #calculate summary function Gest
# mxFDAobject = extract_summary_functions(mxFDAobject, r_vec = 0:100, extract_func = extract_univariate, summary_func = Gest, markvar = 'phenotype_cd8', mark1 = 'CD8+', edge_correction = 'rs')
#
# plot(mxFDAobject, y = 'fundiff', what = 'uni g')
#
# #checking if there are samples that have no computable results
# spat_summ = mxFDAobject@univariate_summaries$Gest
# keep_samples = spat_summ %>%
#   dplyr::group_by(patientImage_id) %>%
#   dplyr::summarise(non_nan_count = sum(!is.nan(fundiff))) %>%
#   dplyr::filter(non_nan_count >= 10)
# #looks like all samples have something at small radii that turns to NaN at larger
#
# #filter clinical and spatial to include only those with enough non-missing values
# mxFDAobject@Metadata = mxFDAobject@Metadata %>%
#   dplyr::filter(patientImage_id %in% keep_samples$patientImage_id)
# mxFDAobject@univariate_summaries$Gest = mxFDAobject@univariate_summaries$Gest %>%
#   dplyr::filter(patientImage_id %in% keep_samples$patientImage_id)
# #remove the spatial data to shrink object
# mxFDAobject@Spatial = data.frame()
#
# #test it works
# test = run_mfpca(mxFDAobject, metric = 'uni g')
#
# lung_FDA = mxFDAobject
usethis::use_data(lung_FDA, overwrite = TRUE)
