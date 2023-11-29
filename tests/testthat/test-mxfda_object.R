#test outputs of creating the mxFDA object
test_that("make_mxfda", {
  clinical = lung_df %>%
    dplyr::select(image_id, patient_id, patientImage_id, gender, age, survival_days, survival_status, stage) %>%
    dplyr::distinct()

  spatial = lung_df %>%
    select(-image_id, -gender, -age, -survival_days, -survival_status, -stage)
  sample_id_column = "patientImage_id"

  mxFDAobject = make_mxfda(metadata = clinical,
                           spatial = spatial,
                           subject_key = "patient_id",
                           sample_key = sample_id_column)

  #should be of right class
  #makes sure that while adding slots something didn't get goofed
  expect_equal(class(mxFDAobject)[1], "mxFDA")
  #clinical and spatial both need the sample_key column
  expect_error(make_mxfda(metadata = clinical %>% dplyr::select(-!!sample_id_column),
                          spatial = spatial,
                          subject_key = "patient_id",
                          sample_key = sample_id_column))
  expect_error(make_mxfda(metadata = clinical,
                          spatial = spatial %>% dplyr::select(-!!sample_id_column),
                          subject_key = "patient_id",
                          sample_key = sample_id_column))
  expect_error(make_mxfda(metadata = clinical,
                          spatial = spatial,
                          sample_key = sample_id_column))
  expect_error(make_mxfda(metadata = clinical,
                          spatial = spatial,
                          subject_key = "patient_id"))
  #should work without the need for spatial information
  expect_no_error(make_mxfda(metadata = clinical %>% dplyr::select(-!!sample_id_column),
                             spatial = NULL,
                             subject_key = "patient_id",
                             sample_key = sample_id_column))
})
