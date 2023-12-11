#test outputs of creating the mxFDA object
#data processesing for tests
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

test_that("output data same dims", {
  mxFDAobject = filter_spatial(mxFDAobject, immune == "other")
  expect_equal(nrow(mxFDAobject@Spatial), nrow(spatial %>% dplyr::filter(immune == "other")))
})

test_that("correct class", {
  mxFDAobject = filter_spatial(mxFDAobject, immune == "other")
  expect_equal(class(mxFDAobject)[1], "mxFDA")
})

test_that("column not in spatial", {
  expect_error(filter_spatial(mxFDAobject, `T cell` == 1))
})

test_that("missing level in column - no force", {
  expect_error(filter_spatial(mxFDAobject, immune == "T cells"))
})

test_that("missing level in column - force", {
  mxFDAobject = filter_spatial(mxFDAobject, immune == "T cell", force = TRUE)
  expect_equal(nrow(mxFDAobject@Spatial), 0)
})
