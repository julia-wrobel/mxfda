#building data for tests
data(lung_df)

clinical = lung_df %>%
  select(image_id, patient_id, patientImage_id, gender, age, survival_days, survival_status, stage) %>%
  distinct()

spatial = lung_df %>%
  select(-image_id, -gender, -age, -survival_days, -survival_status, -stage)

spatial2 = spatial %>%
  select(-X.1, -X) %>%
  mutate(across(phenotype_ck:phenotype_cd4, ~ ifelse(grepl("\\+", .x), 1, 0))) %>%
  relocate(patientImage_id, .before = 1)

mxFDAobject = make_mxfda(metadata = clinical,
                         spatial = spatial,
                         subject_key = "patient_id",
                         sample_key = "patientImage_id")
markers = colnames(spatial) %>%
  grep("phenotype", ., value = TRUE)

test_that("column in spatial", {
  expect_error(extract_spatial_summary(mxFDAobject,
                                       "CD3+CD8+"))
})

test_that("get data frame back", {
  expect_true(is(extract_spatial_summary(mxFDAobject,
                                         markers), "data.frame"))
})

#summary with 1/0 in columns
test_that("1/0 spatial", {
  mxFDAobject = make_mxfda(metadata = clinical,
                           spatial = spatial2,
                           subject_key = "patient_id",
                           sample_key = "patientImage_id")
  expect_error(x)
})
