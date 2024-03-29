#load data
data(lung_df)

clinical = lung_df %>%
  select(image_id, patient_id, patientImage_id, gender, age, survival_days, survival_status, stage) %>%
  distinct()
#make small, just need to make sure it runs
spatial = lung_df %>%
  select(-image_id, -gender, -age, -survival_days, -survival_status, -stage) %>%
  filter(patientImage_id %in% clinical$patientImage_id[1:10])

mxFDAobject = make_mxfda(metadata = clinical,
                         spatial = spatial,
                         subject_key = "patient_id",
                         sample_key = "patientImage_id"
)

#test running metrics
test_that("can run spatial summary", {
  mxFDAobject = extract_summary_functions(mxFDAobject,
                                          extract_func = univariate,
                                          summary_func = Kest,
                                          r_vec = seq(0, 100, by = 1),
                                          edge_correction = "iso",
                                          markvar = "immune",
                                          mark1 = "immune")
  expect_true(is(mxFDAobject, "mxFDA"))
})
#test running metrics
test_that("can run spatial summary", {
  mxFDAobject = extract_summary_functions(mxFDAobject,
                                          extract_func = bivariate,
                                          summary_func = Kcross,
                                          r_vec = seq(0, 100, by = 1),
                                          edge_correction = "iso",
                                          markvar = "immune",
                                          mark1 = "immune",
                                          mark2 = "other")
  expect_true(is(mxFDAobject, "mxFDA"))
})
#test missing spatial
test_that("can run spatial summary", {
  mxFDAobject@Spatial = data.frame()
  expect_error(extract_summary_functions(mxFDAobject,
                                          extract_func = univariate,
                                          summary_func = Kest,
                                          r_vec = seq(0, 100, by = 1),
                                          edge_correction = "iso",
                                          markvar = "immune",
                                          mark1 = "immune"))
})

#test right spot
test_that("can run spatial summary", {
  mxFDAobject = extract_summary_functions(mxFDAobject,
                                          extract_func = bivariate,
                                          summary_func = Kcross,
                                          r_vec = seq(0, 100, by = 1),
                                          edge_correction = "iso",
                                          markvar = "immune",
                                          mark1 = "immune",
                                          mark2 = "other")
  expect_true(is(mxFDAobject@bivariate_summaries$Kcross, "data.frame"))
})
