data("ovarian_FDA")

ovarian_FDA = run_fcm(ovarian_FDA, model_name = "fit_lfcm",
                      formula = survival_time ~ age, event = "event",
                      metric = "uni g", r = "r", value = "fundiff",
                      analysis_vars = c("age", "survival_time"),
                      afcm = FALSE)

#keeping the code but commenting it out.
#theres something with the dplyr selects that still reports i'm not using all_of
test_that("output is correct class", {
  # ovarian_FDA = run_fcm(ovarian_FDA, model_name = "fit_lfcm",
  #                       formula = survival_time ~ age, event = "event",
  #                       metric = "uni g", r = "r", value = "fundiff",
  #                       analysis_vars = c("age", "survival_time"),
  #                       afcm = FALSE)
  expect_equal(class(ovarian_FDA)[1], "mxFDA")
})

test_that("metric needs to be single element vector", {
  expect_error(run_fcm(ovarian_FDA, model_name = "fit_lfcm",
                       formula = survival_time ~ age, event = "event",
                       metric = c("uni g", "uni k"), r = "r", value = "fundiff",
                       analysis_vars = c("age", "survival_time"),
                       afcm = FALSE))
})

test_that("metric needs to have summary function calculate", {
  expect_error(run_fcm(ovarian_FDA, model_name = "fit_lfcm",
                       formula = survival_time ~ age, event = "event",
                       metric = "uni k", r = "r", value = "fundiff",
                       analysis_vars = c("age", "survival_time"),
                       afcm = FALSE))
})

test_that("event needs to be a column in metadata", {
  expect_error(run_fcm(ovarian_FDA, model_name = "fit_lfcm",
                       formula = survival_time ~ age, event = "not_real_column",
                       metric = "uni k", r = "r", value = "fundiff",
                       analysis_vars = c("age", "survival_time"),
                       afcm = FALSE))
})

test_that("event needs to be 0/1", {
  expect_error(run_fcm(ovarian_FDA, model_name = "fit_lfcm",
                       formula = survival_time ~ age, event = "age",
                       metric = "uni k", r = "r", value = "fundiff",
                       analysis_vars = c("age", "survival_time"),
                       afcm = FALSE))
})
