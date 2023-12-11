data("ovarian_FDA")

ovarian_FDA = run_fcm(ovarian_FDA, model_name = "fit_lfcm",
                      formula = survival_time ~ age, event = "event",
                      metric = "uni g", r = "r", value = "fundiff",
                      analysis_vars = c("age", "survival_time"),
                      afcm = FALSE)

test_that("output class is correct", {
  tmp = extract_surface(ovarian_FDA, metric = "uni g", model = "fit_lfcm", analysis_vars = c("age"))
  expect_true(inherits(tmp, 'lfcmSurface'))
})



test_that("value needs to be in summary", {
  expect_error(
    extract_surface(ovarian_FDA, metric = "uni g",
                    r = "r", value = "Degree of Clustering Exact",
                    model = "fit_lfcm", analysis_vars = c("age"))
  )
})

test_that("r needs to be in summary", {
  expect_error(
    extract_surface(ovarian_FDA, metric = "uni g",
                    r = "radius", value = "fundiff",
                    model = "fit_lfcm", analysis_vars = c("age"))
  )
})

test_that("metric needs to be length 1", {
  expect_error(
    extract_surface(ovarian_FDA, metric = c("uni g", 'uni k'),
                    r = "radius", value = "fundiff",
                    model = "fit_lfcm", analysis_vars = c("age"))
  )
})
