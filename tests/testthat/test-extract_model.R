data("ovarian_FDA")

ovarian_FDA = run_fcm(ovarian_FDA, model_name = "fit_lfcm",
                      formula = survival_time ~ age, event = "event",
                      metric = "uni g", r = "r", value = "fundiff",
                      analysis_vars = c("age", "survival_time"),
                      afcm = FALSE)

test_that("model class is right", {
  mod = extract_model(ovarian_FDA, 'uni g', 'cox', 'fit_lfcm')
  expect_equal(class(mod)[1], 'lfcm')
})

test_that("need to provide a model", {
  expect_error(extract_model(ovarian_FDA, 'uni g', 'cox'))
})

#only g exists so error will be about length of 'what' vs missing model
test_that("only one model extraction at a time", {
  expect_error(extract_model(ovarian_FDA, c('uni g', 'uni k'), 'fit_lfcm'))
})
