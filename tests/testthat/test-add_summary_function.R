data("ovarian_FDA")

new_FDA = make_mxfda(metadata = ovarian_FDA@Metadata,
                     spatial = NULL,
                     subject_key = ovarian_FDA@subject_key,
                     sample_key = ovarian_FDA@sample_key)
summ_fun = ovarian_FDA@univariate_summaries$Gest

test_that("function works", {
  expect_true(
    is(add_summary_function(new_FDA,
                            summ_fun,
                            metric = "uni g"), "mxFDA")
  )
})

test_that("new data is data frame",{
  expect_error(add_summary_funciton(new_FDA,
                                    list(summ_fun),
                                    metric = "uni g"))
})

test_that("new data is data frame",{
  expect_error(add_summary_funciton(new_FDA,
                                    summ_fun,
                                    metric = c("uni", "g")))
})
