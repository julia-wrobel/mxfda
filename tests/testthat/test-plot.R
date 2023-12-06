#prep data
data("ovarian_FDA")

ovarian_FDA = run_fpca(ovarian_FDA, metric = "uni g", r = "r", value = "fundiff",
                       lightweight = TRUE,
                       pve = .99)

test_that('plotting summary function', {
  expect_s3_class(plot(ovarian_FDA, y = 'fundiff', what = 'uni g'), "gg")
})

test_that('plotting fpca', {
  expect_s3_class(plot(ovarian_FDA, what = 'uni g fpca', pc_choice = 1), "gg")
})

test_that('plotting fpca', {
  #sneaky way for ggplot to do error handling
  expect_error(plot(ovarian_FDA, what = 'uni g', pc_choice = 1))
})
