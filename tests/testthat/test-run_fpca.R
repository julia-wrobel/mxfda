#data preparation
data("ovarian_FDA")

#univariate fpca tests
test_that("mxFDA object class is right", {
  ovarian_FDA = run_fpca(ovarian_FDA, metric = "uni g", r = "r", value = "fundiff",
                         lightweight = TRUE,
                         pve = .99)
  expect_true(inherits(ovarian_FDA, "mxFDA"))
})

#just looks at the length of 'metric', not wther the summary of metric exists
test_that("metric needs to be length 1", {
  expect_error(run_fpca(ovarian_FDA, metric = c("uni g", "uni k"), r = "r", value = "fundiff",
                        lightweight = TRUE,
                        pve = .99))
})

test_that("metric summary must exist", {
  expect_error(
    run_fpca(ovarian_FDA, metric = "uni k", r = "r", value = "fundiff",
             lightweight = TRUE,
             pve = .99)
  )
})

test_that("value needs to be in summary", {
  expect_error(
    run_fpca(ovarian_FDA, metric = "uni g", r = "r", value = "Degree of Clustering Exact",
             lightweight = TRUE,
             pve = .99)
  )
})

test_that("r needs to be in summary", {
  expect_error(
    run_fpca(ovarian_FDA, metric = "uni g", r = "radius", value = "fundiff",
             lightweight = TRUE,
             pve = .99)
  )
})
