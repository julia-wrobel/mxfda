data("lung_FDA")

lung_FDA <- run_mfpca(lung_FDA,
                         metric = "uni g",
                         r = "r",
                         value = "fundiff",
                         pve = .99)

test_that("mxFDA object class is right", {
  expect_true(inherits(lung_FDA, "mxFDA"))
})
