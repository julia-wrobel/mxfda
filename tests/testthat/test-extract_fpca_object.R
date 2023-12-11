#prep data
data("ovarian_FDA")

ovarian_FDA <- run_fpca(ovarian_FDA,
                        metric = "uni g",
                        r = "r",
                        value = "fundiff",
                        pve = .95)

test_that("length of what", {
  expect_error(extract_fpca_object(ovarian_FDA,
                                   'uni g'))
})

test_that("length of what", {
  out = extract_fpca_object(ovarian_FDA,
                      'uni g fpca')
  expect_true(is(out, "fpca"))
})
