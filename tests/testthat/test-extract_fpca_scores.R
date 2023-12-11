#data preparation
data("ovarian_FDA")

ovarian_FDA = run_fpca(ovarian_FDA, metric = "uni g", r = "r", value = "fundiff",
                       lightweight = TRUE,
                       pve = .99)

#1-to-1 fpca, not mfpca
test_that("class is right", {
  fpc = extract_fpca_scores(ovarian_FDA, 'uni g fpca')
  expect_true(inherits(fpc, "data.frame")) #so much more clean!
})

#for mfpca, should just be a list of scores vs metadata/score table
#this is due to the different lengths of level1 and level2 scores in mfpca

test_that("'what' is of length 1", {
  expect_error(extract_fpca_scores(ovarian_FDA, c('uni g fpca', 'uni k fpca')))
})
