library(testthat)

source("lib/lib.R")
source("lib/geo_annotation.R")


test_that("Test extractTissue", {
  characteristics = "tissue: foo bar"
  expect_equal(extractTissue(characteristics), "foo bar")
  
  characteristics = "other: stuff;  tissue: foo bar   ;  more stuff "
  expect_equal(extractTissue(characteristics), "foo bar")
})