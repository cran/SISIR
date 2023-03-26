library("SISIR")

context("Test that extract works as expected...")

data("truffles")

test_that("extract works as expected without selection.", {
  expected_outputs <- c("dendro", "groups", "summaries", "mse", "importances",
                        "call")
  
  out1 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", range.at = c(5, 7))
  out1b <- extract_at(out1, 5)
  expect_named(out1b, expected_outputs)
  expect_length(out1b$groups, 1)
  expect_length(out1b$summaries, 1)
  expect_length(out1b$importances, 1)
  expect_equal(nrow(out1b$mse), 5)
  
  out2 <- sfcb(rainfall, truffles, group.method = "cclustofvar", 
               summary.method = "cclustofvar", range.at = c(5, 7))
  out2b <- extract_at(out2, 5)
  expect_named(out2b, expected_outputs)
  
  out3 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "basics", range.at = c(5, 7))
  out3b <- extract_at(out3, 5:6)
  expect_named(out3b, expected_outputs)
  expect_length(out3b$groups, 2)
  expect_length(out3b$summaries, 2)
  expect_length(out3b$importances, 2)
  expect_equal(nrow(out3b$mse), 5 * 2)
})

test_that("extract works as expected with selection.", {
  expected_outputs <- c("dendro", "groups", "summaries", "selected", "mse", 
                        "importances", "call")
  
  out4 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", selection.method = "relief", 
               range.at = c(5, 7))
  out4b <- extract_at(out4, 6)
  expect_named(out4b, expected_outputs)
  
  out5 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", selection.method = "relief", 
               range.at = c(5, 12))
  out5b <- extract_at(out5, c(9, 11:12))
  expect_named(out5b, expected_outputs)
  expect_length(out5b$groups, 3)
  expect_length(out5b$summaries, 3)
  expect_length(out5b$importances, 3)
  expect_equal(nrow(out5b$mse), 5 * 3)
})

test_that("extract returns errors as expected.", {
  out1 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", range.at = c(5, 7))
  expect_error({ extract_at(out1, "A") }, "'at' must be a numeric vector",
               fixed = FALSE)
  expect_error({ extract_at(out1, 9) }, 
               "'at' must be included in the range of tested groups for",
               fixed = FALSE)
})
