library("SISIR")

context("Test that `sfcb` works as expected...")

data("truffles")

test_that("`sfcb` works for one `at` and no selection.", {
  expected_outputs <- c("dendro", "groups", "summaries", "mse", "importances",
                        "computational.times", "call")
  
  out1 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", seed = 3)
  expect_named(out1, expected_outputs)
  out1b <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", seed = 3)
  out1$"computational.times" <- NULL
  out1b$"computational.times" <- NULL
  expect_identical(out1, out1b)
  
  out2 <- sfcb(rainfall, truffles, group.method = "cclustofvar", 
               summary.method = "cclustofvar", keep.time = FALSE)
  expect_named(out2, setdiff(expected_outputs, "computational.times"))
  
  out3 <- sfcb(rainfall, truffles, group.method = "cclustofvar", 
               summary.method = "basics", repeats = 6)
  expect_named(out3, expected_outputs)
  expect_equal(nrow(out3$mse), 6)
  expect_equal(ncol(out3$importances[[1]]), 6 * 2)
})

test_that("`sfcb` works for `range_at` and no selection.", {
  expected_outputs <- c("dendro", "groups", "summaries", "mse", "importances",
                        "computational.times", "call")
  
  out4 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", range.at = c(5, 7))
  expect_named(out4, expected_outputs)
  expect_length(out4$summaries, 3)
  
  out5 <- sfcb(rainfall, truffles, group.method = "cclustofvar", 
               summary.method = "cclustofvar", range.at = c(5, 7))
  expect_named(out5, expected_outputs)
  
  out6 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "basics", range.at = c(5, 7))
  expect_named(out6, expected_outputs)
})

test_that("`sfcb` works for one `at` with selection.", {
  expected_outputs <- c("dendro", "groups", "summaries", "selected", "mse", 
                        "importances", "computational.times", "call")
  
  out7 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", selection.method = "relief")
  expect_named(out7, expected_outputs)
  
  out8 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", selection.method = "relief")
  expect_named(out8, expected_outputs)
  
  out9 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "basics", selection.method = "relief")
  expect_named(out9, expected_outputs)
  
  out10 <- sfcb(rainfall, truffles, group.method = "adjclust", 
                summary.method = "pls", selection.method = "boruta")
  expect_named(out10, expected_outputs)
})

test_that("`sfcb` works for one `at` with selection.", {
  expected_outputs <- c("dendro", "groups", "summaries", "selected", "mse", 
                        "importances", "computational.times", "call")
  
  out11 <- sfcb(rainfall, truffles, group.method = "adjclust", 
                summary.method = "pls", selection.method = "relief", 
                range.at = c(5, 7))
  expect_named(out11, expected_outputs)
  expect_length(out11$selected, 3)
  
  out12 <- sfcb(rainfall, truffles, group.method = "adjclust", 
                summary.method = "pls", selection.method = "boruta", seed = 3)
  expect_named(out12, expected_outputs)
  
  out13 <- sfcb(rainfall, truffles, group.method = "adjclust", 
                summary.method = "basics", selection.method = "relief", 
                range.at = c(5, 7))
  expect_named(out13, expected_outputs)
  
  out14 <- sfcb(rainfall, truffles, group.method = "adjclust", 
                summary.method = "pls", selection.method = "relief", 
                range.at = c(5, 12))
  expect_named(out14, expected_outputs)
})

test_that("`sfcb` properly returns error when expected.", {
  expect_error({sfcb(rainfall, truffles, at = "1")},
               "'at' must be a positive integer", fixed = FALSE)
  expect_error({sfcb(rainfall, truffles, at = -3)},
               "'at' must be a positive integer", fixed = FALSE)
  expect_error({sfcb(rainfall, truffles, at = 2.4)},
               "'at' must be a positive integer", fixed = FALSE)
  
  expect_error({sfcb(rainfall, truffles, range.at = 1:3)},
               "'range.at' must be a vector of minimum and maximum numbers", 
               fixed = FALSE)
  expect_error({sfcb(rainfall, truffles, range.at = c("1", "2"))},
               "'range.at' must be a vector of minimum and maximum numbers", 
               fixed = FALSE)
  expect_error({sfcb(rainfall, truffles, range.at = c(-1, 2))},
               "'range.at' must be a vector of minimum and maximum numbers", 
               fixed = FALSE)
  expect_error({sfcb(rainfall, truffles, range.at = c(1.4, 2))},
               "'range.at' must be a vector of minimum and maximum numbers", 
               fixed = FALSE)
  expect_error({sfcb(rainfall, truffles, range.at = c(5, 2))},
               "'range.at' must be a vector of minimum and maximum numbers", 
               fixed = FALSE)
  
  expect_error({sfcb(rainfall, truffles, seed = "1")}, 
               "'seed' must be numeric!", fixed = FALSE)
})
