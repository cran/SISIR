library("SISIR")

context("Test that quality computation and plots for `SFCB` objects work as expected...")

data("truffles")
beta <- c(0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0)

test_that("quality computation works as expected without selection.", {
  expected_outputs <- c("dendro", "groups", "summaries", "mse", "importances",
                        "computational.times", "call", "truth",
                        "quality", "threshold")
  
  out1 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls")
  expect_named(quality(out1, beta, threshold = 0.001), expected_outputs)
  
  out2 <- sfcb(rainfall, truffles, group.method = "cclustofvar", 
               summary.method = "basics")
  expect_named(quality(out2, beta, threshold = 0.001), expected_outputs)
  
  out3 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", range.at = c(5, 7))
  expect_named(quality(out3, beta, threshold = 0.001), expected_outputs)
  
  out4 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "basics", range.at = c(5, 7))
  expect_named(quality(out4, beta, threshold = 0.01), expected_outputs)
})

test_that("quality computation works as expected with selection.", {
  expected_outputs <- c("dendro", "groups", "summaries", "selected", "mse", 
                        "importances", "computational.times", "call", "truth",
                        "quality")
  
  out4 <- sfcb(rainfall, truffles, group.method = "adjclust",
               summary.method = "pls", selection.method = "relief")
  expect_named(quality(out4, beta), expected_outputs)
  
  out5 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", selection.method = "relief", 
               range.at = c(5, 7))
  expect_named(quality(out5, beta), expected_outputs)
  
  out6 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "basics", selection.method = "relief")
  expect_named(quality(out6, beta), expected_outputs)
  
  out7 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "basics", selection.method = "relief", 
               range.at = c(5, 7))
  expect_named(quality(out7, beta), expected_outputs)
  
  out8 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", selection.method = "relief", 
               range.at = c(5, 12))
  expect_named(quality(out8, beta), expected_outputs)
  expect_named(quality(out8, beta, threshold = 0.01), 
               c(expected_outputs, "threshold"))
})

test_that("quality graphics works as expected.", {
  out1 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls")
  out1b <- quality(out1, beta, threshold = 0.001)
  p <- plot(out1b, plot.type = "quality", quality.crit = "mse")
  expect_s3_class(p, "ggplot")
  p <- plot(out1b, plot.type = "quality", quality.crit = "time")
  expect_s3_class(p, "ggplot")
  p <- plot(out1b, plot.type = "quality", quality.crit = "ARI")
  expect_s3_class(p, "ggplot")
  p <- plot(out1b, plot.type = "quality", quality.crit = "NMI")
  expect_s3_class(p, "ggplot")
  p <- plot(out1b, plot.type = "quality", quality.crit = c("mse", "NMI"))
  expect_s3_class(p, "ggplot")
  p <- plot(out1b, plot.type = "quality", 
            quality.crit = c("Precision", "Recall"))
  expect_s3_class(p, "ggplot")
  
  out2 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", range.at = c(5, 7))
  out2b <- quality(out2, beta, threshold = 0.001)
  p <- plot(out2b, plot.type = "quality", quality.crit = "mse")
  expect_s3_class(p, "ggplot")
  p <- plot(out2b, plot.type = "quality", quality.crit = "time")
  expect_s3_class(p, "ggplot")
  p <- plot(out2b, plot.type = "quality", quality.crit = "NMI")
  expect_s3_class(p, "ggplot")
  p <- plot(out2b, plot.type = "quality", quality.crit = c("mse", "ARI"))
  expect_s3_class(p, "ggplot")
  p <- plot(out2b, plot.type = "quality", 
            quality.crit = c("Precision", "Recall"))
  expect_s3_class(p, "ggplot")
  
  out3 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", selection.method = "relief", 
               range.at = c(5, 12))
  out3b <- extract_at(out3, c(9, 11:12))
  out3c <- quality(out3b, beta, threshold = 0.01)
  p <- plot(out3c, plot.type = "quality", quality.crit = "mse")
  expect_s3_class(p, "ggplot")
  p <- plot(out3c, plot.type = "quality", quality.crit = "ARI")
  expect_s3_class(p, "ggplot")
  p <- plot(out3c, plot.type = "quality", quality.crit = c("mse", "NMI"))
  expect_s3_class(p, "ggplot")
  p <- plot(out3c, plot.type = "quality", 
            quality.crit = c("Precision", "Recall"))
  expect_s3_class(p, "ggplot")
})

test_that("quality computation returns errors as expected.", {
  out1 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls")
  expect_error({ quality(out1, c(0, 1), threshold = 0.001) },
               "'ground_truth' must have a length identical to initial number",
               fixed = FALSE)
  expect_error({ quality(out1, beta, threshold = "A") },
               "'threshold' must be a positive number or NULL.",
               fixed = FALSE)
  expect_error({ quality(out1, beta, threshold = -3) },
               "'threshold' must be a positive number or NULL.",
               fixed = FALSE)
  expect_error({ quality(out1, beta) },
               "No selected interval in this 'SFCB' object and no 'threshold'",
               fixed = FALSE)
})

test_that("quality graphics return errors as expected.", {
  out1 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls")
  out1b <- quality(out1, beta, threshold = 0.001)
  expect_error({ plot(out1b, plot.type = "quality", 
                      quality.crit = c("mse", "time")) },
               "'time' is a valid quality criterion to plot only taken alone.",
               fixed = FALSE)
  expect_error({plot(out1, plot.type = "quality", quality.crit = "Precision")},
               "'quality.crit' must be a vector with length at most 2 in",
               fixed = FALSE)
  expect_error({plot(out1, plot.type = "quality", 
                     quality.crit = c("mse", "time", "Precision"))},
               "'quality.crit' must be a vector with length at most 2 in",
               fixed = FALSE)
  expect_error({ plot(out1, plot.type = "quality", quality.crit = "AA") },
               "'quality.crit' must be a vector with length at most 2 in",
               fixed = FALSE)
  
  out2 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", keep.time = FALSE)
  expect_error({ plot(out2, plot.type = "quality", quality.crit = "time") },
               "'quality.crit' must be a vector with length at most 2 in",
               fixed = FALSE)
})
