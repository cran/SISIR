library("SISIR")

context("Test that graphics for `SFCB` objects work as expected...")

data("truffles")

test_that("graphics for 'importance' work as expected.", {
  
  out1 <- sfcb(rainfall, truffles, group.method = "adjclust",
               summary.method = "pls")
  p <- plot(out1, plot.type = "importance")
  expect_s3_class(p, "ggplot")
  p <- plot(out1, plot.type = "importance", shape.imp = "histogram")
  expect_s3_class(p, "ggplot")
  
  out2 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", selection.method = "relief", 
               range.at = c(5, 12))
  p <- plot(out2, plot.type = "importance")
  expect_s3_class(p, "ggplot")
  p <- plot(out2, plot.type = "importance", shape.imp = "histogram")
  expect_s3_class(p, "ggplot")
  
  out3 <- sfcb(rainfall, truffles, group.method = "cclustofvar", 
               summary.method = "basics")
  p <- plot(out3, plot.type = "importance")
  expect_s3_class(p, "ggplot")
  p <- plot(out3, plot.type = "importance", shape.imp = "histogram")
  expect_s3_class(p, "ggplot")
  
  out4 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "basics", range.at = c(5, 7))
  p <- plot(out4, plot.type = "importance")
  expect_s3_class(p, "ggplot")
  p <- plot(out4, plot.type = "importance", shape.imp = "histogram")
  expect_s3_class(p, "ggplot")
  
  out5 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "basics", selection.method = "relief", 
               range.at = c(5, 7))
  p <- plot(out5, plot.type = "importance")
  expect_s3_class(p, "ggplot")
  p <- plot(out5, plot.type = "importance", shape.imp = "histogram")
  expect_s3_class(p, "ggplot")
  
})

test_that("graphics for 'selection' work as expected.", {
  
  out1 <- sfcb(rainfall, truffles, group.method = "adjclust",
               summary.method = "pls")
  p <- plot(out1, plot.type = "selection")
  expect_s3_class(p, "ggplot")
  p <- plot(out1, plot.type = "selection", sel.type = "selection", 
            threshold = 0.00011)
  expect_s3_class(p, "ggplot")
  
  out2 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", selection.method = "relief")
  p <- plot(out2, plot.type = "selection")
  expect_s3_class(p, "ggplot")
  p <- plot(out2, plot.type = "selection", sel.type = "selection")
  expect_s3_class(p, "ggplot")
  
  out3 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", selection.method = "relief", 
               range.at = c(5, 12))
  p <- plot(out3, plot.type = "selection")
  expect_s3_class(p, "ggplot")
  p <- plot(out3, plot.type = "selection", sel.type = "selection")
  expect_s3_class(p, "ggplot")
  
  out4 <- sfcb(rainfall, truffles, group.method = "cclustofvar", 
               summary.method = "basics")
  p <- plot(out4, plot.type = "selection")
  expect_s3_class(p, "ggplot")
  p <- plot(out4, plot.type = "selection", sel.type = "selection", 
            threshold = 0.0015)
  expect_s3_class(p, "ggplot")
  
  out5 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "basics", range.at = c(5, 7))
  suppressWarnings({ p <- plot(out5, plot.type = "selection") })
  expect_s3_class(p, "ggplot")
  p <- plot(out5, plot.type = "selection", sel.type = "selection", 
            threshold = 0.005)
  expect_s3_class(p, "ggplot")
  
  out6 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "basics", selection.method = "relief")
  p <- plot(out6, plot.type = "selection", sel.type = "selection")
  expect_s3_class(p, "ggplot")
  
  out7 <- sfcb(rainfall, truffles, group.method = "adjclust", 
                summary.method = "basics", selection.method = "relief", 
                range.at = c(5, 7))
  p <- plot(out7, plot.type = "selection", sel.type = "selection")
  expect_s3_class(p, "ggplot")
  
  out8 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", range.at = c(5, 7))
  p <- plot(out8, plot.type = "selection", sel.type = "selection", 
            threshold = 0.02)
  expect_s3_class(p, "ggplot")
  
  expect_error({ p <- plot(out8, plot.type = "selection", 
                           sel.type = "selection", threshold = "1") },
               "'threshold' must be numeric.", fixed = FALSE)
  expect_error({ p <- plot(out7, plot.type = "selection", 
                           sel.type = "selection", threshold = 1) },
               "A selection method has already been used", fixed = FALSE)
  expect_warning({ p <- plot(out4, plot.type = "selection",
                             threshold = 0.0015) },
                 "'sel.type' is not 'selection'. Automatically switching", 
                 fixed = FALSE)
  expect_error({ p <- plot(out4, plot.type = "selection", 
                           sel.type = "selection") },
               "Choose sel.type = 'importance'", fixed = FALSE)
  
})

test_that("graphics for 'dendrogram' work as expected.", {
  out1 <- sfcb(rainfall, truffles, group.method = "adjclust",
               summary.method = "pls")
  expect_message({ plot(out1) }, "Reversals detected in the dendrogram.")
  
  out2 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", selection.method = "relief")
  expect_message({ plot(out2) }, "Reversals detected in the dendrogram.")
  
  out3 <- sfcb(rainfall, truffles, group.method = "adjclust", 
               summary.method = "pls", selection.method = "relief", 
               range.at = c(5, 12))
  expect_message({ plot(out3) }, "Only the first 3 selections are represented")
  
  out4 <- sfcb(rainfall, truffles, group.method = "cclustofvar", 
               summary.method = "basics")
  expect_message({ plot(out4) }, "Reversals detected in the dendrogram.")
})
