library(PsigA)
library(GSEABase)
library(breastCancerVDX)
context("Peak Distance")

data(vdx)

data <- exprs(vdx)[1:10, 1:100]
signature <- rownames(data)

distance <- peakDistance2d(signature, data, threshold = 0.005, n = 200,
                           magnitude = FALSE, scale = FALSE)

test_that("peakDistance2d output class is numeric", {
    expect_equal(class(distance), "numeric")
})

test_that("peakDistance2d output type is double", {
    expect_equal(typeof(distance), "double")
})

test_that("peakDistance2d score is consistent", {
    expect_equal(distance, c(1.465414, 10, 10),tolerance=1e-6)
})
