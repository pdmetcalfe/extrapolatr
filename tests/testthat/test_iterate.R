context("testing iterate helper")


test_that("preserves fixed points", {
    expect_equal(iterate(0, function(x) {x}, 10, save=FALSE), 0)
})


test_that("preserves nontrivial fixed points", {
    expect_equal(iterate(1, function(x) { 0.5 * (x + 2/x) }, 10, save=FALSE),
                 sqrt(2))
})


test_that("gives expected answers", {
    trace <- iterate(1, function(x) {0.5 * x}, 5, save=TRUE)
    record <- attr(trace, 'history')
    expect_equal(as.vector(record), 0.5**seq_len(5))
    expect_equal(nrow(record), 1)
    expect_equal(ncol(record), 5)
})


test_that("correctly matches final value", {
    trace <- iterate(1, function(x) {runif(1)}, 5, save=TRUE)
    record <- attr(trace, 'history')
    expect_equal(record[, 5, drop=TRUE], trace, check.attributes=FALSE)
})
