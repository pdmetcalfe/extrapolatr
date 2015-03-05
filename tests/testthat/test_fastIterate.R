context("testing fastIterate helper")


test_that("usefully accelerates #1", {
    tmp <- fastIterate(1, function(x) {1 + 0.99 * x}, rre,
                       5, 3, 5, 1)
    expect_equal(tmp, 100)
})


test_that("usefully accelerates #2", {
    tmp <- fastIterate(1, function(x) {1 - 0.99 * x}, mpe,
                       5, 3, 5, 1)
    expect_equal(tmp, 1/1.99)
})


test_that("correctly passes parameters", {
    my.scale <- 5
    tmp <- fastIterate(1, function(x, scale) {1 + scale * x}, rre,
                       0, 3, 5, 1, my.scale)
    expect_equal(tmp, 1/(1 - my.scale))
})
