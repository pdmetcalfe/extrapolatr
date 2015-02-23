context("RRE")

test_that("rre gives good answers (1)", {
    ind <- seq_len(5) - 1
    xs <- 0.99^ind
    expect_equal(rre(matrix(2 + xs, nrow=1)), 2)
})

test_that("rre gives good answers (2)", {
    ind <- seq_len(5) - 1
    xs <- 0.99^ind
    ys <- 0.5^ind
    zs <- 0.25^ind
    expect_equal(rre(rbind(1 + xs,
                           2 + ys,
                           3 + zs)), c(1, 2, 3))
})

test_that("rre gives good answers (3)", {
    ind <- seq_len(5) - 1
    xs <- 0.99^ind
    ys <- (-0.5)^ind
    zs <- 0.25^ind
    expect_equal(rre(rbind(1 + xs,
                           2 + ys,
                           3 + zs)), c(1, 2, 3))
})

test_that("rre gives good answers (4)", {
    ind <- seq_len(5) - 1
    xs <- 0.99^ind
    ys <- (-2)^ind
    zs <- 0.25^ind
    expect_equal(rre(rbind(1 + xs + ys,
                           2 + ys + zs,
                           3 + zs)), c(1, 2, 3))
})
