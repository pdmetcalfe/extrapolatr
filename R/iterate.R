## Helper code for point iteration

##' Perform a point iteration
##'
##' Run a point iteration, saving intermediate values.  This is
##' a convenient building block for convergence acceleration.
##' @param x the start value of the iteration (not saved)
##' @param f the function that performs the iteration
##' @param n the number of iterations to do
##' @param save do we save the iteration history?
##' @param ... extra arguments for f
##' @return the output of the iteration.  If save is TRUE, it has an
##' attribute 'history', which gives the iteration history as a
##' length(x) x n matrix.
##' @examples
##' iterate(0.2, cos, 100, save=FALSE)
##' iterate(c(-1, 1), function(x) { - x }, 100)
##' @export
iterate <- function(x, f, n, save=TRUE, ...) {
    x <- as.vector(x)
    if (save) {
        result <- matrix(0, nrow=length(x), ncol=n)
        for (ind in seq_len(n)) {
            result[,ind] <- f(x, ...)
            x <- result[, ind, drop=TRUE]
        }
        attr(x, 'history') <- result
    } else {
        for (ind in seq_len(n)) {
            x <- f(x, ...)
        }
    }
    x
}

##' Accelerate a point iteration
##'
##' Use the techniques from this package to accelerate a point iteration.
##' @param start the start value
##' @param f the function to iterate
##' @param method the acceleration function.  must take an m x n
##' matrix and return an m-vector.
##' @param n.prior the number of prior iterations to run before accelerating
##' @param pre.iter how many iterations to do before acceleration (per cycle)
##' @param num.accel how many iterations to accelerate
##' @param num.cycles how many cycles to do?
##' @param ... additional arguments for f
##' @examples
##' fastIterate(c(0.2, 0.6), function(x) {c(1 - sin(x[[1]]), cos(x[[2]]))},
##'             rre, 1, 3, 5, 10)
##' 
##' @export
fastIterate <- function(start, f, method, n.prior, pre.iter,
                        num.accel, num.cycles, ...) {
    ## first do some calming iterations
    start <- iterate(start, f, n.prior, save=FALSE, ...)
    ## now cycle
    iterate(start,
            function(val, method, pre.iter, num.accel, ...) {
                val <- iterate(val, f, pre.iter, save=FALSE, ...)
                val <- iterate(val, f, num.accel, save=TRUE, ...)
                method(attr(val, 'history'))
            },
            num.cycles, save=FALSE, method, pre.iter, num.accel, ...)
}
