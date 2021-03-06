% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{mpe}
\alias{mpe}
\title{Perform minimal polynomial extrapolation}
\usage{
mpe(sequence)
}
\arguments{
\item{sequence}{a matrix, the columns of which are the sequence to
accelerate}
}
\value{
the MPE limit
}
\description{
This method removes any geometric terms from a sequence, leaving
you with the limit.  Note, as in the examples, that it will happily
remove divergent sequences.  This is a feature rather than a bug.
}
\examples{
ind <- seq_len(5)
vals <- 0.5^(ind - 1)
mpe(matrix(vals, nrow=1))
mpe(rbind(vals, 2 + vals))
}
\references{
Cabay, S. and Jackson, L.W. (1976), A polynomial extrapolation method
   for finding limits and antilimits for vector sequences, SIAM J.
   Numer. Anal. 13, 734--752.
}

