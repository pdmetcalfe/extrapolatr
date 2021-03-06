#include <algorithm>
#include <exception>

#include <Rcpp.h>
#include <R_ext/Lapack.h>

namespace {

  Rcpp::NumericMatrix compute_delta(const Rcpp::NumericMatrix& mat) {
    Rcpp::NumericMatrix res(mat.nrow(), mat.ncol() - 1);
    for (R_len_t col=0; col!=mat.ncol() - 1; ++col) {
      res(Rcpp::_, col) = mat(Rcpp::_, col+1) - mat(Rcpp::_, col);
    }
    return res;
  }
}

//' Perform reduced-rank extrapolation
//'
//' This method removes any geometric terms from a sequence, leaving
//' you with the limit.  Note, as in the examples, that it will happily
//' remove divergent sequences.  This is a feature rather than a bug.
//' @param sequence a matrix, the columns of which are the sequence to
//'        accelerate
//' @return the RRE limit
//' @examples
//' ind <- seq_len(5)
//' vals <- 0.5^(ind - 1)
//' rre(matrix(vals, nrow=1))
//' rre(rbind(vals, 2 + vals))
//' @references
//' Eddy, R.P. (1979), Extrapolation to the limit of a vector
//'  sequence, in Wang, P.C.C., ed, Information Linkage Between
//'  Applied Mathematics and Industry, pp 387-396.
//'
//' Mesina, M. (1977), Convergence acceleration for the iterative solution
//'   of x=Ax + f, Methods Appl. Mech. Eng. 10, 165--173.
//'
//' Jbilou, K., and Sadok, H. (1995), Analysis of some vector extrapolation
//'   methods for solving systems of linear equations, Numer. Math.
//'   70, 73--89.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rre(const Rcpp::NumericMatrix& sequence) {

  if (sequence.ncol() <= 2) {
    throw(std::runtime_error("Need at least 3 columns"));
  }

  Rcpp::NumericMatrix delta(compute_delta(sequence));
  Rcpp::NumericMatrix ddelta(compute_delta(delta));

  const int m = sequence.nrow();
  const int n = sequence.ncol() - 2;
  const int lda = m;
  const int ldb = std::max(m, n);
  const int nrhs = 1;

  Rcpp::NumericVector b(ldb);
  // copy into rhs
  std::copy(delta.begin(), delta.begin() + m, b.begin());

  Rcpp::NumericVector s(std::min(m, n));

  double dummy_work;
  int dummy_iwork;
  int lwork = -1;
  int liwork = -1;
  int info;
  int rank;
  const double rcond=-1; // use machine precision
   // now do a workspace query
  F77_CALL(dgelsd)(const_cast<int*>(&m),
		   const_cast<int*>(&n),
		   const_cast<int*>(&nrhs),
		   ddelta.begin(), const_cast<int*>(&lda),
		   b.begin(), const_cast<int*>(&ldb),
		   s.begin(), const_cast<double*>(&rcond), &rank,
		   &dummy_work, &lwork,
		   &dummy_iwork,  &info);

  Rcpp::NumericVector work(1 + (int)dummy_work);
  Rcpp::IntegerVector iwork(dummy_iwork);
  lwork = work.length();
  F77_CALL(dgelsd)(const_cast<int*>(&m),
		   const_cast<int*>(&n),
		   const_cast<int*>(&nrhs),
		   ddelta.begin(), const_cast<int*>(&lda),
		   b.begin(), const_cast<int*>(&ldb),
		   s.begin(), const_cast<double*>(&rcond), &rank,
		   work.begin(), &lwork,
		   iwork.begin(), &info);

  if (info) {
    throw(std::runtime_error("dgelsd failed"));
  }

  // now prepare for output
  Rcpp::NumericVector res(sequence.nrow());
  std::copy(sequence.begin(), sequence.begin() + sequence.nrow(),
	    res.begin());
  const double alpha = -1;
  const double beta = 1;
  const int incx = 1;
  const int incy = 1;
  F77_CALL(dgemv)("N", &m, &n, &alpha,
		  delta.begin(), &m,
		  b.begin(), &incx,
		  &beta,
		  res.begin(), &incy);
  return res;
}

//' Perform minimal polynomial extrapolation
//'
//' This method removes any geometric terms from a sequence, leaving
//' you with the limit.  Note, as in the examples, that it will happily
//' remove divergent sequences.  This is a feature rather than a bug.
//' @param sequence a matrix, the columns of which are the sequence to
//'        accelerate
//' @return the MPE limit
//' @examples
//' ind <- seq_len(5)
//' vals <- 0.5^(ind - 1)
//' mpe(matrix(vals, nrow=1))
//' mpe(rbind(vals, 2 + vals))
//' @references
//' Cabay, S. and Jackson, L.W. (1976), A polynomial extrapolation method
//'    for finding limits and antilimits for vector sequences, SIAM J.
//'    Numer. Anal. 13, 734--752.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mpe(const Rcpp::NumericMatrix& sequence) {

  if (sequence.ncol() <= 2) {
    throw(std::runtime_error("Need at least 3 columns"));
  }

  Rcpp::NumericMatrix delta(compute_delta(sequence));

  const int m = sequence.nrow();
  const int n = delta.ncol() - 1;
  const int lda = m;
  const int ldb = std::max(m, n);
  const int nrhs = 1;

  // ensure that b is always big enough to add on an element.
  Rcpp::NumericVector b(std::max(ldb, sequence.ncol() - 1));
  // copy into rhs
  std::copy(delta.end() - m, delta.end(), b.begin());

  Rcpp::NumericVector s(std::min(m, n));

  double dummy_work;
  int dummy_iwork;
  int lwork = -1;
  int info;
  int rank;
  const double rcond=-1;
  // now do a workspace query
  F77_CALL(dgelsd)(const_cast<int*>(&m),
		   const_cast<int*>(&n),
		   const_cast<int*>(&nrhs),
		   delta.begin(), const_cast<int*>(&lda),
		   b.begin(), const_cast<int*>(&ldb),
		   s.begin(), const_cast<double*>(&rcond), &rank,
		   &dummy_work, &lwork,
		   &dummy_iwork, &info);


  Rcpp::NumericVector work(1 + (int)dummy_work);
  Rcpp::IntegerVector iwork(dummy_iwork);
  lwork = work.length();
  F77_CALL(dgelsd)(const_cast<int*>(&m),
		   const_cast<int*>(&n),
		   const_cast<int*>(&nrhs),
		   delta.begin(), const_cast<int*>(&lda),
		   b.begin(), const_cast<int*>(&ldb),
		   s.begin(), const_cast<double*>(&rcond), &rank,
		   work.begin(), &lwork,
		   iwork.begin(), &info);

  if (info) {
    throw(std::runtime_error("dgelsd failed"));
  }

  // now prepare for output
  Rcpp::NumericVector res(sequence.nrow(), 0.0);
  b[n] = -1;
  const int nn = sequence.ncol() - 1;
  const double scale = 1/std::accumulate(b.begin(), b.begin() + nn, 0.0);
  const double alpha = scale;
  const double beta = 0;
  const int incx = 1;
  const int incy = 1;
  F77_CALL(dgemv)("N", &m, &nn, &alpha,
                  sequence.begin(), &m,
		  b.begin(), &incx,
		  &beta,
		  res.begin(), &incy);
  return res;
}
