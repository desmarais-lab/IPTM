#include <RcppArmadillo.h>
#include <cmath>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
// [[Rcpp::depends(RcppArmadillo)]]

using std::log;
using std::exp;
using std::max;
using std::min;
using std::abs;
using std::sqrt;
using std::pow;
using namespace Rcpp; 

void R_init_markovchain(DllInfo* info) {
	R_registerRoutines(info, NULL, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);	
}

extern void dsyev_( char *jobz, char *uplo, int *n, double *a, int *lda,
                   double *w, double *work, int *lwork, int *info );

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov){
    int n = x.n_rows;
    arma::mat x_cen;
    x_cen.copy_size(x);
    for (int i=0; i < n; i++) {
        x_cen.row(i) = x.row(i) - center;
    }
    return sum((x_cen * cov.i()) % x_cen, 1);
}

// [[Rcpp::export]]
arma::vec ei(arma::mat M) {
    return arma::eig_sym(M);
}

// [[Rcpp::export]]
arma::vec dmvnorm_arma (arma::mat x,  arma::rowvec mean,  arma::mat sigma){
    arma::vec distval = Mahalanobis(x, mean, sigma);
    arma::vec ei_sigma = ei(sigma);
    double logdet = sum(arma::log(ei_sigma));
    arma::vec logretval = -((x.n_cols * log2pi + logdet + distval)/2) ;
    return(logretval);
}

// [[Rcpp::export]]
arma::mat rmvnorm_arma(int n,
                  const arma::vec& mu,
                  const arma::mat& Sigma) {
    unsigned int p = Sigma.n_cols;
    Rcpp::NumericVector draw = Rcpp::rnorm(n*p);
    arma::mat Z = arma::mat(draw.begin(), n, p, false, true);
    arma::mat Y = arma::repmat(mu, 1, n).t()+Z * arma::chol(Sigma);
    return Y;
}

// **********************************************************//
//                  Transpose                                //
// **********************************************************//
// [[Rcpp::export]]
arma::mat transpose (arma::mat x) {
	return x.t();
}

// **********************************************************//
//                 C++ version of tabulate() in R            //
// **********************************************************//
// [[Rcpp::export]]
IntegerVector tabulateC(const IntegerVector& x, const signed max) {
    IntegerVector counts(max);
    std::size_t n = x.size();
    for (std::size_t i=0; i < n; i++) {
        if (x[i] > 0 && x[i] <= max)
            counts[x[i]-1]++;
    }
    return counts;
}

// **********************************************************//
//              Likelihood evaluation of Timepart            //
// **********************************************************//
// [[Rcpp::export]]
double Timepartsum (NumericMatrix mumat, double sigma_tau,
                   IntegerVector senders, NumericVector timestamps){
    int D = senders.size();
    double timesum = 0;
    for (unsigned int d = 0; d < D; d++) {
        int a_d = senders[d] - 1;
        timesum += R::dlnorm(timestamps[d], mumat(d, a_d), sigma_tau, TRUE);
        for (unsigned int i = 0; i < mumat.ncol(); i++) {
        if (i != a_d) {
            timesum += R::plnorm(timestamps[d], mumat(d, i), sigma_tau, FALSE, TRUE);
    	  }
        }
    }
    return timesum;
}

// **********************************************************//
//           Likelihood evaluation of Timepart--exp          //
// **********************************************************//
// [[Rcpp::export]]
double Timepartsum2 (NumericMatrix mumat, IntegerVector senders, NumericVector timestamps){
    int D = senders.size();
    double timesum = 0;
    for (unsigned int d = 0; d < D; d++) {
        int a_d = senders[d] - 1;
        timesum += R::dexp(timestamps[d], exp(mumat(d, a_d)), TRUE);
        for (unsigned int i = 0; i < mumat.ncol(); i++) {
        if (i != a_d) {
            timesum += R::pexp(timestamps[d], exp(mumat(d, i)), FALSE, TRUE);
    	  }
        }
    }
    return timesum;
}

// **********************************************************//
//              Likelihood evaluation of Timepart            //
// **********************************************************//
// [[Rcpp::export]]
NumericVector Timepartindiv (NumericVector mu, double sigma_tau, double timestamp){
    int A = mu.size();
    NumericVector out(A);
    for (unsigned int a = 0; a < A; a++) {
        out[a] += R::dlnorm(timestamp, mu[a], sigma_tau, TRUE);
        for (unsigned int i = 0; i < A; i++) {
            if (i != a) {
                out[a] += R::plnorm(timestamp, mu[i], sigma_tau, FALSE, TRUE);
            }
        }
    }
    return exp(out);
}

// **********************************************************//
//              Likelihood evaluation of Timepart            //
// **********************************************************//
// [[Rcpp::export]]
NumericVector Timepartindiv2 (NumericVector mu, double timestamp){
    int A = mu.size();
    NumericVector out(A);
    for (unsigned int a = 0; a < A; a++) {
        out[a] += R::dexp(timestamp, exp(mu[a]), TRUE);
        for (unsigned int i = 0; i < A; i++) {
            if (i != a) {
                out[a] += R::pexp(timestamp, exp(mu[i]), FALSE, TRUE);
            }
        }
    }
    return exp(out);
}

// **********************************************************//
//                    Multinomial Sampler                    //
// **********************************************************//
// [[Rcpp::export]]
int multinom_vec (NumericVector x) {
	NumericVector x1 = x / sum(x);
	int n = x1.size();
    IntegerVector d(n);
    R::rmultinom(1, x1.begin(), n, d.begin());
	for (unsigned int j = 0; j < n; j++) {
		if (d[j] == 1) {
			return j+1;
		}
	}
	return 0;
}

// **********************************************************//
//               Normalizing constant of Gibbs               //
// **********************************************************//
// [[Rcpp::export]]
double normalizer (arma::vec lambda_da) {
    double product = arma::prod(lambda_da+1)-1;
    if (product == 0) {
        product = exp(-745);
    }
	return log(product);
}

// **********************************************************//
//                   Posterior for Edgepart                  //
// **********************************************************//
// [[Rcpp::export]]
double Edgepartsum (List lambda, List u) {
	double out = 0;
	int D = lambda.size();
	IntegerMatrix samp = u[0];
	int A = samp.nrow();
	for (unsigned int d = 0; d < D; d++) {
		NumericMatrix lambda_d = lambda[d];
		NumericMatrix u_d = u[d];
		for (unsigned int a = 0; a < A; a++) {
			double num = sum(lambda_d(a, _) * u_d(a, _));
			NumericVector lambda_da = exp(lambda_d(a, _));
			lambda_da[a] = 0;
			double denom = normalizer(lambda_da);
			out += num - denom;
		}
	}
	return out;
}

// **********************************************************//
//                     Calculate mu matrix                   //
// **********************************************************//
// [[Rcpp::export]]
arma::mat mu_cpp (arma::cube Y, arma::rowvec eta) {
	int D = Y.n_rows;
	int A = Y.n_cols;
	int Q = Y.n_slices;
	arma::mat mu = arma::zeros(D, A);
	for (unsigned int d = 0; d < D; d++) {
		for (unsigned int a = 0; a < A; a++) {
			arma::vec Y_da = Y.subcube(d, a, 0, d, a, Q-1);
			mu(d, a) = sum(eta * Y_da);
		}
	}
	return mu;
}

// **********************************************************//
//                     Calculate mu matrix                   //
// **********************************************************//
// [[Rcpp::export]]
NumericVector mu_cpp_d (NumericMatrix Y, NumericVector eta) {
  int A = Y.nrow();
  NumericVector mu(A);
    for (unsigned int a = 0; a < A; a++) {
      NumericVector Y_a = Y(a, _);
      mu[a] = sum(eta * Y_a);
    }
  return mu;
}

// **********************************************************//
//                     Calculate lambda list                 //
// **********************************************************//
// [[Rcpp::export]]
arma::mat lambda_cpp (arma::cube X_d, arma::rowvec beta) {
	int A = X_d.n_rows;
	int P = X_d.n_slices;
	arma::mat lambda_d = arma::zeros(A, A);
	for (unsigned int a = 0; a < A; a++) {
		for (unsigned int r = 0; r < A; r++) {
			if (r != a) {
				arma::vec X_dar = X_d.subcube(a, r, 0, a, r, P-1);
				lambda_d(a,r) = sum(beta * X_dar);
			}
		}
	}
	return lambda_d;
}

// [[Rcpp::export]]
Rcpp::IntegerVector sample_int(int n, int min, int max) {
    Rcpp::IntegerVector pool = Rcpp::seq(min, max);
    std::random_shuffle(pool.begin(), pool.end());
    return pool[Rcpp::Range(0, n - 1)];
}

// **********************************************************//
//                       Gibbs update on u                   //
// **********************************************************//
// [[Rcpp::export]]
List u_cpp (List lambda, List u) {
	int D = lambda.size();
	IntegerMatrix samp = u[0];
	int A = samp.nrow();
	List out(D);
	NumericVector prob(2);
	for (unsigned int d = 0; d < D; d++) {
		IntegerMatrix u_d = u[d];
		NumericMatrix lambda_d = lambda[d];
		IntegerVector Asample = sample_int(A, 0, A-1);
		IntegerVector Rsample = sample_int(A, 0, A-1);
		for (unsigned int i = 0; i < A; i++) {
			int a = Asample[i];
			for (unsigned int j = 0; j < A; j++) {
				int r = Rsample[j];
				if (r != a) {
					double prob0 = sum(u_d(a,_)) - u_d(a,r);
					prob[0] = (prob0 > 0);
					prob[1] = exp(lambda_d(a, r));
					u_d(a, r) = multinom_vec(prob)-1;
				}
			}
		}
		out[d] = u_d;
	}
	return out;
}

// **********************************************************//
//                       Gibbs update on u                   //
// **********************************************************//
// [[Rcpp::export]]
IntegerMatrix u_cpp_d (NumericMatrix lambda_d, IntegerMatrix u) {
  int A = u.nrow();
  NumericVector prob(2);
  IntegerMatrix u_d = u;
  IntegerVector Asample = sample_int(A, 0, A-1);
  IntegerVector Rsample = sample_int(A, 0, A-1);
    for (unsigned int i = 0; i < A; i++) {
      int a = Asample[i];
      for (unsigned int j = 0; j < A; j++) {
        int r = Rsample[j];
        if (r != a) {
          double prob0 = sum(u_d(a,_)) - u_d(a,r);
          prob[0] = (prob0 > 0);
          prob[1] = exp(lambda_d(a, r));
          u_d(a, r) = multinom_vec(prob)-1;
        }
      }
    }
  return u_d;
}
