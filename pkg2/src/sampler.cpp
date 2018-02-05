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
NumericVector sortuniq (NumericVector x){
    return sort_unique(x);
}

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
//               Prior calculations for IP sum               //
// **********************************************************//
// [[Rcpp::export]]
double priorsum (arma::mat var, arma::rowvec mu, arma::mat x) {
	double priorsum = sum(dmvnorm_arma(x, mu, var));
	return priorsum;
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
//                  Transpose                                //
// **********************************************************//
// [[Rcpp::export]]
arma::mat transpose (arma::mat x) {
	return x.t();
}

// **********************************************************//
//                    Call rmultinom from R                  //
// **********************************************************//
// [[Rcpp::export]]
IntegerVector callRMultinom (NumericVector x) {
    NumericVector x1 = x / sum(x);
	int n = x1.size();
	IntegerVector d(n);
	R::rmultinom(1, x1.begin(), n, d.begin());
	return d;
}


// **********************************************************//
//                    Multinomial Sampler                    //
// **********************************************************//
// [[Rcpp::export]]
IntegerVector multinom_vec (int nSample, NumericVector props) {
	IntegerVector multinom_vec(nSample);
	for (int i = 0; i < nSample; i++) {
		IntegerVector multinom_i = callRMultinom(props);
		for (unsigned int j = 0; j < props.size(); j++) {
			if (multinom_i[j] == 1) {
				multinom_vec[i] = j+1;
			}
		}
	}
	return multinom_vec;
}

// **********************************************************//
//         Function to search for burn-in documents          //
// **********************************************************//
// [[Rcpp::export]]
int which_int(int value, IntegerVector x) {
	int n = x.size();
	for (int i = 0; i < n; i++) {
		if (x[i] >= value) {
			return i+1;
		}
	}
	return 1;
}

// **********************************************************//
//       Function to search for cutoff point in History      //
// **********************************************************//
// [[Rcpp::export]]
int which_num(double value, NumericVector x) {
  int n = x.size();
  for (int i = 0; i < n; i++) {
    if (x[i] >= value) {
      return i+1;
    }
  }
  return 0;
}

// **********************************************************//
//        Multiple Draws from a Dirichlet Distribution       //
// **********************************************************//
// [[Rcpp::export]]
NumericMatrix rdirichlet_cpp(int num_samples, NumericVector alpha_m) {
	int dist_size = alpha_m.size();
	NumericMatrix distribution(num_samples, dist_size);
	
	for (int i = 0; i < num_samples; ++i) {
		double sum_term = 0;
		for (int j = 0; j < dist_size; ++j) {
			double cur = R::rgamma(alpha_m[j], 1.0);
			distribution(i, j) = cur;
			sum_term += cur;
		}
		for (int j = 0; j < dist_size; ++j) {
			distribution(i, j) = distribution(i, j)/sum_term;
		}
	}
	return(distribution);
}

// **********************************************************//
//                            which in Rcpp                  //
// **********************************************************//
// [[Rcpp::export]]
IntegerVector which_cpp(int value, NumericVector x) {
  IntegerVector v = seq(1, x.size());
  return v[x==value];
}

// **********************************************************//
//          Calculate multinomial sampling probability       //
// **********************************************************//
// [[Rcpp::export]]
NumericVector expconst(NumericVector consts) {
    NumericVector expconsts = exp(consts-max(consts));
    LogicalVector zeros = (expconsts == 0);
    LogicalVector infs = (expconsts == arma::datum::inf);
    if (sum(zeros) > 0) {
        expconsts[zeros] = exp(-745);
    }
    if (sum(infs) > 0) {
        expconsts[infs] = exp(700);
    }
    return expconsts;
}

// **********************************************************//
//              Construct the history of interaction         //
// **********************************************************//
// [[Rcpp::export]]
List History(List edge, NumericVector timestamps, IntegerVector cd, int A, int d, double timeunit) {
    List histlist(3);
    for (unsigned int l = 0; l < 3; l++){
        IntegerMatrix histmat(A, A);
        histlist[l] = histmat;
    }
    double time1 = timestamps[d-2]-384*timeunit;
    double time2 = timestamps[d-2]-96*timeunit;
    double time3 = timestamps[d-2]-24*timeunit;
    int iter = which_num(time1, timestamps);
    
    int cdnow = cd[d-1];
    for (int i = iter-1; i < d-1; i++) {
        if (cd[i] == cdnow) {
            List document2 = edge[i];
            int sender = document2[0];
            IntegerVector receiver = document2[1];
            double time = timestamps[i];
            for (unsigned int r = 0; r < receiver.size(); r++){
                if (time < time2) {
                    IntegerMatrix hist_l = histlist[2];
                    hist_l(sender-1, receiver[r]-1) += 1;
                    histlist[2] = hist_l;
                } else {
                    if (time < time3) {
                        IntegerMatrix hist_l = histlist[1];
                        hist_l(sender-1, receiver[r]-1) += 1;
                        histlist[1] = hist_l;
                    } else {
                            IntegerMatrix hist_l = histlist[0];
                            hist_l(sender-1, receiver[r]-1) += 1;
                            histlist[0] = hist_l;
                    }
                }
            }
        }
    }
    return histlist;
}

// **********************************************************//
//              Construct the history of interaction         //
// **********************************************************//
// [[Rcpp::export]]
List History2(List edge, NumericVector timestamps, IntegerVector cd, int A, IntegerMatrix timeintd, double timeunit) {
    List histlist(3);
    for (unsigned int l = 0; l < 3; l++){
        IntegerMatrix histmat(A, A);
        histlist[l] = histmat;
    }
    int cdnow = cd[timeintd(0,1)];
    for (int i = timeintd(2,0)-1; i < timeintd(0,1); i++) {
        if (cd[i] == cdnow) {
        List document2 = edge[i];
        int sender = document2[0];
        IntegerVector receiver = document2[1];
        for (unsigned int r = 0; r < receiver.size(); r++){
            if (i < timeintd(2,1)) {
                IntegerMatrix hist_l = histlist[2];
                hist_l(sender-1, receiver[r]-1) += 1;
                histlist[2] = hist_l;
            } else {
                if (i < timeintd(1,1)) {
                    IntegerMatrix hist_l = histlist[1];
                    hist_l(sender-1, receiver[r]-1) += 1;
                    histlist[1] = hist_l;
                } else {
                    if (i >= timeintd(0,0)-1) {
                        IntegerMatrix hist_l = histlist[0];
                        hist_l(sender-1, receiver[r]-1) += 1;
                        histlist[0] = hist_l;
                    }
                }
            }
        }
      }
    }
    return histlist;
}

// **********************************************************//
//                  Calculate Degree statistics              //
// **********************************************************//
// [[Rcpp::export]]
NumericMatrix Degree(List history, int A, int sender) {
    NumericMatrix degreemat(A, 6);
    IntegerVector outdegree(A);
  for (unsigned int l = 0; l < 3; l++) {
        IntegerMatrix history_l = history[l];
        int outdegree = sum(history_l(sender, _));
        degreemat(_,l) = rep(outdegree, A); //currently not account for multicast
        degreemat(_,l+3) = colSums(history_l);
      }
  return degreemat / 10;
}

// **********************************************************//
//                  Calculate Outdegree statistics           //
// **********************************************************//
// [[Rcpp::export]]
NumericMatrix Outdegree(IntegerMatrix timeintd, IntegerVector cd,
                        IntegerVector senders, int A) {
    NumericMatrix degreemat(A, 3);
    int cdnow = cd[timeintd(0,1)];
    for (unsigned int l = 0; l < 3; l++) {
        IntegerVector senders_l = senders[Range(timeintd(l,0)-1, timeintd(l,1)-1)];
        IntegerVector cd_l = cd[Range(timeintd(l,0)-1, timeintd(l,1)-1)];
        degreemat(_,l) = tabulateC(senders_l[cd_l==cdnow], A);
    }
    return degreemat / 10;
}

// **********************************************************//
//                  Calculate Indegree statistics              //
// **********************************************************//
// [[Rcpp::export]]
NumericMatrix Indegree(List history, int A, int sender) {
    NumericMatrix degreemat(A, 3);
    for (unsigned int l = 0; l < 3; l++) {
        IntegerMatrix history_l = history[l];
        degreemat(_,l) = colSums(history_l);
    }
    return degreemat / 10;
}

// **********************************************************//
//                  Calculate Dyadic statistics              //
// **********************************************************//
// [[Rcpp::export]]
IntegerMatrix Dyadic(List history, int A, int sender) {
    IntegerMatrix dyadicmat(A, 6);
    IntegerVector dyadic(6);
        for (int receiver = 0; receiver < A; receiver++) {
            if (receiver != sender) {
            for (unsigned int l = 0; l < 3; l++) {
                IntegerMatrix history_l = history[l];
                dyadic[l] = history_l(sender, receiver);
                dyadic[l+3] = history_l(receiver, sender);
            }
            dyadicmat(receiver, _) = dyadic;
            }
        }
    return dyadicmat;
}

// **********************************************************//
//          Calculate (full 3 x 3) Triadic statistic         //
// **********************************************************//
// [[Rcpp::export]]
NumericMatrix Triadic(List history, int A, int sender) {
    IntegerMatrix triadmat(A, 36);
    IntegerVector triadic(36);
        
    for (int receiver = 0; receiver < A; receiver++) {
        if (receiver != sender) {
        IntegerVector twosend(A);
        IntegerVector tworeceive(A);
        IntegerVector sibling(A);
        IntegerVector cosibling(A);
        int iter = 0;
        for (unsigned int l = 0; l < 3; l++) {
            for (unsigned int m = 0; m < 3; m++){
                IntegerMatrix history_l = history[l];
                IntegerMatrix history_m = history[m];
                for (int third = 0; third < A; third++) {
                    if (third != sender && third != receiver) {
                        double stoh = history_l(sender, third);
                        double htos = history_l(third, sender);
                        double rtoh = history_m(receiver, third);
                        double htor = history_m(third, receiver);
                        twosend[third] = stoh*htor;
                        tworeceive[third] = htos*rtoh;
                        sibling[third] = htos*htor;
                        cosibling[third] = stoh*rtoh;
                    }
                }
                    triadic[iter] = sum(twosend);
                    triadic[iter+9] = sum(tworeceive);
                    triadic[iter+18] = sum(sibling);
                    triadic[iter+27] = sum(cosibling);
                    iter = iter+1;
                }
            }
            triadmat.row(receiver) = triadic;
        }
    }
        NumericMatrix triadmat2(A, 12);
        for (int t = 0; t < 4; t++) {
            int add = 9*t;
            for (int i = 0; i < A; i++) {
                triadmat2(i, 3*t+0) = triadmat(i, add+0);
                triadmat2(i, 3*t+1) = triadmat(i, add+1)+triadmat(i, add+3)+triadmat(i, add+4);
                triadmat2(i, 3*t+2) = triadmat(i, add+2)+triadmat(i, add+5)+triadmat(i, add+6)+
                triadmat(i, add+7)+triadmat(i, add+8);
            }                    
        }
    return triadmat2/10;
}

// **********************************************************//
//                    Network statistics                     //
// **********************************************************//
// [[Rcpp::export]]
List Netstats_cpp(List edge, NumericVector timestamps, IntegerMatrix timeintd, IntegerVector senders, IntegerVector cd, int A, int d, double timeunit, IntegerVector netstat) {
    int P = 3*(2*netstat[0]+2*netstat[1]+4*netstat[2]);
    List out(A);
    List history = History(edge, timestamps, cd, A, d, timeunit);
    NumericMatrix outdegree = Outdegree(timeintd, cd, senders, A);
    for (unsigned int a = 0; a < A; a++) {
        NumericMatrix netstatmat(A, P);
        int iter = 0;
        if (netstat[0] == 1) {
            NumericMatrix indegree = Indegree(history, A, a);
            for (int it = 0; it < 3; it++) {
                netstatmat(_, iter) = rep(outdegree(a,it), A);
                netstatmat(_, iter+3) = indegree(_,it);
                iter += 1;
            }
            iter = 6;
        }
        if (netstat[1] == 1) {
            IntegerMatrix dyadic = Dyadic(history, A, a);
            for (int it = 0; it < 6; it++) {
                netstatmat(_, iter) = dyadic(_,it);
                iter += 1;
            }
        }
        if (netstat[2] == 1) {
            NumericMatrix triadic = Triadic(history, A, a);
            for (int it = 0; it < 12; it++) {
                netstatmat(_, iter) = triadic(_,it);
                iter += 1;
            }
        }
        out[a] = netstatmat;
    }
    return out;
}

// **********************************************************//
//                  inner product between two vectors        //
// **********************************************************//
// [[Rcpp::export]]
double inner (arma::vec x, arma::vec y) {
	arma::mat ip = x.t() * y;
	return(ip(0));
}

// **********************************************************//
//      Multiply list of matrix X and list of vector B       //
// **********************************************************//
// [[Rcpp::export]]
arma::mat MultiplyXB(List X, arma::vec B){
    int A = X.size();
    arma::mat XB = arma::zeros(A, A);
    for (unsigned int n = 0; n < A; n++) {
        arma::mat X_n = X[n];
        arma::vec rows = X_n * B;
        XB.row(n) = rows.t();
    }
    return XB;
}



// **********************************************************//
//        Topic and Word contribution in update of Z         //
// **********************************************************//
// [[Rcpp::export]]
NumericVector TopicWord(int K, NumericVector tabledk, NumericVector tableWd, 
                        NumericVector tableCd, NumericVector tablek, int N,
                        NumericVector alphas, double beta, int V){
    NumericVector num1 = tablek + alphas[0]/K;
    double denom1 = N + alphas[0];
    NumericVector num2 = log(tableWd+beta/V);
    NumericVector denom2 = log(tablek+beta);
    NumericVector num3 = tableCd + alphas[1] * num1/denom1;
    double denom3 = sum(tableCd) + alphas[1];
    NumericVector consts = num2-denom2+log(tabledk + alphas[2] * num3/denom3);
    return consts;
}

// **********************************************************//
//        Topic and Word contribution in update of Z         //
// **********************************************************//
// [[Rcpp::export]]
NumericVector TopicWord0(int K, NumericVector tableCd, NumericVector tablek, int N,
                        NumericVector alphas, double beta, int V){
    NumericVector num1 = tablek + alphas[0]/K;
    double denom1 = N + alphas[0];
    NumericVector num3 = tableCd + alphas[1] * num1/denom1;
    double denom3 = sum(tableCd) + alphas[1];
    NumericVector consts = log(alphas[2] * num3/denom3);
    return consts;
}

// **********************************************************//
//        Topic and Word contribution in update of Z         //
// **********************************************************//
// [[Rcpp::export]]
double Topicpart(int K, IntegerVector z_d, NumericVector tableCd,
                  NumericVector tablek, NumericVector alphas){
    IntegerVector table_topics = tabulateC(z_d, K);
    double denom3 = sum(tableCd) + alphas[1] -1;
    double denom1 = sum(tablek) -1;
    double consts = 0;
    for (unsigned int k = 0; k < z_d.size(); k++) {
    	int zdn = z_d[k]-1;
        double dz = table_topics[zdn]-1;
        double num1 = tablek[zdn] - 1 + alphas[0]/K;
        double num3 = tableCd[zdn]-1+alphas[1]*num1/denom1;
        double ratio = alphas[2] * num3/ denom3;
        consts += log(dz+ratio);
    }
    return consts;
}

// **********************************************************//
//                    Multinomial Sampler                    //
// **********************************************************//
// [[Rcpp::export]]
int lmultinom (NumericVector lprops) {
    NumericVector props = expconst(lprops);
    IntegerVector samp = callRMultinom(props);
    for (unsigned int j = 0; j < props.size(); j++) {
        if (samp[j] == 1) {
            return j+1;
        }
    }
    return 1;
}

// **********************************************************//
//         Resampling the augmented data J_a (Sec 3.1)       //
// **********************************************************//
// [[Rcpp::export]]
NumericVector u_Gibbs(NumericVector u_di, NumericVector lambda_di, double delta, int j) {
    NumericVector prob(2);
    NumericVector u_di0 = u_di;
    u_di0[j-1] = 0;
    double sumu0 = sum(u_di0);
    prob[1] = delta+lambda_di[j-1];
    if (sumu0 == 0) {
        prob[0] = -1.0/0.0;
    }
    return prob;
}

// **********************************************************//
//              Likelihood evaluation of Edgepart            //
// **********************************************************//
// [[Rcpp::export]]
double Edgepart(arma::mat u, arma::mat lambda, double delta){
    double edgesum = 0;
    for (unsigned int i = 0; i < u.n_rows; i++) {
        arma::vec normal = arma::zeros(u.n_rows-1);
        double prob = 0;
        int iter = 0;
        for (unsigned int j = 0; j < u.n_rows; j++) {
            if (i != j) {
                double pre = delta+lambda(i, j);
                if (pre > 35) {
                    normal[iter] = pre;
                } else {
                    if (pre < -10) {
                        normal[iter] = exp(pre);
                    } else {
                        normal[iter] = log(exp(pre)+1);
                    }
                }
                prob += (delta+lambda(i, j))*u(i, j);
                iter = iter+1;
            }
        }
        double sumnorm = sum(normal);
        double normalizer = 0;
        if (sumnorm >= 13) {
            normalizer = sumnorm;
        } else {
            if (exp(sumnorm) <= 1) {
                normalizer = -745;
            } else {
                normalizer = log(exp(sumnorm)-1);
            }
        }
        edgesum += prob-normalizer;
    }
    return edgesum;
}


// **********************************************************//
//              Likelihood evaluation of Edgepart            //
// **********************************************************//
// [[Rcpp::export]]
double Edgepartsum(List X, arma::vec B, arma::mat u, double delta){
    arma::mat lambda = MultiplyXB(X, B);
    double edgesum = Edgepart(u, lambda, delta);
    return edgesum;
}


// **********************************************************//
//                Multiply matrix X and vector B             //
// **********************************************************//
// [[Rcpp::export]]
NumericVector mu_vec(NumericVector timemat, int A, NumericVector eta) {
    NumericVector mu(A);
    double etatime = timemat[0] * eta[A] + timemat[1] * eta[A+1];
    for (unsigned int a = 0; a < A; a++) {
      mu[a] = eta[a];
  }
  return mu + etatime;
}

// **********************************************************//
//                Multiply matrix X and vector B             //
// **********************************************************//
// [[Rcpp::export]]
NumericMatrix mu_mat(NumericMatrix timemat, int A, NumericMatrix eta, IntegerVector cd) {
    int D = cd.size();
    NumericMatrix mu(D, A);
    for (unsigned int d = 0; d < D; d++) {
        mu(d,_) = mu_vec(timemat(d,_), A, eta(cd[d]-1,_));
    }
    return mu;
}

// **********************************************************//
//              Likelihood evaluation of Timepart            //
// **********************************************************//
// [[Rcpp::export]]
double Timepart(NumericVector mu, double sigma_tau, double a_d, double t_d){
    unsigned int observed = a_d-1;
    double timesum = R::dnorm(log(t_d), mu[observed], sigma_tau, TRUE);
    for (unsigned int i = 0; i < mu.size(); i++) {
        if (i != observed) {
            timesum += R::pnorm(log(t_d), mu[i], sigma_tau, FALSE, TRUE);
        }
    }
    return timesum;
}


// **********************************************************//
//              Likelihood evaluation of Timepart            //
// **********************************************************//
// [[Rcpp::export]]
double Timepartsum(NumericMatrix mumat, double sigma_tau,
                   IntegerVector senders, NumericVector timestamps){
    int D = senders.size();
    double timesum = 0;
    for (unsigned int d = 0; d < D; d++) {
        double a_d = senders[d];
        timesum += Timepart(mumat(d,_), sigma_tau, a_d, timestamps[d]);
    }
    return timesum;
}

// **********************************************************//
//                    Cache time intervals                   //
// **********************************************************//
// [[Rcpp::export]]
List timefinder (NumericVector timestamps, IntegerVector edgetrim, double timeunit) {
    int D = timestamps.size();
    List out(D);
    for (int d = min(edgetrim)-1; d < D; d++) {
        double time0 = timestamps[d-1];
        double time1 = time0-24*timeunit;
        double time2 = time0-96*timeunit;
        double time3 = time0-384*timeunit;
        int id1 = which_num(time1, timestamps);
        int id2 = which_num(time2, timestamps);
        int id3 = which_num(time3, timestamps);
        arma::mat intervals(3, 2);
        intervals(0,0) = id1;
        intervals(0,1) = d;
        intervals(1,0) = id2;
        intervals(1,1) = id1-1;
        intervals(2,0) = id3;
        intervals(2,1) = id2-1;
        out[d] = intervals;		
    }
    return out;
}

