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
//               Prior calculations for IP sum               //
// **********************************************************//
// [[Rcpp::export]]
double priorsum (arma::mat var, arma::rowvec mu, arma::mat x) {
	double priorsum = sum(dmvnorm_arma(x, mu, var));
	return priorsum;
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
arma::mat rdirichlet_cpp(int num_samples, arma::vec alpha_m) {
	int dist_size = alpha_m.n_elem;
	arma::mat distribution = arma::zeros(num_samples, dist_size);
	
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
//                       Construct p.d matrix                //
// **********************************************************//
// [[Rcpp::export]]
NumericMatrix pdmat(List z, NumericVector l, int nIP) {
	int nDocs = z.size();
	NumericMatrix pd(nDocs, nIP);
	for (int IP = 1; IP < (nIP+1); IP++) {
		  IntegerVector IPvec = which_cpp(IP, l);
		for (int d = 0; d < nDocs; d++) {
			IntegerVector z_d = z[d];
		  for (unsigned int k = 0; k < IPvec.size(); k++) {
			pd(d, IP-1) += sum(z_d == IPvec[k]);
		  }
		  pd(d, IP-1) = pd(d, IP-1)/z_d.size();
		}
	}
	return pd;
}


// **********************************************************//
//              Construct the history of interaction         //
// **********************************************************//
// [[Rcpp::export]]
List History(List edge, NumericVector timestamps, NumericMatrix p_d, IntegerVector node, int d, double timeunit) {
  int nIP = p_d.ncol();
  int A = node.size();
  List IPmat(nIP);
  for (int IP = 0; IP < nIP; IP++) {
  	List IPlist_IP(3);
  	for (unsigned int l = 0; l < 3; l++){
  		NumericMatrix IP_l(A, A);
  		IPlist_IP[l] = IP_l;
  	}
  		IPmat[IP] = IPlist_IP;
  }
    double time1 = timestamps[d-2]-384*timeunit;
	double time2 = timestamps[d-2]-96*timeunit;
    double time3 = timestamps[d-2]-24*timeunit;
    
	  for (int i = 0; i < (d-1); i++) {
	    List document2 = edge[i];
	    int sender = document2[0];
	    IntegerVector receiver = document2[1];
	    double time = timestamps[i];
	    for (unsigned int r = 0; r < receiver.size(); r++){
	       for (int IP = 0; IP < nIP; IP++) {
  			  List IPlist_IP = IPmat[IP];
               if (time >= time1 && time < time2) {
                   NumericMatrix IP_l = IPlist_IP[2];
                   IP_l(sender-1, receiver[r]-1) += p_d(i, IP);
                   IPlist_IP[2] = IP_l;
               } else {
                   if (time >= time2 && time < time3) {
                       NumericMatrix IP_l = IPlist_IP[1];
                       IP_l(sender-1, receiver[r]-1) += p_d(i, IP);
                       IPlist_IP[1] = IP_l;
                   } else {
                       if (time >= time3) {
                       NumericMatrix IP_l = IPlist_IP[0];
                       IP_l(sender-1, receiver[r]-1) += p_d(i, IP);
                       IPlist_IP[0] = IP_l;
                       }
                   }
               }
			    IPmat[IP] = IPlist_IP;
  			}
         }
	  }
	return IPmat;
}


// **********************************************************//
//              Construct the history of interaction         //
// **********************************************************//
// [[Rcpp::export]]
List History2(List edge, NumericMatrix p_d, IntegerVector node, double when, double timeunit) {
    int nIP = p_d.ncol();
    List IPmat(nIP);
    for (int IP = 1; IP < (nIP+1); IP++) {
        List IPlist_IP(3);
        for (unsigned int l = 0; l < 3; l++){
            NumericMatrix IP_l(node.size(), node.size());
            IPlist_IP[l] = IP_l;
        }
        IPmat[IP-1] = IPlist_IP;
    }
    NumericVector timestamps(edge.size());
    for (unsigned int d = 0; d < edge.size(); d++) {
        List document = edge[d];
        timestamps[d] = Rcpp::as<double>(document[2]);
    }
    int iter = which_num(when, timestamps);
    
        for (int i = 0; i < iter; i++) {
            List document2 = edge[i];
            int sender = document2[0];
            IntegerVector receiver = document2[1];
            double time = Rcpp::as<double>(document2[2]);
            double time1 = when-384*timeunit;
            double time2 = when-96*timeunit;
            double time3 = when-24*timeunit;
            for (unsigned int r = 0; r < receiver.size(); r++){
                for (int IP = 0; IP < nIP; IP++) {
                    List IPlist_IP = IPmat[IP];
                    if (time >= time3) {
                        NumericMatrix IP_l = IPlist_IP[0];
                        IP_l(sender-1, receiver[r]-1) += p_d(i, IP);
                        IPlist_IP[0] = IP_l;
                    }
                    if (time >= time2 && time < time3) {
                        NumericMatrix IP_l = IPlist_IP[1];
                        IP_l(sender-1, receiver[r]-1) += p_d(i, IP);
                        IPlist_IP[1] = IP_l;
                    }  				
                    if (time >= time1 && time < time2) {
                        NumericMatrix IP_l = IPlist_IP[2];
                        IP_l(sender-1, receiver[r]-1) += p_d(i, IP);
                        IPlist_IP[2] = IP_l;
                    } 		
                    IPmat[IP] = IPlist_IP;
                }
            }
    }
    return IPmat;
}

// **********************************************************//
//                  Calculate Degree statistics              //
// **********************************************************//
// [[Rcpp::export]]
List Degree(List history, IntegerVector node, int sender) {
  int nIP = history.size();
  int A = node.size();
  List IPmat(nIP);
 
  for (int IP = 0; IP < nIP; IP++) {
  	 NumericMatrix degreemat_IP(A, 6);
  	 List historyIP = history[IP];
     NumericVector degree(6); 
	 
   for (unsigned int receiver = 0; receiver < A; receiver++) {
    for (unsigned int l = 0; l < 3; l++) {
    	NumericMatrix historyIP_l = historyIP[l];
    double send = historyIP_l(sender, receiver);
    	
   	NumericVector indegree(A);
    	for (unsigned int h = 0; h < A; h++) {
     	double htor = historyIP_l(h, receiver);
    	 indegree[h] = htor;
     	}
    	degree[l+3] = sum(indegree);
    	degree[l] = degree[l]+send;
      }
    degreemat_IP(receiver,_) = degree;
    }
   for (unsigned int l2 = 0; l2 < 3; l2++) {
    degreemat_IP(_, l2) = rep(max(degreemat_IP(_,l2)), A);
   }
  IPmat[IP] = degreemat_IP;
  }
  return IPmat;
}

// **********************************************************//
//                  Calculate Dyadic statistics              //
// **********************************************************//
// [[Rcpp::export]]
List Dyadic(List history, IntegerVector node, int sender) {
  int nIP = history.size();
  int A = node.size();
  List IPmat(nIP);  
  
  for (int IP = 0; IP < nIP; IP++) {
  	 NumericMatrix dyadicmat_IP(A, 6);
  	 List historyIP = history[IP];
     NumericVector dyadic(6); 
     for (unsigned int receiver = 0; receiver < A; receiver++) {
        for (unsigned int l = 0; l < 3; l++) {
    		NumericMatrix historyIP_l = historyIP[l];
    		dyadic[l] = historyIP_l(sender, receiver);
    		dyadic[l+3] = historyIP_l(receiver, sender);
    	}
      dyadicmat_IP(receiver, _) = dyadic;
    }
    IPmat[IP] = dyadicmat_IP;
  } 
  return IPmat;
}  

// **********************************************************//
//          Calculate (full 3 x 3) Triadic statistic         //
// **********************************************************//
// [[Rcpp::export]]
List Triadic(List history, IntegerVector node, int sender) {
   int nIP = history.size();
   int A = node.size();
   List IPmat(nIP);
   for (int IP = 0; IP < nIP; IP++) {
      NumericMatrix triadmat_IP(node.size(), 36);
  	  List historyIP = history[IP];
  	  NumericVector triadic(36); 
       
        for (unsigned int receiver = 0; receiver < A; receiver++) {
        NumericVector twosend(A);
        NumericVector tworeceive(A);
        NumericVector sibling(A);
        NumericVector cosibling(A); 
        int iter = 0;
        for (unsigned int l = 0; l < 3; l++) {
          for (unsigned int m = 0; m < 3; m++){
            NumericMatrix historyIP_l = historyIP[l];
            NumericMatrix historyIP_m = historyIP[m];
       	    for (unsigned int third = 0; third < A; third++) {
     	       double stoh = historyIP_l(sender, third);
      	       double htos = historyIP_l(third, sender); 
     	       double rtoh = historyIP_m(receiver, third);
      	       double htor = historyIP_m(third, receiver); 
      	        twosend[third] = stoh*htor;
      	        tworeceive[third] = htos*rtoh;
      	        sibling[third] = htos*htor;
      	        cosibling[third] = stoh*rtoh;
       	    }
       	    triadic[iter] = sum(twosend);
       	    triadic[iter+9] = sum(tworeceive);
       	    triadic[iter+18] = sum(sibling);
       	    triadic[iter+27] = sum(cosibling);
            iter = iter+1;
          }
        }
        triadmat_IP(receiver,_) = triadic;
      }
      
      NumericMatrix triadmat_IP2(A, 12);
   	  for (int t = 0; t < 4; t++) {
		int add = 9*t;
		for (int i = 0; i < A; i++) {
   	  		triadmat_IP2(i, 3*t+0) = triadmat_IP(i, add+0);
	    	triadmat_IP2(i, 3*t+1) = triadmat_IP(i, add+1)+triadmat_IP(i, add+3)+triadmat_IP(i, add+4);
	    	triadmat_IP2(i, 3*t+2) = triadmat_IP(i, add+2)+triadmat_IP(i, add+5)+triadmat_IP(i, add+6)+
	                        triadmat_IP(i, add+7)+triadmat_IP(i, add+8);
	    }                    
	}
      IPmat[IP] = triadmat_IP2/10;
 }
  return IPmat;
}


// **********************************************************//
//                    Network statistics                     //
// **********************************************************//
// [[Rcpp::export]]
List Netstats_cpp(List historyIP, IntegerVector node, IntegerVector netstat) {
	int A = node.size();
	int P = 3*(2*netstat[0]+2*netstat[1]+4*netstat[2]);
	int nIP = historyIP.size();
	List out(A);	
	for (int a = 0; a < A; a++) {
		List aout(nIP);
		for (int IP = 0; IP < nIP; IP++) {
		arma::mat netstatIP(A, P);
		aout[IP] = netstatIP;
		}
		int iter = 0;
		if (netstat[0] == 1) {
			List degree = Degree(historyIP, node, a);
			for (int IP = 0; IP < nIP; IP++){
				arma::mat aoutIP = aout[IP];
				arma::mat degreeIP = degree[IP];
				aoutIP.cols(iter, iter+5) = degreeIP;
			    aout[IP] = aoutIP;
			}
			iter += 6;
		}
		if (netstat[1] == 1) {
			List dyadic = Dyadic(historyIP, node, a);
			for (int IP = 0; IP < nIP; IP++){
				arma::mat aoutIP = aout[IP];
				arma::mat dyadicIP = dyadic[IP];
				aoutIP.cols(iter, iter+5) = dyadicIP;
			    aout[IP] = aoutIP;
			}
			iter += 6;			
		}	
		if (netstat[2] == 1) {
			List triadic = Triadic(historyIP, node, a);
			for (int IP = 0; IP < nIP; IP++){
				arma::mat aoutIP = aout[IP];
				arma::mat triadicIP = triadic[IP];
				aoutIP.cols(iter, iter+11) = triadicIP;
			    aout[IP] = aoutIP;
			}		
		}
		out[a] = aout;
	}
	return out;
}

// [[Rcpp::export]]
double inner (arma::vec x, arma::vec y) {
	arma::mat ip = x.t() * y;
	return(ip(0));
}

// **********************************************************//
//                Multiply matrix X and vector B             //
// **********************************************************//
// [[Rcpp::export]]
NumericMatrix ximat(NumericVector timemat, NumericMatrix eta1, NumericMatrix eta2) {
    int n_node = eta1.ncol();
    NumericMatrix xi(n_node, eta1.nrow());
    for (int IP = 0; IP < eta1.nrow(); IP++) {
      NumericVector etanode = eta1(IP, _);
      NumericVector etatime = eta2(IP, _);
      xi(_, IP) = etanode+inner(timemat, etatime);
  }
  return xi;
}

// **********************************************************//
//            Calculate xi over the entire corpus            //
// **********************************************************//
// [[Rcpp::export]]
List xi_all(NumericMatrix timemat, NumericMatrix eta1,NumericMatrix eta2, IntegerVector edgetrim) {
  List xi(timemat.nrow());
  for (IntegerVector::iterator it = edgetrim.begin(); it != edgetrim.end(); ++it) {
 		xi[*it-1] = ximat(timemat(*it-2, _), eta1, eta2);
	}
  return xi;
}

// **********************************************************//
//                Multiply matrix Y and vector B             //
// **********************************************************//
// [[Rcpp::export]]
NumericVector MultiplyYeta(NumericVector Y, NumericMatrix eta){
  NumericVector Yeta(eta.nrow());
    for (int IP = 0; IP < eta.nrow(); IP++) {
        arma::vec eta_IP = eta(IP,_);
        double sum = 0;
        for (unsigned int j = 0; j < Y.size(); j++) {
            sum = sum+Y[j]*eta_IP[j];
        }
        Yeta[IP] = sum;
    }
  return Yeta;
}

// **********************************************************//
//      Multiply list of matrix X and list of vector B       //
// **********************************************************//
// [[Rcpp::export]]
List MultiplyXB(List X, NumericMatrix B){
	List XB(B.nrow());
	for (int IP = 0; IP < B.nrow(); IP++) {
		arma::mat XB_IP(X.size(), X.size());
		arma::vec B_IP = B(IP, _);
		for (unsigned int n = 0; n < X.size(); n++) {
			List X_n = X[n];
			arma::mat X_n_IP = X_n[IP];
			arma::vec rows = X_n_IP*B_IP;
			XB_IP.row(n) = rows.t();
		}
		XB[IP] = XB_IP;
	}	
  return XB;
}

// **********************************************************//
//     Hyperparameter optimization of alpha-denominator    //
// **********************************************************//
// [[Rcpp::export]]
double UpdateDenom(double alpha, IntegerVector nwordtable){
 double D = 0;
 double S = 0;
  for (unsigned int n = 1; n < (nwordtable.size()+1); n++) {
    D += 1/(n-1+alpha);
    S += nwordtable[n-1]*D;
  }
  return S;
}

// **********************************************************//
//      Hyperparameter optimization of alpha-numerator     //
// **********************************************************//
// [[Rcpp::export]]
NumericVector UpdateNum(NumericVector vec, List nKwordtable) {
  NumericVector s(vec.size());
  for (unsigned int k = 0; k < vec.size(); k++){
   double d = 0;
    IntegerVector newtable = nKwordtable[k];
    for (unsigned int n = 1; n < (newtable.size()+1); n++) {
      if (vec[k] > 0) { 
      d += 1/(n-1+vec[k]);
      s[k] += newtable[n-1]*d;
      }
    }
  }
  return s;
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
//      Calculate lambda matrix for document d (Eq. (4))     //
// **********************************************************//
// [[Rcpp::export]]
arma::mat lambda_cpp(arma::vec p_d, List XB) {
  int nIP = XB.size();
  arma::mat example = XB[0];
  int node = example.n_rows;
  arma::mat lambdamat = arma::zeros(node, node);
  for (int IP = 0; IP < nIP; IP++) {
  	if (p_d[IP] > 0) {
    arma::mat XB_IP = XB[IP];
    lambdamat += p_d[IP]*XB_IP;
  	}
  }
  lambdamat.diag().zeros();
  return lambdamat;
}

// **********************************************************//
//              Calculate mu matrix for document d           //
// **********************************************************//
// [[Rcpp::export]]
double mu_cpp(arma::vec p_d, NumericVector xi) {
    int nIP = xi.size();
    double ximat = 0;
    for (int IP = 0; IP < nIP; IP++) {
        double pdIP = p_d[IP];
        if (pdIP > 0) {
           ximat += p_d[IP]*xi[IP];
        }
    }
    return ximat;
}

// **********************************************************//
//              Calculate mu matrix for document d           //
// **********************************************************//
// [[Rcpp::export]]
NumericVector mu_vec(arma::vec p_d, NumericMatrix xi) {
    int nIP = xi.ncol();
    NumericVector muvec(xi.nrow());
    for (int IP = 0; IP < nIP; IP++) {
        double pdIP = p_d[IP];
        for (int i = 0; i < xi.nrow(); i++) {
        	if (pdIP > 0) {
           		muvec[i] += p_d[IP]*xi(i,IP);
        	}
        }
    }
    return muvec;
}

// **********************************************************//
//            Calculate mu matrix for entire document        //
// **********************************************************//
// [[Rcpp::export]]
NumericMatrix mu_mat(NumericMatrix p_d, List xi, IntegerVector edgetrim) {
	NumericMatrix sample = xi[max(edgetrim)-1];
	NumericMatrix mumat(xi.size(), sample.nrow());
	for (IntegerVector::iterator it = edgetrim.begin(); it != edgetrim.end(); ++it) {
        int it2 = *it-1;
		mumat(it2, _) = mu_vec(p_d(it2, _), xi[it2]);
	}
    return mumat;
}


// **********************************************************//
//        Topic and Word contribution in update of Z         //
// **********************************************************//
// [[Rcpp::export]]
NumericVector TopicWord(int K, IntegerVector z_d, IntegerVector textlistd, IntegerMatrix tableW,
                        NumericVector alphamvec, double beta, int V, int w){
    IntegerVector table_topics = tabulateC(z_d, K);
    NumericVector num = log(tableW(_, textlistd[w-1]-1)+beta/V);
    NumericVector denom = log(rowSums(tableW)+beta);
    NumericVector consts = num - denom + log(as<NumericVector>(table_topics) + alphamvec);
    return consts;
}

// **********************************************************//
//        Topic and Word contribution in update of Z         //
// **********************************************************//
// [[Rcpp::export]]
NumericVector TopicWord0(int K, IntegerMatrix tableW, NumericVector alphamvec, double beta, int V){
    double num = log(beta/V);
    NumericVector denom = log(rowSums(tableW) + beta);
    NumericVector consts = num - denom + log(alphamvec);
    return consts;
}

// **********************************************************//
//         Resampling the augmented data J_a (Sec 3.1)       //
// **********************************************************//
// [[Rcpp::export]]
arma::vec u_Gibbs(arma::vec u_di, arma::vec lambda_di, double delta, int j) {
	arma::vec prob = arma::zeros(2);
	arma::vec u_di0 = u_di;
	u_di0[j-1] = 0;
	double sumu0 = sum(u_di0);
	prob[1] = delta+lambda_di[j-1];
	if (sumu0 == 0) {
		prob[0] = -arma::datum::inf;
	}
	return prob;
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
double Edgepartsum(List X, NumericVector p_d, NumericMatrix B, arma::mat u, double delta){
    List XB = MultiplyXB(X, B);
    arma::mat lambda = lambda_cpp(p_d, XB);
    double edgesum = Edgepart(u, lambda, delta);
    return edgesum;
}


// **********************************************************//
//              Likelihood evaluation of Timepart            //
// **********************************************************//
// [[Rcpp::export]]
double Timepart(arma::vec mu, double sigma_tau, double a_d, double t_d){
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
NumericVector Timepartindiv(arma::vec mu, double sigma_tau, double t_d){
    NumericVector timesum(mu.size());
    for (unsigned int i = 0; i < mu.size(); i++) {
    timesum[i] = R::dnorm(log(t_d), mu[i], sigma_tau, TRUE);
    for (unsigned int j = 0; j < mu.size(); j++) {
        if (j != i) {
            timesum[i] += R::pnorm(log(t_d), mu[j], sigma_tau, FALSE, TRUE);
        }
    }
    }
    return timesum;
}

// **********************************************************//
//              Likelihood evaluation of Timepart            //
// **********************************************************//
// [[Rcpp::export]]
double Timepartsum(NumericMatrix mumat, double sigma_tau, IntegerVector senders, NumericVector timeinc, IntegerVector edgetrim){
   double timesum = 0;
    for (IntegerVector::iterator it = edgetrim.begin(); it != edgetrim.end(); ++it) {
        int it2 = *it-1;
        double a_d = senders[it2];
		timesum += Timepart(mumat(it2,_), sigma_tau, a_d, timeinc[it2]);
	}
    return timesum;
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
//                    Cache time intervals                   //
// **********************************************************//
// [[Rcpp::export]]
List timefinder (NumericVector timestamps, IntegerVector edgetrim, double timeunit) {
    int D = timestamps.size();
    List out(D);
    for (unsigned int d = min(edgetrim)-1; d < D; d++) {
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

// **********************************************************//
//                 Use Array in Rcpp                         //
// **********************************************************//
// [[Rcpp::export]]
arma::cube histcache (int A, arma::vec senders, List edge) {
	int D = edge.size();
	arma::cube Array = arma::zeros<arma::cube>(A, A, D);
	for (unsigned int d = 0; d < D; d++) {
		List edged = edge[d];
	    int a = edged[0];
	    IntegerVector receiver = edged[1];
	   	for (unsigned int r = 0; r < receiver.size(); r++) {
	   		Array(a-1, receiver[r]-1, d) = 1;
	   	}	   					
	}	
	return Array;
}

// **********************************************************//
//                  Calculate Dyadic statistics              //
// **********************************************************//
// [[Rcpp::export]]
arma::cube dyadicstat (arma::cube cache, arma::mat intervals_d, int a, int A, arma::mat pd) {
	int nIP = pd.n_cols;
	arma::cube out = arma::zeros<arma::cube>(A, nIP, 6);
	for (unsigned int i = 0; i < 3; i++) {
		int mins = intervals_d(i, 0)-1;
		int maxs = intervals_d(i, 1)-1;
		arma::mat pdnew = pd.rows(mins, maxs);
		for (unsigned int r = 0; r < A; r++) {
			arma::colvec docsar = cache(arma::span(a), arma::span(r), arma::span(mins, maxs));
			arma::colvec docsra = cache(arma::span(r), arma::span(a), arma::span(mins, maxs));
			arma::rowvec sendr = docsar.t() * pdnew;
			arma::rowvec receiver = docsra.t() * pdnew;
			for (int IP = 0; IP < nIP; IP++) {
				out(r, IP, i) = sendr[IP];
				out(r, IP, i+3) = receiver[IP];
			}
		}		
	}	
	return out;
}

// **********************************************************//
//                  Calculate Dyadic statistics              //
// **********************************************************//
// [[Rcpp::export]]
arma::cube degreestat (arma::cube cache, arma::mat intervals_d, int a, int A, arma::mat pd) {
	int nIP = pd.n_cols;
	arma::cube out = arma::zeros<arma::cube>(A, nIP, 6);
	for (unsigned int i = 0; i < 3; i++) {
		int mins = intervals_d(i, 0)-1;
		int maxs = intervals_d(i, 1)-1;
		arma::rowvec sendr = arma::zeros<arma::rowvec>(nIP);
		arma::mat receiver = arma::zeros<arma::mat>(A, nIP);
		arma::mat pdnew = pd.rows(mins, maxs);
		for (unsigned int r = 0; r < A; r++) {
			arma::colvec docsar = cache(arma::span(a), arma::span(r), arma::span(mins, maxs));
			sendr += docsar.t() * pdnew;
			for (unsigned int h = 0; h < A; h++) {
				arma::colvec docsra = cache(arma::span(h), arma::span(r), arma::span(mins, maxs));
				receiver.row(r) += docsra.t() * pdnew;	
				}	
		}
		for (int IP = 0; IP < nIP; IP++) {
			for (unsigned int r = 0; r < A; r++) {
				out(r, IP, i) = sendr[IP];
				out(r, IP, i+3) = receiver(r,IP);
			}
		}		
	}	
	return out;
}

// **********************************************************//
//                  Calculate Triadic statistics              //
// **********************************************************//
// [[Rcpp::export]]
arma::cube triadicstat (arma::cube cache, arma::mat intervals_d, int a, int A, arma::mat pd) {
	int nIP = pd.n_cols;
	arma::cube out = arma::zeros<arma::cube>(A, nIP, 36);
	int it = 0;
	for (unsigned int i = 0; i < 3; i++) {
		for (unsigned int j = 0; j < 3; j++) {
		int mins1 = intervals_d(i, 0)-1;
		int maxs1 = intervals_d(i, 1)-1;
		int mins2 = intervals_d(j, 0)-1;
		int maxs2 = intervals_d(j, 1)-1;
		arma::mat tsend = arma::zeros<arma::mat>(A, nIP);
		arma::mat treceive = arma::zeros<arma::mat>(A, nIP);
		arma::mat sibling = arma::zeros<arma::mat>(A, nIP);
		arma::mat cosibling = arma::zeros<arma::mat>(A, nIP);
		arma::mat pdnew1 = pd.rows(mins1, maxs1);
		arma::mat pdnew2 = pd.rows(mins2, maxs2);
		for (unsigned int r = 0; r < A; r++) {
			for (unsigned int h = 0; h < A; h++) {
				arma::colvec docsah = cache(arma::span(a), arma::span(h), arma::span(mins1, maxs1));
				arma::colvec docsrh = cache(arma::span(r), arma::span(h), arma::span(mins2, maxs2));
				arma::colvec docsha = cache(arma::span(h), arma::span(a), arma::span(mins1, maxs1));
				arma::colvec docshr = cache(arma::span(h), arma::span(r), arma::span(mins2, maxs2));
				tsend.row(r) += (docsah.t() * pdnew1) % (docshr.t() * pdnew2);				
				treceive.row(r) += (docsha.t() * pdnew1) % (docsrh.t() * pdnew2);	
				sibling.row(r) += (docsha.t() * pdnew1) % (docshr.t() * pdnew2);				
				cosibling.row(r) += (docsah.t() * pdnew1) % (docsrh.t() * pdnew2);	
				}	
		}
		for (int IP = 0; IP < nIP; IP++) {
			for (unsigned int r = 0; r < A; r++) {
				out(r, IP, it) = tsend(r, IP);
				out(r, IP, it+9) = treceive(r,IP);
				out(r, IP, it+18) = sibling(r, IP);
				out(r, IP, it+27) = cosibling(r,IP);
			}
		}
		it +=1;		
	}
	}
	arma::cube out2 = arma::zeros<arma::cube>(A, nIP, 12);
	for (int t = 0; t < 4; t++) {
		int add = 9*t;
	for (int IP = 0; IP < nIP; IP++) {
		for (unsigned int r = 0; r < A; r++) {
			out2(r, IP, 3*t+0) = out(r, IP, add+0);
			out2(r, IP, 3*t+1) = out(r, IP, add+1) + out(r, IP, add+3)+ out(r, IP, add+4);
			out2(r, IP, 3*t+2) = out(r, IP, add+2) + out(r, IP, add+5)+ out(r, IP, add+6)+ out(r, IP, add+7)+ out(r, IP, add+8);
			}
		}
	}		
	return out2 / 10;
}




// **********************************************************//
//                    Network statistics                     //
// **********************************************************//
// [[Rcpp::export]]
List Netstats(arma::cube cache, arma::mat intervals_d, int A, IntegerVector netstat, arma::mat pd) {
	List out(A);
	int P = 3*(2*netstat[0]+2*netstat[1]+4*netstat[2]);
	int nIP = pd.n_cols;
	for (int a = 0; a < A; a++) {
		List aout(nIP);
		for (int IP = 0; IP < nIP; IP++) {
		arma::mat netstatIP(A, P);
		aout[IP] = netstatIP;
		}
	int iter = 0;	
	if (netstat[0] == 1) {
		arma::cube degree = degreestat(cache, intervals_d, a, A, pd);
			for (int IP = 0; IP < nIP; IP++){
				arma::mat aoutIP = aout[IP];
				arma::vec deg = degree(arma::span::all, arma::span(IP), arma::span::all);
				aoutIP.cols(iter, iter+5) = deg;
			    aout[IP] = aoutIP;
			}
			iter += 6;
	}	
	if (netstat[1] == 1) {
		arma::cube dyadic = dyadicstat(cache, intervals_d, a, A, pd);
		for (int IP = 0; IP < nIP; IP++){
				arma::mat aoutIP = aout[IP];
				arma::vec dyad = dyadic(arma::span::all, arma::span(IP), arma::span::all);
				aoutIP.cols(iter, iter+5) = dyad;
			    aout[IP] = aoutIP;
			}
			iter += 6;		
	}	
	if (netstat[2] == 1) {
		arma::cube triadic = triadicstat(cache, intervals_d, a, A, pd);
		for (int IP = 0; IP < nIP; IP++){
			arma::mat aoutIP = aout[IP];
			arma::vec triad = triadic(arma::span::all, arma::span(IP), arma::span::all);
			aoutIP.cols(iter, iter+11) = triad;
			aout[IP] = aoutIP;
		}
	}
	out[a] = aout;
	}
	return out;
}

