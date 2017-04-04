#include <Rcpp.h>

using namespace Rcpp; 

// [[Rcpp::export]]
List History(List edge, NumericMatrix pd, IntegerVector node, double when) {
  // Calculate the weighted time difference from previous interactions to certain time 'when'
  //
  // Args:
  //  edge: list of document information with 3 elements (sender, receiver, time)
  //  pd: distribution of interaction patterns for each document in the corpus
  //  node: nodelist containing the ID of nodes
  //  when: specific timepoint that we are calculating the time difference from
  //
  // Returns:
  //  Matrix of time differences between all nodes
  int nIP = pd.ncol();
  List IPmat(nIP);
  for (int IP = 1; IP < (nIP + 1); IP++) {
  	List IPlist_IP(4);
  	for (int l = 0; l < 4; l++){
  		NumericMatrix IP_l(node.size(), node.size());
  		IPlist_IP[l] = IP_l;
  	}
  		IPmat[IP - 1] = IPlist_IP;
  }
	int iter = 0;
	for (int d = 0; d < edge.size(); d++) {
		List document = edge[d];
		double time = document[2];
		if (time < when) {
		  iter = iter + 1;
		}
	}
	if (iter > 0) {
	  for (int i = 0; i < iter; i++) {
	    List document2 = edge[i];
	    int sender = document2[0];
	    IntegerVector receiver = document2[1];
	    double time = document2[2];
	    double time1 = when - 384;
		double time2 = when - 96;
		double time3 = when - 24; 
	    for (int r = 0; r < receiver.size(); r++){
	       for (int IP = 1; IP < (nIP + 1); IP++) {
  			List IPlist_IP = IPmat[IP - 1];
  			if (time < time1) {
  				NumericMatrix IP_l = IPlist_IP[0];
  				IP_l(sender - 1, receiver[r] -1) += pd(i, IP -1);	
  				IPlist_IP[0] = IP_l;
			}
  			if (time >= time1 && time < time2) {
  				NumericMatrix IP_l = IPlist_IP[1];
				IP_l(sender - 1, receiver[r] -1) += pd(i, IP -1);
				IPlist_IP[1] = IP_l;
			} 
  			if (time >= time2 && time < time3) {
  				NumericMatrix IP_l = IPlist_IP[2];
				IP_l(sender - 1, receiver[r] -1) += pd(i, IP -1);
				IPlist_IP[2] = IP_l;
			}  				
			if (time >= time3) {
				NumericMatrix IP_l = IPlist_IP[3];
				IP_l(sender - 1, receiver[r] -1) += pd(i, IP -1);
				IPlist_IP[3] = IP_l;
			}		
			IPmat[IP - 1] = IPlist_IP;
  			}
	   }
	 }
  }
	return IPmat;
}

// [[Rcpp::export]]
List Degree(List history, IntegerVector node, int sender) {
  // Calculate degree statistics (indegree and outdegree) given the history of interactions
  //
  // Args:
  //  history: list of document information with 3 elements (sender, receiver, time)
  //  node: nodelist containing the ID of nodes
  //  sender: specific timepoint that we are calculating the time difference from
  //
  // Returns:
  //  Degree network statistic for specific sender and all possible receivers
  int nIP = history.size();
  List IPmat(nIP);
 
 
  for (int IP = 1; IP < (nIP + 1); IP++) {
  	 NumericMatrix degreemat_IP(node.size(), 8);
  	 List historyIP = history[IP - 1];
     NumericVector degree(8); 
	 
   for (int b = 0; b < node.size(); b++) {
    int receiver = node[b];
    for (int l = 0; l < 4; l++) {
    	NumericMatrix historyIP_l = historyIP[l];
    	double send = historyIP_l(sender - 1, receiver - 1);
    	
   		NumericVector indegree(node.size());
    	for (int h = 0; h < node.size(); h++) {
     	 int third = node[h];
     	 double htor = historyIP_l(third - 1, receiver - 1);
    	  indegree[h] = htor;
     	 }
    	 degree[l + 4] = sum(indegree);			
    	 degree[l] = degree[l] + send;
      }
      degreemat_IP(b,_) = degree;
      }
      for (int l = 0; l < 4; l++) {
      degreemat_IP(_, l) = rep(max(degreemat_IP(_,l)), node.size());
      }
    IPmat[IP - 1] = degreemat_IP;
  }
  return IPmat;
}

// [[Rcpp::export]]
List Dyadic(List history, IntegerVector node, int sender) {
  // Calculate dyadic network statistics (send and receive) given the history of interactions
  //
  // Args:
  //  history: list
  //  node: nodelist containing the ID of nodes
  //  sender: sender of the document which we exclued from possible receiver set
  //
  // Returns: 
  //  Dyadic network statistic for specific sender and all possible receivers
  int nIP = history.size();
  List IPmat(nIP);  
  
  for (int IP = 1; IP < (nIP + 1); IP++) {
  	 NumericMatrix dyadicmat_IP(node.size(), 8);
  	 List historyIP = history[IP - 1];
     NumericVector dyadic(8); 
     for (int b = 0; b < node.size(); b++) {
    	int receiver = node[b];
        for (int l = 0; l < 4; l++) {
    		NumericMatrix historyIP_l = historyIP[l];
    		dyadic[l] = historyIP_l(sender - 1, receiver - 1);
    		dyadic[l + 4] = historyIP_l(receiver - 1, sender - 1);
    	}
      dyadicmat_IP(b, _) = dyadic;
      }
      IPmat[IP - 1] = dyadicmat_IP;
   } 
  return IPmat;
}  

// [[Rcpp::export]]
List Triadic(List history, IntegerVector node, int sender) {
  // Calculate Triadic network statistics (2-send, 2-receive, sibling, cosibling) given the history of interactions
  //
  // Args:
  //  history: list of weighted time differences between all nodes
  //  node: nodelist containing the ID of nodes
  //  sender: specific timepoint that we are calculating the time difference from
  //
  // Returns:
  //  Triadic network statistic for specific sender and all possible receivers
   int nIP = history.size();
   List IPmat(nIP);
   for (int IP = 1; IP < (nIP+1); IP++) {
      NumericMatrix triadmat_IP(node.size(), 64);
  	  List historyIP = history[IP - 1];
  	  NumericVector triadic(64); 
       
        for (int b = 0; b < node.size(); b++) {
        int receiver = node[b];
 
        NumericVector twosend(node.size());
        NumericVector tworeceive(node.size());
        NumericVector sibling(node.size());
        NumericVector cosibling(node.size()); 
        int iter1 = 0;
        for (int l = 0; l < 4; l++) {
        for (int m = 0; m < 4; m++){
       
         NumericMatrix historyIP_l = historyIP[l];
         NumericMatrix historyIP_m = historyIP[m];
       
       	 for (int h = 0; h < node.size(); h++) {
     	   int third = node[h];	
     	   double stoh = historyIP_l(sender - 1, third - 1);
      	   double htos = historyIP_l(third - 1, sender - 1); 
     	   double rtoh = historyIP_m(receiver - 1, third - 1);
      	   double htor = historyIP_m(third - 1, receiver - 1); 
      	   twosend[h] = stoh * htor;
      	   tworeceive[h] = htos * rtoh;
      	   sibling[h] = htos * htor;
      	   cosibling[h] = stoh * rtoh;
       	}
       	triadic[iter1] = sum(twosend);
       	triadic[iter1 + 16] = sum(tworeceive);
       	triadic[iter1 + 32] = sum(sibling);
       	triadic[iter1 + 48] = sum(cosibling);
        iter1 = iter1 + 1;
       }
     }
        triadmat_IP(b,_) = triadic;
        }
    IPmat[IP - 1] = triadmat_IP;
  }
  return IPmat;
}



// [[Rcpp::export]]
NumericVector MultiplyXB(NumericMatrix X, NumericVector beta){
  // Multiply the (N by P) matrix X and vector beta (of length P)
  //
  // Args:
  //  X: (N by P) matrix X
  //  beta: vector of length P
  //
  // Returns:
  //  The vector (of length N) with ith element correspond to X[i, ] %*% beta
  NumericVector XB(X.nrow());
  for (int i = 0; i < X.nrow(); i++) {
    double sum = 0;
    for (int j = 0; j < beta.size(); j++) {
      sum = sum + X(i, j) * beta[j];
    }
    XB[i] = sum;
  }
  return XB;
}

// [[Rcpp::export]]
NumericMatrix MultiplyXBList(List X, NumericVector beta){
  // Multiply the (D by N by P) list X and vector beta (of length P)
  //
  // Args:
  //  X: (D by N by P) list 
  //  beta: vector of length P
  //
  // Returns:
  //  The matrix (D by N) with (i,j)th element correspond to X[i,j, ] %*% beta
  NumericMatrix example = X[0];	
	NumericMatrix XB(X.size(), example.nrow());
	for (int d = 0; d < X.size(); d++) {
		NumericVector XB_d(example.nrow());
		NumericMatrix X_d = X[d];
		 for (int i = 0; i < example.nrow(); i++) {
  		  double sum = 0;
    	for (int j = 0; j < beta.size(); j++) {
     	 sum = sum + X_d(i, j) * beta[j];
    	  }
    	XB_d[i] = sum;
  		}
		XB(d, _) = XB_d;
	}
  return XB;
}

// [[Rcpp::export]]
double UpdateDenom(double alpha, IntegerVector nwordtable){
  // Denominator part used for parameter optimization
  //
  // Args:
  //  alpha: Dirichlet concentration prior for topic distribution
  //  nwordtable: table of the number of words over the corpus
  //
  // Returns:
  //  The number which goes into the denominator in the equation of parameter optimization
  double D = 0;
  double S = 0;
  for (int n = 1; n < (nwordtable.size() + 1); n++) {
    D += 1 / (n - 1 + alpha);
    S += nwordtable[n - 1] * D;
  }
  return S;
}

// [[Rcpp::export]]
NumericVector UpdateNum(NumericVector vec, List nKwordtable) {
  // Numerator part used for parameter optimization
  //
  // Args:
  //  vec: alpha*mvec from Dirichlet priors for topic distribution
  //  nKwordtable: table of the number of words in each topic (K) over the corpus
  //
  // Returns:
  //  The number which goes into the numerator in the equation of parameter optimization
  NumericVector s(vec.size());
  for (int k = 0; k < vec.size(); k++){
    double d = 0;
    IntegerVector newtable = nKwordtable[k];
    for (int n = 1; n < (newtable.size() + 1); n++) {
      d += 1 / (n - 1 + vec[k]);
      s[k] += newtable[n - 1] * d;
    }
  }
  return s;
}


// [[Rcpp::export]]
NumericMatrix WordInEqZ(int K, IntegerVector textlistd, List tableW, 
                        double beta, NumericVector nvec){
  // Calculate word-topic part of the equation used in multinomial sampling of Z
  //
  // Args:
  //  K: total number of topics specified by the user
  //  textlistd: list of text containing the words in specific document d
  //  tableW: summary table of topic-word assignments
  //  beta: Dirichlet concentration prior for topic-word distribution
  //  nvec: Dirichlet base prior for topic-word distribution
  //
  // Returns:
  //  The vector of constants representing word part of each K
  NumericMatrix consts(textlistd.size(), K);
	for (int k = 0; k < K; k++){
		NumericVector tablek = tableW[k];
		NumericVector num(textlistd.size());
		NumericVector denom(textlistd.size());
		for (int w = 0; w < textlistd.size(); w++){
			num[w] = log(tablek[textlistd[w] - 1] - (tablek[textlistd[w] - 1] > 0) + 
			         beta * nvec[textlistd[w] - 1]);  
	 		denom[w] = log(sum(tablek) - (tablek[textlistd[w] - 1] > 0) + beta);
		}
		consts(_,k) = num - denom;
	}
	return consts;
}

// [[Rcpp::export]]
double EdgeInEqZ(IntegerMatrix iJi, NumericMatrix lambda, double delta) {
	double edges = 0;
	for (int i = 0; i < iJi.nrow(); i++) {
		for (int j = 0; j < iJi.ncol(); j++) {
		  double deltalambda = exp(-delta * lambda(i, j));
		  if (deltalambda < 0.0000001) { deltalambda += 0.0000001;}
		  double deltalambda2 = 1 / deltalambda - 1;
		  if (deltalambda2 < 0.0000001) { deltalambda2 += 0.0000001;}
			edges = edges + iJi(i, j) * log(deltalambda2) - delta * lambda(i, j);
		}
	}
	return edges;
}

// [[Rcpp::export]]
double EdgeInEqZ2(IntegerMatrix iJi, NumericMatrix lambda, double delta) {
	double edges = 0;
	for (int i = 0; i < iJi.nrow(); i++) {
		for (int j = 0; j < iJi.ncol(); j++) {
		  double deltalambda = delta * lambda(i, j);
		  if (deltalambda < 0.0000001) { deltalambda += 0.0000001;}
			edges = edges + iJi(i, j) * log(deltalambda) - log(deltalambda + 1);
		}
	}
	return edges;
}

// [[Rcpp::export]]
double TimeInEqZ(NumericVector LambdaiJi, NumericVector LambdaiJi2, double tdiff, double delta) {
	double times = sum(log(delta * LambdaiJi));
	double obs = tdiff * delta * (sum(LambdaiJi) + sum(LambdaiJi2));
	double consts = times - obs;
	return consts;
}


// [[Rcpp::export]]
IntegerVector tabulateC(const IntegerVector& x, const unsigned max) {
  // C++ version of R function: tabulate()
  //
  // Args:
  //  x: a numeric vector, or a factor
  //  max: the number of bins to be used
  //
  // Returns:
  //  An integer vector counting the number of times each integer occurs in it
  IntegerVector counts(max);
  std::size_t n = x.size();
  for (std::size_t i=0; i < n; i++) {
    if (x[i] > 0 && x[i] <= max)
      counts[x[i] - 1]++;
  }
  return counts;
}


