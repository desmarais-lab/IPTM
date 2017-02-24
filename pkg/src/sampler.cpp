#include <Rcpp.h>
using namespace Rcpp; 

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
List SortedZ(int nIP, IntegerVector currentC, List currentZ){
  // Sort the current assignment of topics corpus by current assignment of IPs 
  //
  // Args:
  //  nIP: total number of interaction patterns specified by the user
  //  currentC: current state of the assignment of interaction patterns
  //  currentZ: current state of the assignment of topics 
  //
  // Returns:
  //  Current assignment of topics across the documents under each interaction patern
  List sortedZ(nIP);
   	for (int IP = 1; IP < (nIP + 1); IP++){
   		IntegerVector currentClist(currentC.size());
  		int iter = 0;
  		for (int d = 0; d < currentC.size(); d++){
			  if (currentC[d] == IP) {
			  currentClist[iter] = d;
			  }
  		  if (currentC[d] == IP) {
  			iter = iter + 1;
  		  }
  		}
    	List sortedZ_IP(iter);
    	for (int dnew = 0; dnew < iter; dnew++){
    	  sortedZ_IP[dnew] = currentZ[currentClist[dnew]];
  		}
  		sortedZ[IP - 1] = sortedZ_IP;	
   }		
  return sortedZ;
}

// [[Rcpp::export]]
List SortedC(int nIP, IntegerVector currentC, List edge){
  // Sort the documents by current assignment of IPs
  //
  // Args:
  //  nIP: total number of interaction patterns specified by the user
  //  currentC: current state of the assignment of interaction patterns
  //  edge: list of document information with 3 elements (sender, receiver, time)
  //
  // Returns:
  //  The documents assigned to each interaction patern
  
	List sortedC(nIP);
   	for (int IP = 1; IP < (nIP+1); IP++){
   		IntegerVector currentClist(currentC.size());
  		int iter = 0;
  		for (int d = 0; d < currentC.size(); d++){
  		  if (currentC[d] == IP) {
			    currentClist[iter] = d;
			  }
  		  if (currentC[d] == IP) {
  		    iter = iter + 1;
  		    }
  		}
    	List sortedC_IP(iter);
    	for (int dnew = 0; dnew < iter; dnew++){
    	  sortedC_IP[dnew] = edge[currentClist[dnew]];
  		}
  		sortedC[IP - 1] = sortedC_IP;	
   }		
  return sortedC;
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
NumericVector BetaInEqB(int nIP, List lambdamat, List edgebyC){
  // Calculate Beta part of the equation used in Metropolis-Hasting sampling of B
  //
  // Args:
  //  nIP: total number of interaction patterns specified by the user
  //  lambdamat: list of lambdas (stochastic intensities) between all nodes by each IP
  //  edgebyC: list of edges sorted by each interaction patern
  //
  // Returns:
  //  The vector of constants representing beta part of each interaction pattern
  NumericVector consts(nIP);
  for (int i = 0; i < nIP; i++) {
  	double sum = 0;
  	List edge_IP = edgebyC[i];
  	NumericMatrix lambda_IP = lambdamat[i];
  	for (int d = 0; d < edge_IP.size(); d++){
  		List edge_d = edge_IP[d];
  		IntegerVector receiver_d = edge_d[1];
  		int sender_d = edge_d[0];
  		IntegerVector receiver_dnew(receiver_d.size());
  		NumericVector all_receiver = exp(lambda_IP(d,_));
  		double sum_receiver = 0;
  		for (int r = 0; r < all_receiver.size(); r++) {
  		  sum_receiver += all_receiver[r];
  		}
  		double sum_r = 0; 
  		for (int r2 = 0; r2 < receiver_d.size(); r2++){
  			if (sender_d > receiver_d[r2]) {
  			  receiver_dnew[r2] = receiver_d[r2];
  			} else {
  				receiver_dnew[r2] = receiver_d[r2] - 1;
  			}
  			sum_r += lambda_IP(d, receiver_dnew[r2] - 1);
  		}
  		sum = sum + sum_r - receiver_d.size() * log(sum_receiver);
  	}
  	consts[i] = sum;
    }
  return consts;
}

// [[Rcpp::export]]
NumericMatrix Timediff(List edge, IntegerVector node, double when, double lambda) {
  // Calculate the weighted time difference from previous interactions to certain time 'when'
  //
  // Args:
  //  edge: list of document information with 3 elements (sender, receiver, time)
  //  node: nodelist containing the ID of nodes
  //  when: specific timepoint that we are calculating the time difference from
  //  lambda: parameter of response speed with larger values indicating faster response
  //
  // Returns:
  //  Matrix of weighted time differences between all nodes
  NumericMatrix histmat(node.size(), node.size());
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
	    for (int r = 0; r < receiver.size(); r++){
	      histmat(sender - 1, receiver[r] - 1) = histmat(sender - 1, receiver[r] - 1) + exp(- lambda * (when - time));
	   }
	 }
  }
	return histmat;
}

// [[Rcpp::export]]
NumericMatrix Dyadic(NumericMatrix history, IntegerVector node, int sender) {
  // Calculate dyadic network statistics (send and receive) given the history of interactions
  //
  // Args:
  //  history: matrix of weighted time differences between all nodes
  //  node: nodelist containing the ID of nodes
  //  sender: sender of the document which we exclued from possible receiver set
  //
  // Returns: 
  //  Dyadic network statistic for specific sender and all possible receivers
  NumericMatrix dyadicmat(node.size() - 1, 2);
  IntegerVector node_i(node.size() - 1);
  int iter = 0;
  for (int a = 0; a < node.size(); a++) {
    if (node[a] != sender) {node_i[iter] = node[a];}
    if (node[a] != sender) {iter = iter + 1;}
  }
  for (int b = 0; b < node_i.size(); b++) {
    int receiver = node_i[b];
    double send = history(sender - 1, receiver - 1);
    double receive = history(receiver - 1, sender - 1);
    dyadicmat(b, _) = NumericVector::create(send, receive);
  }
  return dyadicmat;
}  

// [[Rcpp::export]]
NumericMatrix Triadic(NumericMatrix history, IntegerVector node, int sender) {
  // Calculate Triadic network statistics (2-send, 2-receive, sibling, cosibling) given the history of interactions
  //
  // Args:
  //  history: matrix of weighted time differences between all nodes
  //  node: nodelist containing the ID of nodes
  //  sender: specific timepoint that we are calculating the time difference from
  //
  // Returns:
  //  Triadic network statistic for specific sender and all possible receivers
  NumericMatrix triadmat(node.size() - 1, 4);
  IntegerVector node_i(node.size() - 1);
  int iter = 0;
  for (int a = 0; a < node.size(); a++) {
    if (node[a] != sender) {node_i[iter] = node[a];}
    if (node[a] != sender) {iter = iter + 1;}
  }
  for (int b = 0; b < node_i.size(); b++) {
    int receiver = node_i[b];
    IntegerVector node_ij(node_i.size() - 1);		
    int iter2 = 0;
    for (int d = 0; d < node_i.size(); d++) {
      if (node_i[d] != receiver) {node_ij[iter2] = node_i[d];}
      if (node_i[d] != receiver) {iter2 = iter2 + 1;}
    }
    NumericVector twosend(node_ij.size());
    NumericVector tworeceive(node_ij.size());
    NumericVector sibling(node_ij.size());
    NumericVector cosibling(node_ij.size());
    for (int h = 0; h < node_ij.size(); h++) {
      int third = node_ij[h];
      double stoh = history(sender - 1, third - 1);
      double htos = history(third - 1, sender - 1); 
      double rtoh = history(receiver - 1, third - 1);
      double htor = history(third - 1, receiver - 1);
      twosend[h] = stoh * htor;
      tworeceive[h] = htos * rtoh;
      sibling[h] = htos * htor;
      cosibling[h] = stoh * rtoh;
      }			
    triadmat(b,_) = NumericVector::create(sum(twosend), sum(tworeceive), sum(sibling), sum(cosibling));
  }
  return triadmat;
}

// [[Rcpp::export]]
NumericMatrix Degree(NumericMatrix history, IntegerVector node, int sender) {
  // Calculate degree statistics (indegree and outdegree) given the history of interactions
  //
  // Args:
  //  history: list of document information with 3 elements (sender, receiver, time)
  //  node: nodelist containing the ID of nodes
  //  sender: specific timepoint that we are calculating the time difference from
  //
  // Returns:
  //  Degree network statistic for specific sender and all possible receivers
  NumericMatrix degreemat(node.size() - 1, 2);
  IntegerVector node_i(node.size() - 1);
  int iter = 0;
  for (int a = 0; a < node.size(); a++) {
    if (node[a] != sender) {node_i[iter] = node[a];}
    if (node[a] != sender) {iter = iter + 1;}
  }
  double outdegree = 0;
  for (int b = 0; b < node_i.size(); b++) {
    int receiver = node_i[b];
    double send = history(sender - 1, receiver - 1);
    IntegerVector node_ij(node_i.size() - 1);		
    int iter2=0;
    for (int d = 0; d < node_i.size(); d++) {
      if (node_i[d]!=receiver) {node_ij[iter2] = node_i[d];}
      if (node_i[d]!=receiver) {iter2 = iter2 + 1;}
    }
    NumericVector indegree(node_ij.size());
    for (int h = 0; h < node_ij.size(); h++) {
      int third = node_ij[h];
      double htor = history(third - 1, receiver - 1);
      indegree[h] = htor;
      }			
    outdegree = outdegree + send;
    degreemat(b,_) = NumericVector::create(0, send+sum(indegree));
  }
  degreemat(_,0) = rep(outdegree, node_i.size());
  return degreemat;
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
		for (int w = 0; w < textlistd.size(); w++){
			num[w] = log(tablek[textlistd[w] - 1]-(tablek[textlistd[w] - 1] > 0) + 
			         beta*nvec[textlistd[w] - 1]);  
		}
		double denom = log(sum(tablek) - sum(tablek > 0)+beta);
		consts(_,k) = num - rep(denom, textlistd.size());
	}
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

// [[Rcpp::export]]
double logWZ(int nIP, int K, IntegerVector currentC, List currentZ, List textlist, 
             List tableW, NumericVector alpha, NumericMatrix mvec, double beta, 
             NumericVector nvec){
  // Calculate the log of unnormalized constant to check the convergence
  //
  // Args:
  //  nIP:
  //  K: the number of bins to be used
  //  currentC: current state of the assignment of interaction patterns 
  //  currentZ: current state of the assignment of topics 
  //  textlist: list of text containing the words in each document
  //  tableW: summary table of topic-word assignments
  //  alpha: Dirichlet concentration prior for document-topic distribution
  //  mvec: Dirichlet base prior for document-topic distribution
  //  beta: Dirichlet concentration prior for topic-word distribution
  //  nvec: Dirichlet base prior for topic-word distribution
  //
  // Returns:
  //  The log value of unnormalized constant
  double finalsum = 0;
	for (int d = 0; d < currentC.size(); d++) {
		IntegerVector currentZ_d = currentZ[d];
		IntegerVector k_table = tabulateC(currentZ_d, K);
		IntegerVector textlist_d = textlist[d];
		double part4 = log(sum(k_table) - 1 + alpha[currentC[d] - 1]);
		int iter = 0;
		for (int k = 0; k < currentZ_d.size(); k++) {
			IntegerVector table_Wk = tableW[currentZ_d[iter] - 1];
			double part1 = log(table_Wk[textlist_d[iter] - 1] - 1 + beta * nvec[textlist_d[iter] - 1]);
			double part2 = log(sum(table_Wk) - sum(table_Wk>0) + beta);
			double part3 = log(k_table[currentZ_d[iter] - 1] - 1 + 
                     alpha[currentC[d] - 1] * mvec(currentZ_d[iter] - 1,currentC[d] - 1));
			finalsum += part1 - part2 + part3 - part4;
			iter = iter + 1;
			}
		}
	return finalsum;
}
