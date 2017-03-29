#include <Rcpp.h>

using namespace Rcpp; 
// [[Rcpp::export]]
List History(List edge, NumericMatrix pd, IntegerVector node, double when, int nIP) {
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
  NumericMatrix histmat(node.size(), node.size());
  NumericMatrix countmat(node.size(), node.size());
  List IPmat(nIP);
  for (int IP = 1; IP < (nIP + 1); IP++) {
  	NumericMatrix IPmat_IP(node.size(), node.size());
  		IPmat[IP - 1] = IPmat_IP;
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
	    for (int r = 0; r < receiver.size(); r++){
	    	countmat(sender - 1, receiver[r] - 1) = countmat(sender - 1, receiver[r] - 1) + 1;
	      	histmat(sender - 1, receiver[r] - 1) = histmat(sender - 1, receiver[r] - 1) + log(when - time);
	       for (int IP = 1; IP < (nIP + 1); IP++) {
  			NumericMatrix IPmat_IP = IPmat[IP - 1];
  			IPmat_IP(sender - 1, receiver[r] - 1) = IPmat_IP(sender - 1, receiver[r] - 1) + pd(i, IP -1);
  			IPmat[IP - 1] = IPmat_IP;
  			}
	   }
	 }
  }
	return List::create(countmat, histmat, IPmat);
}


// [[Rcpp::export]]
List Dyadic(List history, IntegerVector node, int sender, int nIP) {
  // Calculate dyadic network statistics (send and receive) given the history of interactions
  //
  // Args:
  //  history: list
  //  node: nodelist containing the ID of nodes
  //  sender: sender of the document which we exclued from possible receiver set
  //
  // Returns: 
  //  Dyadic network statistic for specific sender and all possible receivers

  List IPmat(nIP);
  for (int IP = 1; IP < (nIP + 1); IP++) {
  	NumericMatrix dyadicmat(node.size() - 1, 2);
  	IPmat[IP - 1] = dyadicmat;
  }
  IntegerVector node_i(node.size() - 1);
  int iter = 0;
  for (int a = 0; a < node.size(); a++) {
    if (node[a] != sender) {node_i[iter] = node[a];}
    if (node[a] != sender) {iter = iter + 1;}
  }
  for (int b = 0; b < node_i.size(); b++) {
    int receiver = node_i[b];
    for (int IP = 1; IP < (nIP+1); IP++) {
      NumericMatrix historyIP = history[IP - 1];
      NumericMatrix dyadicmat_IP = IPmat[IP - 1];
      double send = historyIP(sender - 1, receiver - 1);
      double receive = historyIP(receiver - 1, sender - 1);
      dyadicmat_IP(b, _) = NumericVector::create(send, receive);
      IPmat[IP - 1] = dyadicmat_IP;
    }
   } 
  return IPmat;
}  

// [[Rcpp::export]]
List Triadic(List history, IntegerVector node, int sender, int nIP) {
  // Calculate Triadic network statistics (2-send, 2-receive, sibling, cosibling) given the history of interactions
  //
  // Args:
  //  history: list of weighted time differences between all nodes
  //  node: nodelist containing the ID of nodes
  //  sender: specific timepoint that we are calculating the time difference from
  //
  // Returns:
  //  Triadic network statistic for specific sender and all possible receivers
  List IPmat(nIP);
  for (int IP = 1; IP < (nIP + 1); IP++) {
  	NumericMatrix triadmat(node.size() - 1, 4);
   	IPmat[IP - 1] = triadmat;
  }
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
    for (int IP = 1; IP < (nIP+1); IP++) {
      NumericMatrix historyIP = history[IP - 1];
      NumericMatrix triadmat_IP = IPmat[IP - 1];
    for (int h = 0; h < node_ij.size(); h++) {
      int third = node_ij[h];
      double stoh = historyIP(sender - 1, third - 1);
      double htos = historyIP(third - 1, sender - 1); 
      double rtoh = historyIP(receiver - 1, third - 1);
      double htor = historyIP(third - 1, receiver - 1);
      twosend[h] = stoh * htor;
      tworeceive[h] = htos * rtoh;
      sibling[h] = htos * htor;
      cosibling[h] = stoh * rtoh;
      }			
    triadmat_IP(b,_) = NumericVector::create(sum(twosend), sum(tworeceive), sum(sibling), sum(cosibling));
    IPmat[IP - 1] = triadmat_IP;
  }
  }
  return IPmat;
}

// [[Rcpp::export]]
List Degree(List history, IntegerVector node, int sender, int nIP) {
  // Calculate degree statistics (indegree and outdegree) given the history of interactions
  //
  // Args:
  //  history: list of document information with 3 elements (sender, receiver, time)
  //  node: nodelist containing the ID of nodes
  //  sender: specific timepoint that we are calculating the time difference from
  //
  // Returns:
  //  Degree network statistic for specific sender and all possible receivers
   List IPmat(nIP);
  for (int IP = 1; IP < (nIP + 1); IP++) {
  	NumericMatrix degreemat(node.size() - 1, 2);
   	IPmat[IP - 1] = degreemat;
  }
  IntegerVector node_i(node.size() - 1);
  int iter = 0;
  for (int a = 0; a < node.size(); a++) {
    if (node[a] != sender) {node_i[iter] = node[a];}
    if (node[a] != sender) {iter = iter + 1;}
  }
  for (int IP = 1; IP < (nIP + 1); IP++) {
      NumericMatrix historyIP = history[IP - 1];
      NumericMatrix degreemat_IP = IPmat[IP - 1];
  double outdegree = 0;
   for (int b = 0; b < node_i.size(); b++) {
    int receiver = node_i[b];
    double send = historyIP(sender - 1, receiver - 1);
    IntegerVector node_ij(node_i.size() - 1);		
    int iter2 = 0;
    for (int d = 0; d < node_i.size(); d++) {
      if (node_i[d]!=receiver) {node_ij[iter2] = node_i[d];}
      if (node_i[d]!=receiver) {iter2 = iter2 + 1;}
    }
    NumericVector indegree(node_ij.size());
    for (int h = 0; h < node_ij.size(); h++) {
      int third = node_ij[h];
      double htor = historyIP(third - 1, receiver - 1);
      indegree[h] = htor;
      }			
    outdegree = outdegree + send;
    degreemat_IP(b,_) = NumericVector::create(0, send + sum(indegree));
  }
  degreemat_IP(_,0) = rep(outdegree, node_i.size());
  IPmat[IP - 1] = degreemat_IP;
  }
  return IPmat;
}

// [[Rcpp::export]]
NumericMatrix DyadicTW(List history, IntegerVector node, int sender) {
  // Calculate dyadic network statistics (send and receive) given the history of interactions
  //
  // Args:
  //  history: list
  //  node: nodelist containing the ID of nodes
  //  sender: sender of the document which we exclued from possible receiver set
  //
  // Returns: 
  //  Dyadic network statistic for specific sender and all possible receivers
  NumericMatrix countmat = history[0];
  NumericMatrix timemat = history[1];
  NumericMatrix dyadicmat(node.size() - 1, 2);

  IntegerVector node_i(node.size() - 1);
  int iter = 0;
  for (int a = 0; a < node.size(); a++) {
    if (node[a] != sender) {node_i[iter] = node[a];}
    if (node[a] != sender) {iter = iter + 1;}
  }
  for (int b = 0; b < node_i.size(); b++) {
    int receiver = node_i[b];
    double sendtime = 0;
    double receivetime = 0;
    if (countmat(sender - 1, receiver - 1) > 0) {
    	sendtime = timemat(sender - 1, receiver - 1) / countmat(sender - 1, receiver - 1);
    } 
     if (countmat(receiver - 1, sender - 1) > 0) {
    	receivetime = timemat(receiver - 1, sender - 1) / countmat(receiver - 1, sender - 1);
    } 
    dyadicmat(b, _) = NumericVector::create(sendtime, receivetime);
   } 
  return dyadicmat;
}  

// [[Rcpp::export]]
NumericMatrix TriadicTW(List history, IntegerVector node, int sender) {
  // Calculate Triadic network statistics (2-send, 2-receive, sibling, cosibling) given the history of interactions
  //
  // Args:
  //  history: list of weighted time differences between all nodes
  //  node: nodelist containing the ID of nodes
  //  sender: specific timepoint that we are calculating the time difference from
  //
  // Returns:
  //  Triadic network statistic for specific sender and all possible receivers
  NumericMatrix countmat = history[0];
  NumericMatrix timemat = history[1];
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
    IntegerVector twosendTW(node_ij.size());
    NumericVector tworeceive(node_ij.size());
    IntegerVector tworeceiveTW(node_ij.size());
    NumericVector sibling(node_ij.size());
    IntegerVector siblingTW(node_ij.size());
    NumericVector cosibling(node_ij.size());
    IntegerVector cosiblingTW(node_ij.size());
   
    for (int h = 0; h < node_ij.size(); h++) {
      int third = node_ij[h];
      double shtime = timemat(sender - 1, third - 1);
      double hstime = timemat(third - 1, sender - 1);
      double rhtime = timemat(receiver - 1, third - 1);
      double hrtime = timemat(third - 1, receiver - 1);    
      int shcount = countmat(sender - 1, third - 1);
      int hscount = countmat(third - 1, sender - 1);
      int rhcount = countmat(receiver - 1, third - 1);
      int hrcount = countmat(third - 1, receiver - 1);
      twosend[h] = shtime + hrtime;
      twosendTW[h] = shcount + hrcount;
      tworeceive[h] = hstime + rhtime;
      tworeceiveTW[h] = hscount + rhcount;
      sibling[h] = hstime + hrtime;
      siblingTW[h] = hscount + hrcount;
      cosibling[h] = shtime + rhtime;
      cosiblingTW[h] = shcount + rhcount;
      }	
      NumericVector triadic(4);
      if (sum(twosendTW) > 0) {
      	triadic[0] = sum(twosend) / sum(twosendTW);
      	}
      if (sum(tworeceiveTW) > 0) {
      	triadic[1] = sum(tworeceive) / sum(tworeceiveTW);
      	}
      if (sum(siblingTW) > 0) {
      	triadic[2] = sum(sibling) / sum(siblingTW);
      	}
      if (sum(cosiblingTW) > 0) {
      	triadic[3] = sum(cosibling) / sum(cosiblingTW);
      	}
    triadmat(b,_) = triadic;
      }
  return triadmat;
}

// [[Rcpp::export]]
NumericMatrix DegreeTW(List history, IntegerVector node, int sender) {
  // Calculate degree statistics (indegree and outdegree) given the history of interactions
  //
  // Args:
  //  history: list of document information with 3 elements (sender, receiver, time)
  //  node: nodelist containing the ID of nodes
  //  sender: specific timepoint that we are calculating the time difference from
  //
  // Returns:
  //  Degree network statistic for specific sender and all possible receivers
  NumericMatrix countmat = history[0];
  NumericMatrix timemat = history[1];
  NumericMatrix degreemat(node.size() - 1, 2);
  
  IntegerVector node_i(node.size() - 1);
  int iter = 0;
  for (int a = 0; a < node.size(); a++) {
    if (node[a] != sender) {node_i[iter] = node[a];}
    if (node[a] != sender) {iter = iter + 1;}
  }
   double outtime = 0;
   int outcount = 0;
   for (int b = 0; b < node_i.size(); b++) {
    int receiver = node_i[b];
    double sendtime = timemat(sender - 1, receiver - 1);
    int sendcount = countmat(sender - 1, receiver - 1);
    IntegerVector node_ij(node_i.size() - 1);		
    int iter2 = 0;
    for (int d = 0; d < node_i.size(); d++) {
      if (node_i[d] != receiver) {node_ij[iter2] = node_i[d];}
      if (node_i[d] != receiver) {iter2 = iter2 + 1;}
    }
    NumericVector intime(node_ij.size());
    IntegerVector incount(node_ij.size());
    for (int h = 0; h < node_ij.size(); h++) {
      int third = node_ij[h];
      double htortime = timemat(third - 1, receiver - 1);
      double htorcount = countmat(third - 1, receiver - 1);
      intime[h] = htortime;
      incount[h] = htorcount;
      }			
    outtime = outtime + sendtime;
    outcount = outcount + sendcount;
    double indegree = 0;
    if (sendcount + sum(incount) > 0) {
    	indegree = (sendtime + sum(intime)) / (sendcount + sum(incount));
    }
    degreemat(b,_) = NumericVector::create(0, indegree);
  }
  double outdegree = 0;
  if (outcount > 0) {
  	outdegree = outtime / outcount;
  }
  degreemat(_,0) = rep(outdegree, node_i.size());
  return degreemat;
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
double LambdaInEqZ(List iJi, NumericMatrix lambda, NumericVector LambdaiJi, NumericVector LambdaiJi2, double delta, double tdiff) {
	double edges = 0;
	double times = sum(log(LambdaiJi));
	double obs = tdiff * (sum(LambdaiJi) + sum(LambdaiJi2));
	double consts = 0;
	for (int i = 0; i < iJi.size(); i++) {
		IntegerVector Ji = iJi[i];
		for (int j = 0; j < iJi.size() - 1; j++) {
			edges = edges + Ji[j] * log(exp(delta * lambda(i, j) - 1)) - delta * lambda(i, j);
		}
	}
	consts = edges + times - obs;
	return consts;
}
