#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double Supdate(double alpha, IntegerVector ctable){
  // denominator part used for parameter optimization
  double D = 0;
  double S = 0;
  for (int n = 1; n < (ctable.size() + 1); n++) {
    D += 1 / (n - 1 + alpha);
    S += ctable[n - 1] * D;
  }
  return S;
}

// [[Rcpp::export]]
NumericVector Skupdate(NumericVector vec, List cktable) {
  // numerator part used for parameter optimization
  NumericVector s(vec.size());
  for (int k = 0; k < vec.size(); k++){
    double d = 0;
    IntegerVector nowtable = cktable[k];
    for (int n = 1; n < (nowtable.size() + 1); n++) {
      d += 1 / (n - 1 + vec[k]);
      s[k] += nowtable[n-1] * d;
    }
  }
  return s;
}

// [[Rcpp::export]]
List sortedZ(int nIP, IntegerVector currentC, List currentZ){
	List out(nIP);
   	for (int IP=1; IP<(nIP+1); IP++){
   		IntegerVector currentClist(currentC.size());
  		int it=0;
  		for (int d=0; d<currentC.size(); d++){
			if (currentC[d]==IP) {currentClist[it]=d;}
  			if (currentC[d]==IP) {it=it+1;}
  		}
    	List out2(it);
    	for (int d2=0; d2<it; d2++){
     		out2[d2]=currentZ[currentClist[d2]];
  		}
  		out[IP-1]=out2;	
   }		
  return out;
}

// [[Rcpp::export]]
List sortedC(int nIP, IntegerVector currentC, NumericMatrix edge){
	List out(nIP);
   	for (int IP=1; IP<(nIP+1); IP++){
   		IntegerVector currentClist(currentC.size());
  		int it=0;
  		for (int d=0; d<currentC.size(); d++){
			if (currentC[d]==IP) {currentClist[it]=d;}
  			if (currentC[d]==IP) {it=it+1;}
  		}
    	NumericMatrix out2(it, 3);
    	for (int d2=0; d2<it; d2++){
     		out2(d2,_)=edge(currentClist[d2], _);
  		}
  		out[IP-1]=out2;	
   }		
  return out;
}

// [[Rcpp::export]]
NumericVector multiplyXB(NumericMatrix allxmatlist, NumericVector beta){
  // construct beta %*% x
  NumericVector out(allxmatlist.nrow());
  for (int i = 0; i < allxmatlist.nrow(); i++) {
    double sum = 0;
    for (int j = 0; j < beta.size(); j++) {
      sum = sum + allxmatlist(i, j) * beta[j];
    }
    out[i] = sum;
  }
  return out;
}



// [[Rcpp::export]]
NumericVector betapartB(int nIP, List lambdai, List edgeC){
  NumericVector out(nIP);
  for (int i = 0; i < nIP; i++) {
  	double sum=0;
  	NumericMatrix edge = edgeC[i];
  	NumericMatrix lambda = lambdai[i];
  	IntegerVector receiver(edge.nrow());
  	for (int d=0; d<edge.nrow(); d++){
  		if (edge(d,0)>edge(d,1)) {
  			receiver[d]=edge(d,1);
  			} else {
  			receiver[d]=edge(d,1)-1;
  			}
  	NumericVector allreceiver = exp(lambda(d,_));
  	double sumreceiver=0;
  	for (int l = 0; l<allreceiver.size(); l++) {
  		sumreceiver = sumreceiver + allreceiver[l];
  	} 
  	sum = sum + lambda(d,receiver[d]-1) - log(sumreceiver);
  	}
  	out[i] = sum;
    }
  return out;
}

// [[Rcpp::export]]
NumericMatrix netstats(List allxmat, IntegerVector node, int sender) {
	NumericMatrix netstatmat(node.size()-1, 6);
	IntegerVector node2(node.size()-1);
	int iter=0;
	for (int a = 0; a < node.size(); a++) {
		if (node[a]!=sender) {node2[iter] = node[a];}
		if (node[a]!=sender) {iter = iter + 1;}
	}
	 	double outdegree = 0;
	for (int b = 0; b < node2.size(); b++) {
		int receiver = node2[b];
		List allxmats = allxmat[sender-1];
		List allxmatr = allxmat[receiver-1];
		double send = allxmats[receiver-1];
		double receive = allxmatr[sender-1];
		IntegerVector node3(node2.size()-1);		
		int iter2=0;
		for (int d = 0; d < node2.size(); d++) {
			if (node2[d]!=receiver) {node3[iter2] = node2[d];}
			if (node2[d]!=receiver) {iter2 = iter2 + 1;}
		}
		NumericVector triangle(node3.size());
		NumericVector indegree(node3.size());
		for (int h = 0; h < node3.size(); h++) {
			int third = node3[h];
			List allxmath = allxmat[third-1];
			double stoh = allxmats[third-1];
			double htos = allxmath[sender-1];
			double rtoh = allxmatr[third-1];
			double htor = allxmath[receiver-1];
			triangle[h] = stoh*htor+	htos*rtoh+htos*htor+	stoh*rtoh;
			indegree[h] = htor;}			
		outdegree = outdegree + send;
		netstatmat(b,_) = NumericVector::create(1, send, receive, sum(triangle), 0, send+sum(indegree));
	}
	netstatmat(_,4) = rep(outdegree, node2.size());
	return netstatmat;
}

// [[Rcpp::export]]
List nullmat(IntegerVector node) {
	List out2(node.size());
	List out(node.size());
	for (int i=0; i < node.size(); i++){
	 out2[i] = out;		
	}
	return out2;
}
	

// [[Rcpp::export]]
NumericMatrix wordpartZ(int K, IntegerVector textlistd, List tableW, double delta, NumericVector nvec){
	NumericMatrix out(textlistd.size(), K);
	for (int k=0; k<K; k++){
		NumericVector tablek = tableW[k];
		NumericVector num(textlistd.size());
		for (int l=0; l<textlistd.size(); l++){
			num[l] = log(tablek[textlistd[l]-1]- (tablek[textlistd[l]-1]>0)+delta*nvec[textlistd[l]-1]);  
		}
		double denom = log(sum(tablek) - sum(tablek>0)+delta);
		out(_,k) = num-rep(denom, textlistd.size());
	}
	return out;
}




