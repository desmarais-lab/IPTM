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


