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
NumericMatrix allxmat(NumericMatrix edge, IntegerVector node, List histlist,
		      int sender, double lambda) {
  // using the list of history, calculate the four time-weighted
  // network statistics (intercept, send, receive, triangles) for specific sender
  // and all the rest possible receivers
NumericMatrix out(node.size()-1, 4);
  IntegerVector A1(node.size()-1);
  int iter=0;
  for (int b=0; b<node.size(); b++){		
    if (node[b]!=sender){A1[iter]=node[b];}
    if (node[b]!=sender){iter=iter+1;}
  }
  for (int a=0; a<A1.size(); a++){
    int receiver=A1[a];
    IntegerVector A2(A1.size()-1);
    int it=0;
    for (int d=0; d<A1.size(); d++){		
      if (A1[d]!=receiver){A2[it]=A1[d];}
      if (A1[d]!=receiver){it=it+1;}
    }
    List ilist=histlist[sender-1];
    List jlist=histlist[receiver-1];
    NumericVector ijlist=ilist[receiver-1];
    NumericVector jilist=jlist[sender-1];
    double send=0;
    double receive=0;
    NumericVector triangle(A2.size());
    for (int ii=0; ii<ijlist.size(); ii++){send=send+exp(-lambda*ijlist[ii]);}
    for (int jj=0; jj<jilist.size(); jj++){receive=receive+exp(-lambda*jilist[jj]);}
    for (int h=0; h<A2.size(); h++){
      double triangle1=0;
      double triangle2=0;
      double triangle3=0;
      double triangle4=0;
      List alist= histlist[A2[h]-1];
      NumericVector ailist=alist[sender-1];
      NumericVector ajlist=alist[receiver-1];
      NumericVector ialist=ilist[A2[h]-1];
      NumericVector jalist=jlist[A2[h]-1];
      for (int i=0; i<ialist.size(); i++){triangle1+=exp(-lambda*ialist[i]);}
      for (int j=0; j<ajlist.size(); j++){triangle2+=exp(-lambda*ajlist[j]);}
      for (int i=0; i<ailist.size(); i++){triangle3+=exp(-lambda*ailist[i]);}
      for (int j=0; j<jalist.size(); j++){triangle4+=exp(-lambda*jalist[j]);}
      triangle[h]= triangle1*triangle2+triangle3*triangle4+triangle3*triangle2+triangle1*triangle4;
    }
    out(a,_)=NumericVector::create(1, send, receive, sum(triangle));	
  }
  return out;
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

