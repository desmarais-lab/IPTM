#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double logPLC(List finalZ, List zcnew, IntegerVector currentC, List tableW, List corpusW,
	      List textlist, double alpha, NumericVector mvec, double delta,
	      NumericVector nvec) {
  NumericVector out(finalZ.size());
  for (int d = 0; d < finalZ.size(); d++) {
    IntegerVector finalZd = finalZ[d];
    NumericVector out2 = finalZd.size();
    IntegerVector zcnewd = zcnew[currentC[d] - 1];
    IntegerVector textlistd = textlist[d];
    int ifelse1 = 0;
    int ifelse2 = 0;
    for (int n = 0; n < finalZd.size(); n++) {
      IntegerVector tableWdn = tableW[finalZd[n] - 1];
      IntegerVector corpusWdn = corpusW[finalZd[n] - 1];
      IntegerVector ucorpusWdn = unique(corpusWdn);
      if (zcnewd[finalZd[n]-1] > 0) {
	ifelse1 = 1;
      } else {
	ifelse1 = 0;
      }
      if (tableWdn[textlistd[n] - 1] > 0) {
	ifelse2 = 1;
      } else {
	ifelse2 = 0;
      }
      out2[n] = log(zcnewd[finalZd[n] - 1] - ifelse1 + alpha * mvec[finalZd[n] - 1]) -
	log(finalZd.size() - 1 + alpha) +
	log(tableWdn[textlistd[n] - 1] - ifelse2 + delta * nvec[n]) -
	log(corpusWdn.size() - ucorpusWdn.size() + delta);
    }
    out[d] = sum(out2);
  }
  return sum(out);
}

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
List currentZ2(int c, IntegerVector currentClist, List currentZ) {
  List out(currentClist.size());
  for (int d = 0; d < currentClist.size(); d++) {
    IntegerMatrix currentZd = currentZ[currentClist[d] - 1];
    out[d] = currentZd(c - 1, _);
  }
  return out;
}
// [[Rcpp::export]]
List finalZ2(IntegerVector currentC, List currentZ) {
  List out(currentC.size());
  for (int d = 0; d < currentC.size(); d++){
    IntegerMatrix currentZd = currentZ[d];
    out[d] = currentZd(currentC[d] - 1, _);
  }
  return out;
}
// [[Rcpp::export]]
NumericMatrix allxmatC(NumericMatrix edge, List pre, int i, IntegerVector A) {
  NumericMatrix out(A.size(), 4);
  IntegerVector A1(A.size() - 1);
  int iter = 0;
  for (int b = 0; b < A.size(); b++) {
    if (A[b] != i) {
      A1[iter] = A[b];
    }
    if (A[b] != i) {
      iter += 1;
    }
  }
  for (int a = 0; a < A1.size(); a++) {
    int j = A1[a];
    IntegerVector A2(A1.size() - 1);
    int it = 0;
    for (int c = 0; c < A1.size(); c++) {
      if (A1[c] != j) {
	A2[it] = A1[c];
      }
      if (A1[c]!=j) {
	it += 1;
      }
    }
    List ilist = pre[i - 1];
    List jlist = pre[j - 1];
    NumericVector ijlist = ilist[j - 1];
    NumericVector jilist = jlist[i - 1];
    double send = 0;
    double receive = 0;
    NumericVector triangle(A2.size());
    for (int ii = 0; ii < ijlist.size(); ii++) {
      send += exp(-0.05 * ijlist[ii]);
    }
    for (int jj = 0; jj < jilist.size(); jj++) {
      receive += exp(-0.05 * jilist[jj]);
    }
    for (int h = 0; h < A2.size(); h++) {
      double triangle1 = 0;
      double triangle2 = 0;
      double triangle3 = 0;
      double triangle4 = 0;
      List alist = pre[A2[h] - 1];
      NumericVector ailist = alist[i - 1];
      NumericVector ajlist = alist[j - 1];
      NumericVector ialist = ilist[A2[h] - 1];
      NumericVector jalist = jlist[A2[h]-1];
      for (int i = 0; i<ialist.size(); i++) {
	triangle1 += exp(-0.05 * ialist[i]);
	triangle3 += exp(-0.05 * ailist[i]);
      }
      for (int j = 0; j < ajlist.size(); j++) {
	triangle2 += exp(-0.05 * ajlist[j]);
	triangle4 += exp(-0.05 * jalist[j]);
      }
      triangle[h] = triangle1 * triangle2 + triangle3 * triangle4 +
	triangle3 * triangle2 + triangle1 * triangle4;
    }
    out(a, _) = NumericVector::create(1, send, receive, sum(triangle));
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix allxmat(NumericMatrix edge, IntegerVector node, List histlist,
		      int sender, double lambda) {
  // using the list of history, calculate the four time-weighted
  // network statistics (intercept, send, receive, triangles) for specific sender
  // and all the rest possible receivers
  NumericMatrix out(node.size() - 1, 4);
  IntegerVector A1(node.size() - 1);
  int iter = 0;
  for (int b = 0; b < node.size(); b++) {
    if (node[b]!=sender) {
      A1[iter] = node[b];
    }
    if (node[b]!=sender) {
      iter += 1;
    }
  }
  for (int a = 0; a < A1.size(); a++) {
    int receiver = A1[a];
    IntegerVector A2(A1.size() - 1);
    int it = 0;
    for (int d = 0; d < A1.size(); d++) {
      if (A1[d] != receiver) {
	A2[it] = A1[d];
      }
      if (A1[d]!=receiver) {
	it += 1;
      }
    }
    List ilist = histlist[sender - 1];
    List jlist = histlist[receiver - 1];
    NumericVector ijlist = ilist[receiver - 1];
    NumericVector jilist = jlist[sender - 1];
    double send = 0;
    double receive = 0;
    NumericVector triangle(A2.size());
    for (int ii=0; ii < ijlist.size(); ii++) {
      send += exp(-lambda*ijlist[ii]);
    }
    for (int jj=0; jj<jilist.size(); jj++) {
      receive = receive + exp(-lambda * jilist[jj]);
    }
    for (int h=0; h < A2.size(); h++) {
      double triangle1 = 0;
      double triangle2 = 0;
      double triangle3 = 0;
      double triangle4 = 0;
      List alist = histlist[A2[h] - 1];
      NumericVector ailist = alist[sender - 1];
      NumericVector ajlist = alist[receiver - 1];
      NumericVector ialist = ilist[A2[h] - 1];
      NumericVector jalist = jlist[A2[h] - 1];
      for (int i = 0; i < ialist.size(); i++) {
	triangle1 += exp(-lambda * ialist[i]);
      }
      for (int j=0; j<ajlist.size(); j++) {
	triangle2 += exp(-lambda * ajlist[j]);
      }
      for (int i = 0; i < ailist.size(); i++) {
	triangle3 += exp(-lambda * ailist[i]);
      }
      for (int j=0; j < jalist.size(); j++) {
	triangle4 += exp(-lambda * jalist[j]);
      }
      triangle[h] = triangle1 * triangle2 + triangle3 * triangle4 + triangle3 *
	triangle2 + triangle1 * triangle4;
    }
    out(a,_) = NumericVector::create(1, send, receive, sum(triangle));
  }
  return out;
}

// [[Rcpp::export]]
List sortedZ(int nIP, IntegerVector currentC, List currentZ) {
  // sort the documents and corresponding topics according to the IPs
  List out(nIP);
  for (int IP = 1; IP < (nIP + 1); IP++){
    IntegerVector currentClist(currentC.size());
    int it = 0;
    for (int d = 0; d < currentC.size(); d++){
      if (currentC[d] == IP) {
	currentClist[it] = d;
      }
      if (currentC[d]==IP) {
	it = it + 1;
      }
    }
    List out2(it);
    for (int d2 = 0; d2 < it; d2++) {
      IntegerMatrix currentZd = currentZ[currentClist[d2]];
      out2[d2] = currentZd(IP - 1, _);
    }
    out[IP - 1] = out2;
  }
  return out;
}

// [[Rcpp::export]]
List finalZ(IntegerVector currentC, List currentZ) {
  // extract the finalized topic assignments for every word in every document
  // using the IP assignments
  List out(currentC.size());
  for (int d = 0; d < currentC.size(); d++){
    IntegerMatrix currentZd = currentZ[d];
    out[d] = currentZd(currentC[d] - 1, _);
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
