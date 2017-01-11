library(mvtnorm)
library(MCMCpack)
library(inline)
library(Rcpp)

setwd('/Users/bomin8319/Desktop/SU16/Rcode/Vance')
#Inputs
########################
#1. edgelist : matrix of 3 columns (col1:sender, col2:receiver, col3=time in unix.time format)
load("Vance_edge.RData") 
edge <- edge
edge[,3]<-edge[,3]/3600

#2. nodelist : matrix with the first column being the ID of nodes (ID's starting from 1), NOTE: rest columns containing node info are optional -later
load("Vance_node.RData") 
node <- as.numeric(as.vector(node[,1]))

#3. textlist: list (of length=number of edges in total) containing the words in each edge (i.e. document)
load("Vance_text.RData") 
textlist <- textlist

#4. vocabulary: all vocabularies used over the corpus
load("Vance_vocab.RData")
vocabulary <- vocabulary
########################




cppFunction('
')

#3. sortedZ: sort the documents and corresponding topics according to the IPs
cppFunction('
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
     		IntegerMatrix currentZd=currentZ[currentClist[d2]];
     		out2[d2]=currentZd(IP-1,_);
  		}
  		out[IP-1]=out2;	
   }		
  return out;
}
')

#4. finalZ: extract the finalized topic assignments for every word in every document using the IP assignments
cppFunction('
List finalZ(IntegerVector currentC, List currentZ){
  List out(currentC.size());
  for (int d=0; d<currentC.size(); d++){
    IntegerMatrix currentZd=currentZ[d];
    out[d]=currentZd(currentC[d]-1,_);
  }
  return out;
}
')

#4. logWZ: calculate the log of unnormalized constant (corresponding to product of Eq.13 in Section 3.2) to check the convergence
logWZ <- function(nIP, K, textlist2, W, zcnew, tableW, alpha, mvec, delta, nvec){
	
	out1=0
	for (IP in 1:nIP){
	num1 <- zcnew[[IP]]+alpha*mvec	
	denom1 <- sum(zcnew[[IP]])+alpha
	out1 <- out1+sum(lgamma(num1))-lgamma(denom1)
	}
	out2=0
	for (k in 1:K){
	num2 <- tableW[[k]]+delta*nvec	
	denom2 <- sum(tableW[[k]])+delta
	out2 <- out2+sum(lgamma(num2))-lgamma(denom2)
	}
	return(out1+out2)
}

#5. Supdate: denominator part used for parameter optimization
cppFunction('
double Supdate(double alpha, IntegerVector ctable){
  double D=0;
  double S=0;
  for (int n = 1; n < (ctable.size()+1); n++){
    D+=1/(n-1+alpha);
    S+=ctable[n-1]*D;
  }
  return S;
}
')

#6. Skupdate: numerator part used for parameter optimization
cppFunction('
NumericVector Skupdate(NumericVector vec, List cktable){
  NumericVector s(vec.size());
  for (int k=0; k<vec.size(); k++){
    double d=0;
    IntegerVector nowtable=cktable[k];
    for (int n = 1; n < (nowtable.size()+1); n++){
      d+=1/(n-1+vec[k]);
      s[k]+=nowtable[n-1]*d;
    }
  }
  
  return s;
}
')

#13. multiplyXB: construct beta%*%x
cppFunction('
NumericVector multiplyXB(NumericMatrix allxmatlist, NumericVector beta){
	NumericVector out(allxmatlist.nrow());
	for (int i=0; i<allxmatlist.nrow(); i++){
		double sum=0;
		for (int j=0; j<beta.size(); j++){
			sum=sum+allxmatlist(i,j)*beta[j];
		}
	out[i]=sum;
	}
	return out;
}
')
###########################################################

unix.time(VanceMCMC<- MCMC(edge, node, textlist, vocabulary, 2, 10, 0.5, outer=2, n1=3, n2=3, n3=3300, burn=300, thin=3, opt=TRUE, seed=1))
save(VanceMCMC, file="VanceMCMC_asym1.RData")

