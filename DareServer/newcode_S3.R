library(mvtnorm)
library(MCMCpack)
library(inline)
library(Rcpp)

#Inputs
########################
#1. edgelist : matrix of 3 columns (col1:sender, col2:receiver, col3=time in unix.time format)
load("Dare_edge.RData") 
edge <- edge
edge[,3]<-edge[,3]/3600

#2. nodelist : matrix with the first column being the ID of nodes (ID's starting from 1), NOTE: rest columns containing node info are optional -later
load("Dare_node.RData") 
node <- as.numeric(as.vector(node[,1]))

#3. textlist: list (of length=number of edges in total) containing the words in each edge (i.e. document)
load("Dare_text.RData") 
textlist <- textlist

#4. vocabulary: all vocabularies used over the corpus
load("Dare_vocab.RData")
vocabulary <- vocabulary
########################


#Functions
########################
#1. history: calculate the time difference from previous interactions to certain time 'beforetime'-> need to be faster... don't know how to do this in Rcpp
history <- function(edge, node, when){
	histlist <- lapply(node, function(v){lapply(node, function(u){numeric(0)})})
	edge2 <- matrix(edge[edge[,3] < when, ], byrow=FALSE, ncol=3)
	if (nrow(edge2)>0){
	for (d in 1:nrow(edge2)){
		histlist[[edge2[d,1]]][[edge2[d,2]]]<-c(histlist[[edge2[d,1]]][[edge2[d,2]]], when-(edge2[d ,3]))
	}
	}		
return(histlist)
}


#2. allxmat: using the list of history, calculate the four time-weighted network statistics (intercept, send, receive, triangles) for specific sender and all the rest possible receivers
cppFunction('
NumericMatrix allxmat(NumericMatrix edge, IntegerVector node, List histlist, int sender, double lambda){
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
}')

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

#7. parupdate: parameter optimization of alpha and mvec at the same time -> need to be faster
parupdate<- function(nIP, K, currentC, currentZ, alpha, mvec){
	vec <- alpha*mvec
	corpusCnew <- sortedZ(nIP,currentC, currentZ)
	zcnew <- lapply(1:nIP, function(IP){tabulate(unlist(corpusCnew[[IP]]), nbins=K)})
	clist <- vapply(1:nIP, function(IP){sum(zcnew[[IP]])}, c(1))
	
	iter=1
	ctable <- rep(0, max(clist))
	ctable[clist]<-tabulate(currentC)
	
	cklist <- list()
	for (IP in 1:nIP){
		cklist[[IP]] <- matrix(NA, nrow=length(corpusCnew[[IP]]), ncol=K)
		for (d in 1:length(corpusCnew[[IP]])){
			cklist[[IP]][d,]<-tabulate(corpusCnew[[IP]][[d]], nbins=K)
			}
	}
	cktable <- lapply(1:K, function(k){tabulate(unlist(sapply(1:nIP, function(IP){cklist[[IP]][,k]})))})
		
			
	while ((abs(alpha-sum(vec))>0.001) | (iter==1)){
	alpha <- sum(vec)
	S <- Supdate(alpha, ctable)		
	s <- Skupdate(vec, cktable)	
	vec<- vec*s/S
	iter=iter+1
	}
	return(vec)	
}	

#8. selected: find out the chosen category from the multinomial distribution
selected <- function(samples, proportions){
	chosen <-rmultinom(samples, 1, proportions)
	out <- vapply(1:samples, function(s){which(chosen[,s]==TRUE)}, c(1))
	return(out)
}

#9. betapartC: calculate the log of beta part used to obtain the constants for multinomial sampling of IP
betapartC <- function(nIP, lambdai, specificedge){
	const <- rep(NA, nIP)
	for (IP in 1:nIP){
	edge <- matrix(specificedge, ncol=3)
	receiver <- ifelse(edge[,1] > edge[,2], edge[,2],  edge[,2]-1)
	const[IP]<-lambdai[[IP]][receiver]-log(sum(exp(lambdai[[IP]])))
	}
    return(const)
}

#10. topicpartC: calculate the log of topic-IP part used to obtain the constants for multinomial sampling of IP
topicpartC <- function(nIP, currentZ, zc, alpha, mvec, document){
	const <- rep(NA, nIP)
	for (IP in 1:nIP){
	topics <- currentZ[[document]][IP,]
	num <- sum(log(zc[[IP]][topics]-ifelse(zc[[IP]][topics]>0, 1,0)+alpha*mvec[topics]))
	denom <-length(topics)*log(sum(zc[[IP]])-sum(ifelse(zc[[IP]]>0, 1, 0))+alpha)
	const[IP]<- num-denom
	}
	return(const)
} 

#11. topicpartZ: calculate the log of topic-IP part used to obtain the constants for multinomial sampling of K
topicpartZ <- function(currentIP, zcnew, alpha, mvec){
	const <- log(zcnew[[currentIP]]-ifelse(zcnew[[currentIP]]>0, 1, 0)+alpha*mvec)
	return(const)
}

#12. wordpartZ: calculate the log of word-topic part used to obtain the constants for multinomial sampling of K
wordpartZ <- function(K, textlistd, tableW, delta, nvec){
	const <- matrix(NA, nrow=length(textlistd), ncol=K)
	for (k in 1:K){
		num <-log(tableW[[k]][textlistd]-ifelse(tableW[[k]][textlistd]>0, 1, 0)+delta*nvec[textlistd])
		denom <-log(sum(tableW[[k]])-sum(ifelse(tableW[[k]]>0, 1, 0))+delta) 	
		const[,k]<- num-denom
	}
	return(const)
}

#13. betapartB: calculate the log of beta part used to obtain the constants for multinomial sampling of Beta
betapartB <- function(nIP, lambdai, edgeC){
	const <- rep(NA, nIP)
	for (IP in 1:nIP){
	edge <- matrix(edgeC[[IP]], ncol=3)
	receiver <- ifelse(edge[,1] > edge[,2], edge[,2],  edge[,2]-1)
	const[IP]<-sum(vapply(1:length(receiver), function(r){
			   lambdai[[IP]][r, receiver[r]]-log(sum(exp(lambdai[[IP]][r,])))
			   }, c(1)))
	}
    return(const)
}

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

MCMC =function(edge, node, textlist, vocabulary, nIP, K, delta_B, outer=200, n1=3, n2=3, n3=3300, burn=300, thin=3, opt=TRUE, seed=1){

	set.seed(seed)         	   
	
	#initialize alpha, mvec, delta, nvec, eta, lvec, and gammas
	alpha<- 50/K				
	mvec <- rep(1, K)/K       		
	delta=0.01					
	W = length(vocabulary)		 
	nvec <- rep(1, W)/W			
	eta = 50/nIP					
	lvec <- rep(1, nIP)/nIP			
	gammas <- rdirichlet(1,  eta*lvec)	

	#initialize theta 
	theta <- list()
	for (IP in 1:nIP){
	theta[[IP]]<- rdirichlet(1,alpha*mvec)							
	}

    #initialize C and Z 
    currentC <- selected(nrow(edge),gammas)
	currentZ <- lapply(1:nrow(edge), function(d){
				matrix(vapply(1:nIP, function(IP){
				selected(length(textlist[[d]]), theta[[IP]])
				}, rep(0, length(textlist[[d]]))), byrow=TRUE,nrow=nIP)
				})	  
	#initialize beta
	P=4							
	sigma=0.1					
	bmat <- list()
	for (IP in 1:nIP){
	  bmat[[IP]] <- matrix(NA, nrow=P, ncol=(n3-burn)/thin)
	  bmat[[IP]][,]<-rmvnorm(1, rep(0, P), sigma^2*diag(P))
	}
	  
	logWZmat <- c()							
	  
	if (opt==TRUE){
	  alphamat <- alpha
	  mvecmat <- mvec
	 }
  	  
	#start outer iteration	
	for(o in 1:outer){	
	cat("outer iteration = ", o, "\n") 

	  	if (opt==TRUE) {
	  	 	op <- options(warn = (-1))
	  	 	vec <- parupdate(nIP, K, currentC, currentZ, alpha, mvec)
	  	 	alpha <- sum(vec)
	  	    alphamat <- c(alphamat, alpha)
	  	 	mvec <- vec/alpha	
	  	    mvecmat <- rbind(mvecmat, mvec)
	  	 	options(op)
	  	 }						

	 cat("inner iteration 1", "\n") 
	 lambdai<-list() 
	 corpusC <-sortedZ(nIP, currentC, currentZ)
	 zc <-  lapply(1:nIP, function(IP){tabulate(unlist(corpusC[[IP]]), nbins=K)})
    	 for (i1 in 1:n1){
	 	#C update given Z and B - within each document d 	
	    for (d in 1:nrow(edge)){
	    	currentC2 <- currentC[-d]
			edgeC <- lapply(1:nIP, function(IP){edge[which(currentC2==IP),]})
			zc[[currentC[d]]] <- zc[[currentC[d]]]-tabulate(currentZ[[d]][currentC[d],], nbins=K)
			for (IP in 1:nIP){
	    		histlist <- history(edgeC[[IP]], node, edge[d,3])
	      		allxmatlist<-allxmat(edgeC[[IP]], node, histlist, edge[d,1], 0.05)
	      		lambdai[[IP]]<- multiplyXB(allxmatlist, rowMeans(bmat[[IP]]))
	      	}
	       	const<-log(gammas)+betapartC(nIP, lambdai, edge[d,])+topicpartC(nIP, currentZ, zc, alpha, mvec, d)
	        currentC[d]<- selected(1, exp(const))
	        zc[[currentC[d]]] <- zc[[currentC[d]]]+tabulate(currentZ[[d]][currentC[d],], nbins=K)
	        }
	    }
	    
	 cat("inner iteration 2", "\n")
	 textlist2 <- unlist(textlist)
	 corpusCnew <- sortedZ(nIP, currentC, currentZ)
	 zcnew <- lapply(1:nIP, function(IP){tabulate(unlist(corpusCnew[[IP]]), nbins=K)})
	 finalZlist <- finalZ(currentC, currentZ)
	 finalZlist2 <- unlist(finalZlist)
	 tableW <- lapply(1:K, function(k){tabulate(textlist2[which(finalZlist2==k)], nbins=W)})

	 for (i2 in 1:n2){	 	
	 	for (d in 1:nrow(edge)){ 
	 		currentCd <- currentC[d]
	 		textlistd <-textlist[[d]]
	 		topicpartd <- topicpartZ(currentCd, zcnew, alpha, mvec)
	 		wordpartd <- wordpartZ(K, textlistd, tableW, delta, nvec)
	 		for (w in 1:nrow(wordpartd)){
	 			const2 <- topicpartd + wordpartd[w,]
	    		zwold <- currentZ[[d]][currentCd,w]
	    		zwnew <- selected(1, exp(const2))
	    		if (zwnew!=zwold){
	    			zcnew[[currentCd]][zwold] <- zcnew[[currentCd]][zwold]-1
	    			zcnew[[currentCd]][zwnew] <- zcnew[[currentCd]][zwnew]+1
	    			tableW[[zwold]][textlistd[w]]<-tableW[[zwold]][textlistd[w]]-1
	    			tableW[[zwnew]][textlistd[w]]<-tableW[[zwnew]][textlistd[w]]+1
	    			topicpartd <- topicpartZ(currentCd, zcnew, alpha, mvec)
	 				wordpartd <- wordpartZ(K, textlistd, tableW, delta, nvec)
	    			currentZ[[d]][currentCd,w] <- zwnew	    			
	    			}
	    		}
	    	}
	    }

	logWZmat <-c(logWZmat, logWZ(nIP, K, textlist2, W, zcnew, tableW, alpha, mvec, delta, nvec))
	
	cat("inner iteration 3", "\n")
	bold <- list()
    allxmatlist2<-list()
    lambdaiold <-list()
	for (IP in 1:nIP){
		bold[[IP]] <- rowMeans(bmat[[IP]])
	    edgeC[[IP]]<- edge[which(currentC==IP),]
	    allxmatlist2[[IP]]<- list()
			
		for (d in 1:nrow(edgeC[[IP]])){
			histlist2 <- history(edgeC[[IP]], node, edgeC[[IP]][d,3])
	      	allxmatlist2[[IP]][[d]]<-allxmat(edgeC[[IP]], node, histlist2, edgeC[[IP]][d,1], 0.05)
	    }
	    lambdaiold[[IP]]<-t(vapply(1:nrow(edgeC[[IP]]), function(d){
	      					  multiplyXB(allxmatlist2[[IP]][[d]], bold[[IP]])
	      					  }, rep(1, length(node)-1)))
	}
	
	lambdainew <-list()	
	bnew <- list()		

	for (i3 in 1:n3){ 
		for (IP in 1:nIP){
	    	bnew[[IP]]<-rmvnorm(1, bold[[IP]], (delta_B)^2*diag(P))
	      	lambdainew[[IP]]<-t(vapply(1:nrow(edgeC[[IP]]), function(d){
	      					  multiplyXB(allxmatlist2[[IP]][[d]], bnew[[IP]])
	      					  }, rep(1, length(node)-1)))
	    }
	      	      		      			
		prior <- vapply(1:nIP, function(IP){
				 dmvnorm(bnew[[IP]],rep(0,P), sigma^2*diag(P), log=TRUE)-dmvnorm(bold[[IP]],rep(0,P), sigma^2*diag(P), log=TRUE)
				 }, c(1))
		post <- betapartB(nIP, lambdainew, edgeC)-betapartB(nIP, lambdaiold, edgeC)	  
		loglikediff <- prior+post
		
		u <- log(runif(nIP, 0, 1))
		for (IP in which((u<loglikediff)==TRUE)){
				bold[[IP]] <- bnew[[IP]] 
				lambdaiold[[IP]]<-t(vapply(1:nrow(edgeC[[IP]]), function(d){
	      					  multiplyXB(allxmatlist2[[IP]][[d]], bold[[IP]])
	      					  }, rep(1, length(node)-1)))
		}
		 
		if (i3 > burn && i3%%(thin)==0){
		  for (IP in 1:nIP){
		  	bmat[[IP]][,(i3-burn)/thin]<-bold[[IP]]
		  }
		}
	}
}
	
	if (opt==TRUE){
		out = list(C=currentC, Z=finalZ(currentC, currentZ), B=bmat, L=logWZmat, A=alphamat, M=mvecmat)
		}
	else {
		out = list(C=currentC, Z=finalZ(currentC, currentZ), B=bmat, L=logWZmat)
	}	
	return(out)
}
unix.time(DareMCMC<- MCMC(edge, node, textlist, vocabulary, 3, 10, 0.02, outer=150, n1=1, n2=1, n3=3300, burn=300, thin=3, opt=FALSE, seed=3))
save(DareMCMC, file="DareMCMC_sym3.RData")

