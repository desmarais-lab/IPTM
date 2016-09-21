install.packages("mvtnorm")
install.packages("gtools")
install.packages("MCMCpack")
install.packages("MCMCglmm")
install.packages("reshape2")
install.packages("entropy")
install.packages("tmvtnorm")
library(mvtnorm)
library(gtools)
library(MCMCpack)
library(MCMCglmm)
library(reshape2)
library(entropy)
library(tmvtnorm)
load("Temporal_Email_Data.Rdata")
Vance<-Temporal_Email_Data$Vance					 
Vancetext <-load("Vance_Internal.Rdata")
attach(Vance)
Vance_edge<-matrix(NA, nrow= sum(colSums(email_data[,-1:-2])), ncol=3)
iter=1
for (i in 1:nrow(email_data)){
	for (j in 3:ncol(email_data)){
		if (email_data[i,j]==1){Vance_edge[iter,]=c(email_data[i,1], email_data[i,2], j-2); iter=iter+1}
	}
}
colnames(Vance_edge)<-c("Timestamp", "Sender_ID", "Receiver_ID")
Vance_node<-matrix(NA, nrow= length(manager_gender), ncol=3)
Vance_node[,1]<-1:length(manager_gender); Vance_node[,2]<-manager_gender; Vance_node[,3]<-manager_department

time<-strptime(Vance_edge[,1],  "%d %b %Y %H:%M:%S")
unix_time<-as.numeric(as.POSIXct(time, format="%d %b %Y"), origin="1970-01-01")
Vance_edge<-cbind(Vance_edge,unix_time)

edge<-as.data.frame(Vance_edge[,-1])
node<-as.data.frame(Vance_node)						
A<- 1:18											


edge<-matrix(as.numeric(as.matrix(edge)), nrow=nrow(edge), ncol=ncol(edge)) 		
Vance_text<-matrix(NA, nrow=nrow(edge), ncol=ncol(document_word_matrix))
iter=1
	for (i in 1:nrow(document_edge_matrix)){
	for (j in 2:ncol(document_edge_matrix)){
		if (document_edge_matrix[i,j]==1){Vance_text[iter,]=document_word_matrix[i,]; iter=iter+1}
	}
}
text<- Vance_text									
				



edge <- edge[c(-1, -81),]
text <- text[c(-1, -81),]

text <- text[order(edge[,3]),]						
edge <- edge[order(edge[,3]),] 			


textlist <- list()
for (d in 1:nrow(text)){
	words<- which(!text[d,]==0)
	textlist[[d]]<-unlist(sapply(1:length(words), function(w){rep(words[w], text[d,words[w]])}))
}

wordlist<-list()
for (d in 1:nrow(text)){
	wordlist[[d]]<-vocabulary[textlist[[d]]]
	}
	

options(digits=8)
edge[,3]<-edge[,3]/(3600)


cauch <- list()
x <- function(edge){
	iter=1
for (t in unique(edge[,3])){  cauch[[iter]]<-list()
	for (i in A){ cauch[[iter]][[i]]<-list()
	for (j in A){
		edgeij <- t-(edge[edge[,1]==i & edge[,2]==j & edge[,3] < t,3])
		cauch[[iter]][[i]][[j]]<- edgeij 
		}
		}
		iter=iter+1
}
return(cauch)
}

cauch <- list()
x1 <- function(edge, t){
	for (i in A){ cauch[[i]]<-list()
	for (j in A){
		edgeij <- t-(edge[edge[,1]==i & edge[,2]==j & edge[,3] < t,3])
		cauch[[i]][[j]]<- edgeij 
		}
		}
return(cauch)
}

stats<- function(edge, pre, i, j, t, lambda=0.05){
	mat <- pre[[which(unique(edge[,3])==t)]]
 	send <- sum(exp(-lambda*mat[[i]][[j]]))
	receive <-sum(exp(-lambda*mat[[j]][[i]]))
	triangle <- sum(sapply(A[!A==i & !A==j], function(h){sum(exp(-lambda*mat[[i]][[h]]))*sum(exp(-lambda*mat[[h]][[j]]))}))
	twosend <-sum(sapply(A[!A==i & !A==j], function(h){sum(exp(-lambda*mat[[i]][[h]]))*sum(exp(-lambda*mat[[h]][[j]]))}))
	tworeceive <- sum(sapply(A[!A==i & !A==j], function(h){sum(exp(-lambda*mat[[h]][[i]]))*sum(exp(-lambda*mat[[j]][[h]]))}))
	sibling <- sum(sapply(A[!A==i & !A==j], function(h){sum(exp(-lambda*mat[[h]][[i]]))*sum(exp(-lambda*mat[[h]][[j]]))}))
	cosibling <- sum(sapply(A[!A==i & !A==j], function(h){sum(exp(-lambda*mat[[i]][[h]]))*sum(exp(-lambda*mat[[j]][[h]]))}))
	return(c(1, send, receive, sum(twosend, tworeceive, sibling, cosibling)))
}

stats1<- function(edge, pre, i, j, t, lambda=0.05){
 	send <- sum(exp(-lambda*pre[[i]][[j]]))
	receive <-sum(exp(-lambda*pre[[j]][[i]]))
	triangle <- sum(sapply(A[!A==i & !A==j], function(h){sum(exp(-lambda*pre[[i]][[h]]))*sum(exp(-lambda*pre[[h]][[j]]))}))
	twosend <-sum(sapply(A[!A==i & !A==j], function(h){sum(exp(-lambda*pre[[i]][[h]]))*sum(exp(-lambda*pre[[h]][[j]]))}))
	tworeceive <- sum(sapply(A[!A==i & !A==j], function(h){sum(exp(-lambda*pre[[h]][[i]]))*sum(exp(-lambda*pre[[j]][[h]]))}))
	sibling <- sum(sapply(A[!A==i & !A==j], function(h){sum(exp(-lambda*pre[[h]][[i]]))*sum(exp(-lambda*pre[[h]][[j]]))}))
	cosibling <- sum(sapply(A[!A==i & !A==j], function(h){sum(exp(-lambda*pre[[i]][[h]]))*sum(exp(-lambda*pre[[j]][[h]]))}))
	return(c(1, send, receive, sum(twosend, tworeceive, sibling, cosibling)))
}




set.seed(100)
K=20; 						
C=3							
		
alpha<- 50/20					
mvec <- rep(1, K)			
	
P=4						
	
sigma=1						

delta=0.01				
W = length(vocabulary)		 
nvec <- rep(1, W)			

eta = 50/C						
lvec <- rep(1, C)			
gamma <- rdirichlet(1,  eta*lvec)	


betainit<-list()
theta <- list()
for (c in 1:C){
	betainit[[c]]<-rmvnorm(1, rep(0, P), sigma^2*diag(P))			
	theta[[c]]<- rdirichlet(1,alpha*mvec)						
}



cinit<- sample(1:C, nrow(edge), replace=TRUE)
zinit<-list()
for (d in 1:nrow(text)){
 	Md <- sum(text[d,])				
 	zinit[[d]]<-matrix(0, nrow=C, ncol=Md)
 	for(c in 1:C){
 	zinit[[d]][c,]<-sapply(1:Md, function(m){which(rmultinom(1, 1, theta[[c]])==TRUE)})		
}
}


linit<- lapply(1:C, function(c){rgamma(3, 1, 5)})

set <- A


MCMC =function(outer, n1, n2, n3, delta_B){
 
 	  #define the output matrix and assign initial values
	  currentC <- cinit	
	  currentZ <- zinit
	  bmat <- list()																			
	  for (c in 1:C){bmat[[c]] <- matrix(NA, nrow=P, ncol=n3); bmat[[c]][,]<-betainit[[c]]}
	  entropy <- matrix(NA, nrow=outer, ncol=n1)
	  
	  for(o in 1:outer){	
	  	 print(o)							
	  	
	  	 for(i in 1:n1){ 
	  	  if(i%%10==0) print(i)	
	      zc<-list()
		  allxmat<-list()
	       for (d in 1:nrow(edge)){
	       	currentC2<-currentC[-d]
	       	edgeC<- lapply(1:C, function(c){edge[which(currentC2==c),]})	
	       	corpusC<- lapply(1:C, function(c){sapply(which(currentC2==c), function(d){currentZ[[d]][c,]})})	
	       	for (c in 1:C){
	       	if(length(corpusC[[c]])>0){zc[[c]]<-tabulate(unlist(corpusC[[c]]), nbins=K)} else {zc[[c]] <- rep(0, K)}
	      	edgeC[[c]]<- rbind(edgeC[[c]], edge[d,]); edgeC[[c]]<- edgeC[[c]][order(edgeC[[c]][,3]),]
	      	pre<-x1(edgeC[[c]], edge[d,3])
	      	allxmat[[c]]<-t(sapply(set, function(j){stats1(edgeC[[c]], pre, edge[d,1],j, edge[d,3])}))
	      	for (z in currentZ[[d]][c,]){zc[[c]][z]<-zc[[c]][z]+1} 
	      	}
	       	lambdai<- lapply(1:C, function(c){rowSums(sweep(allxmat[[c]], 2, bmat[[c]][, n3], '*'))})	      
	       	const<-sapply(1:C, function(c){log(gamma[c])+lambdai[[c]][edge[d,2]]-log(sum(exp(lambdai[[c]])))-length(currentZ[[d]][c,])*log(sum(zc[[c]])-length(zc[[c]])+alpha)+sum(sapply(currentZ[[d]][c,], function(k){log(zc[[c]][k]-ifelse(zc[[c]][k]>0, 1, 0)+alpha*mvec[k])}))})											 
	       currentC[d]<-which(rmultinom(1, 1, exp(const))==TRUE)
	       }
	       entropy[o,i]<- entropy.empirical(tabulate(currentC)) }
	      	       
	      for (i in 1:n2){	
	      if(i%%10==0) print(i)
	      corpusCnew<-list()
	      zcnew<-list()
	      for (c in 1:C){corpusCnew[[c]] <- sapply(which(currentC==c), function(d){currentZ[[d]][c,]})			
	      if(length(corpusCnew[[c]])>0){zcnew[[c]]<-tabulate(unlist(corpusCnew[[c]]), nbins=K)} else {zcnew[[c]] <- rep(0, K)}	
	       }													     											
	      corpusW <- list()
	      tableW <- list()
	      for (k in 1:K){corpusW[[k]] <- unlist(sapply(1:nrow(edge), function(d){textlist[[d]][which(currentZ[[d]][c,]==k)]}))
	      	tableW[[k]] <-tabulate(corpusW[[k]], nbins=W)}	
	      for (d in 1:nrow(edge)){
	      for (w in 1:length(currentZ[[d]][currentC[d],])){	     																									
	      const2<-sapply(1:K, function(k){(zcnew[[currentC[d]]][k]-ifelse(zcnew[[currentC[d]]][k]>0, 1, 0)+alpha*mvec[k])*(tableW[[k]][textlist[[d]][w]]-ifelse(tableW[[k]][textlist[[d]][w]]>0, 1, 0)+delta*nvec[w])/(length(corpusW[[k]])-length(unique(corpusW[[k]]))+delta)}) 																											
	      currentZ[[d]][currentC[d],w] <- which(rmultinom(1, 1, const2)==TRUE)															
	     }
         }
	     }
	     
	        bold<-lapply(1:C, function(c){bmat[[c]][, n3]})				
	      	edgeC<- lapply(1:C, function(c){edge[which(currentC==c),]})	
			allxmat<-list()
	      	for(c in 1:C){
	      	allxmat[[c]]<- list()
			for (d in 1:nrow(edgeC[[c]])){
			 pre<-x1(edgeC[[c]], edgeC[[c]][d,3])
	      	 allxmat[[c]][[d]]<-t(sapply(set, function(j){stats1(edgeC[[c]],pre, edge[d,1],j, edge[d,3])}))
	      	}	      	
	      	}			
			for (i in 1:n3){
	      	if(i%%1000==0) print(i)	
	      	bnew<-lapply(1:C, function(c){rmvnorm(1, bold[[c]], (delta_B)^2*diag(P))})	   		      	      		      								
	      	lambdaiold <- lapply(1:C, function(c){t(sapply(1:nrow(edgeC[[c]]), function(d){rowSums(sweep(allxmat[[c]][[d]], 2, bold[[c]], '*'))}))}) 	
         	lambdainew<- lapply(1:C, function(c){t(sapply(1:nrow(edgeC[[c]]), function(d){rowSums(sweep(allxmat[[c]][[d]], 2, bnew[[c]], '*'))}))}) 	
		  	loglikediff<- sapply(1:C, function(c){log(dmvnorm(bnew[[c]],rep(0,P), sigma^2*diag(P)))-log(dmvnorm(bold[[c]],rep(0,P), sigma^2*diag(P)))+	
		  				  sum(sapply(1:nrow(edgeC[[c]]), function(d){lambdainew[[c]][d, edgeC[[c]][d,2]]-log(sum(exp(lambdainew[[c]][d,])))-lambdaiold[[c]][d, edgeC[[c]][d,2]]+log(sum(exp(lambdaiold[[c]][d,])))}))}) 
		    u<-log(runif(5)) 		 																										
		   for (c in 1:C){if (u[c] < loglikediff[c]){bmat[[c]][,i]<-bnew[[c]]; bold[[c]] <-bnew[[c]]} else {bmat[[c]][,i]<-bold[[c]]}}				
	      }																													
		  }	
	   return(list(C=currentC, E=entropy, Z=sapply(1:nrow(edge), function(d){currentZ[[d]][currentC[d],]}), B=bmat))
       }


unix.time(VanceMCMC2<- MCMC(200, 10, 10, 3500, 0.5)) #about 250 seconds for 1 outer iteration
VanceMCMC<-VanceMCMC2
save(VanceMCMC, file="VanceMCMCserver.RData")