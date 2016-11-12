library(mvtnorm)
library(gtools)
library(MCMCpack)
library(MCMCglmm)
library(reshape2)
library(entropy)
library(tmvtnorm)
library(Rcpp)
library(inline)

###################################################
#load the originals
setwd("/Users/bomin8319/Desktop/SU16/Rcode/IPTMcode")
load("Temporal_Email_Data.Rdata")
Dare<-Temporal_Email_Data$Dare				 
Daretext <-load("Dare_Internal.Rdata")

#pre-processing document_edge_matrix -> duplicate multicast edges
attach(Dare)
Dare_edge<-matrix(NA, nrow= sum(colSums(email_data[,-1:-2])), ncol=3)
iter=1
for (i in 1:nrow(email_data)){
	for (j in 3:ncol(email_data)){
		if (email_data[i,j]==1){Dare_edge[iter,]=c(email_data[i,1], email_data[i,2], j-2); iter=iter+1}
	}
}
colnames(Dare_edge)<-c("Timestamp", "Sender_ID", "Receiver_ID")
Dare_node<-matrix(NA, nrow= length(manager_gender), ncol=3)
Dare_node[,1]<-1:length(manager_gender); Dare_node[,2]<-manager_gender; Dare_node[,3]<-manager_department

time<-strptime(Dare_edge[,1],  "%d %b %Y %H:%M:%S")
unix_time<-as.numeric(as.POSIXct(time, format="%d %b %Y"), origin="1970-01-01")
Dare_edge<-cbind(Dare_edge,unix_time)

edge<-as.data.frame(Dare_edge[,-1])
node<-as.data.frame(Dare_node)					#node info (in case if we want to use this e.g. gender)	
A<- 1:27										#total actors over the corpus


edge<-matrix(as.numeric(as.matrix(edge)), nrow=nrow(edge), ncol=ncol(edge)) 		

#Duplicating the text information for multicast emails
Dare_text<-matrix(NA, nrow=nrow(edge), ncol=ncol(document_word_matrix))
iter=1
	for (i in 1:nrow(document_edge_matrix)){
	for (j in 2:ncol(document_edge_matrix)){
		if (document_edge_matrix[i,j]==1){Dare_text[iter,]=document_word_matrix[i,]; iter=iter+1}
	}
}
text<- Dare_text									
#vocabulary 	 				

#delete two documents with zero word : 1 and 81 -> total 269 documents between 18 actors, 620 vocabulary
edge <- edge[-which(rowSums(text)==0),]
text <- text[-which(rowSums(text)==0),]
text <- text[order(edge[,3]),]			#final document_word_matrix				
edge <- edge[order(edge[,3]),] 			#final document_edge_matrix 

#also need to preprocess the observed word matrix W
textlist <- list()
for (d in 1:nrow(text)){
	words<- which(!text[d,]==0)
	textlist[[d]]<-unlist(sapply(1:length(words), function(w){rep(words[w], text[d,words[w]])}))
}

#for reference, we can also look at the actual words! (but will use index (textlist) for inference)
wordlist<-list()
for (d in 1:nrow(text)){
	wordlist[[d]]<-vocabulary[textlist[[d]]]
	}
	
options(digits=8)
edge[,3]<-edge[,3]/(3600)
edge<-cbind(edge,1:nrow(edge))

#pre-processing over: will use 'edge=(i,j,t)', 'textlist=Z', 'wordlist=W', 'A' (and 'node' for later)

#############################
#Initial assignment of the parameters as in Section 2.2. 
#initialize C, B, and Z 
set.seed(1)                 #change seed for different runs

K=10; 						#assume 20 number of topics
C=3							#assume 5 number of interaction patterns
		
alpha<- 50/K				#Dirichlet concentration prior for theta
mvec <- rep(1, K)/K			#Dirichlet base prior for theta - assign equal probability to every topic
	
P=4							#dimension of x(i,j,t) - we use 4 dynamic network statistics (intercept, send, receive, triangle)
	
sigma=1						#variance of Normal prior for beta

delta=0.01					#Dirichlet concentration prior for phi
W = length(vocabulary)		#total number of words 
nvec <- rep(1, W)/W			#Dirichlet base prior for phi - assign equal probability to every word

eta = 50/C					#Dirichlet concentration prior for gamma
lvec <- rep(1, C)/C			#Dirichlet base prior for gamma - assign equal probability to every interaction pattern
gamma <- rdirichlet(1,  eta*lvec)	#multinomial prior for gamma

#initialize beta and theta 
betainit<-list()
theta <- list()
for (c in 1:C){
	betainit[[c]]<-rmvnorm(1, rep(0, P), sigma^2*diag(P))		#beta follows multivariate normal
	theta[[c]]<- rdirichlet(1,alpha*mvec)						#theta follows dirichlet (alpha, m)
}

#initialize C 
cinit<- vapply(1:nrow(edge), function(d){which(rmultinom(1, 1, gamma)==TRUE)}, c(1))

#initialize Z - based on the initial assignment of interactin pattern C
zinit<-list()
for (d in 1:nrow(text)){
 	Md <- sum(text[d,])				#assign topic for each word in the document
 	zinit[[d]]<-matrix(0, nrow=C, ncol=Md)
 	for(c in 1:C){
 	zinit[[d]][c,]<-vapply(1:Md, function(m){which(rmultinom(1, 1, theta[[c]])==TRUE)}, c(1))		#topic should be assigned jointly with c
}}

################################################
#functions used in MCMC 
#most of the functions written in Rcpp -> refer to sourceCpp file

sourceCpp("/Users/bomin8319/Desktop/SU16/Rcode/IPTMcode/samplercpp.cpp")

#log-partial likelihood used in Section 5.2 of the draft (corresponding to product of Eq.13 in Section 3.2)
#same function written in Rcpp -> same function name with C (much faster but still needs to be faster)
logPL <- function(finalZ, zcnew, currentC, tableW, corpusW){
	sum(vapply(1:nrow(edge), function(d){sum(sapply(1:length(finalZ[[d]]), function(w){log(zcnew[[currentC[d]]][finalZ[[d]][w]]-ifelse(zcnew[[currentC[d]]][finalZ[[d]][w]]>0, 1, 0)+alpha*mvec[finalZ[[d]][w]])-log(length(finalZ[[d]])-1+alpha)+log(tableW[[finalZ[[d]][w]]][textlist[[d]][w]]-ifelse(tableW[[finalZ[[d]][w]]][textlist[[d]][w]] > 0, 1, 0)+delta*nvec[w])-log(length(corpusW[[finalZ[[d]][w]]])-length(unique(corpusW[[finalZ[[d]][w]]]))+delta) 	}))
}, c(1)))
	}

#hyperparameter optimization of topic distribution theta as in Section 5.1.		
#uses rcpp functions currentZ2, Supdate and Skupdate
#Note that can get warning messages as before when the number of counts in topic k is zero -> not a problem
#Warning messages:
#1: In cktable[[k]][cklist[[k]]] <- tabulate(currentC) :
#  number of items to replace is not a multiple of replacement length

vec <- alpha*mvec
parupdateZ <- function(currentZ, currentC){
	corpusCnew <- list()
	zcnew <- list()
	clist <- c()
	for (c in 1:C){
	corpusCnew[[c]] <- currentZ2(c, which(currentC==c), currentZ)
	zcnew[[c]] <- tabulate(unlist(corpusCnew[[c]]), nbins=K)
	clist[c]<-	sum(zcnew[[c]])
	}
	iter=1
	ctable <- rep(0, max(clist))
	ctable[clist]<-tabulate(currentC)
	cklist<-list();
	cktable<-list();
	for (k in 1:K){
	cklist[[k]] <- vapply(1:C, function(c){zcnew[[c]][k]}, c(1))
	cktable[[k]] <-rep(0, max(cklist[[k]]))
   	cktable[[k]][cklist[[k]]] <- tabulate(currentC)}

	while ((abs(alpha-sum(vec))>0.001) | (iter==1)){
	alpha <- sum(vec)
	S <- Supdate(alpha, ctable)		
	s <- Skupdate(vec, cktable)	
	vec<- vec*s/S
	iter=iter+1
	}
	return(vec)	
	}
	
#saving all the previous interaction histories between all actors, before time t -> used for beta part
#could be possilby re-written using Rcpp, since currently quite slow.
x1 <- function(edge, t){
	cauch<- lapply(1:length(A), function(v){lapply(1:length(A), function(u){numeric(0)})})
	edge2 <- matrix(edge[edge[,3] < t,1:3], ncol=3)
	if (nrow(edge2)>0){
	for (d in 1:nrow(edge2)){
		cauch[[edge2[d,1]]][[edge2[d,2]]]<-c(cauch[[edge2[d,1]]][[edge2[d,2]]],t-(edge2[d ,3]))
	}
	}		
return(cauch)
}

###########################################################
#run MCMC and save the results

MCMC =function(outer, n1, n2, n3, delta_B){
 	  currentCmat <- matrix(NA, nrow=outer, ncol=nrow(edge))
 	  currentC <- cinit	
	  currentZ <- zinit
	  alphamat <- alpha
	  mvecmat <- mvec/length(mvec)
	  bmat <- list()
	  logPLmat <- c()								
	  for (c in 1:C){bmat[[c]] <- matrix(NA, nrow=P, ncol=n3); bmat[[c]][,]<-betainit[[c]]}
	  #entropy <- matrix(NA, nrow=outer, ncol=n1)    #Let's skip entropy at this point since it seems to be not helpful
	  
	  for(o in 1:outer){	
	  	 cat("outer iteration = ", o, "\n") 
	  	 if (o%%1==0) {
	  	 vec <- parupdateZ(currentZ, currentC)
	  	 alpha <- sum(vec)
	  	 alphamat <- c(alphamat, alpha)
	  	 mvec <- vec/alpha	
	  	 mvecmat <- rbind(mvecmat, mvec)
	  	 }						
	  	 
	  	 edgeC <-list()
	  	 corpusC<-list()
	  	  zc<-list()
		  allxmat<-list()
		  lambdai<-list()
		
		  cat("inner iteration 1", "\n") 
	  	 for(i in 1:n1){ #inter iteration 1
	      #C update given Z and B - within each document d 	
	       for (d in 1:nrow(edge)){
	       	currentC2<-currentC[-d]
	       	for (c in 1:C){
	       	edgeC[[c]]<- edge[which(currentC2==c),]	
	       	corpusC[[c]]<-currentZ2(c, which(currentC2==c), currentZ)

	       	zc[[c]]<- tabulate(unlist(corpusC[[c]]), nbins=K)
	      	edgeC[[c]]<- rbind(edgeC[[c]], edge[d,]); 
	      	zc[[c]]<-zc[[c]]+tabulate(currentZ[[d]][c,], nbins=K)
	      	pre <- x1(edgeC[[c]], edge[d,3])
	      	allxmat[[c]]<-allxmatC(edgeC[[c]], pre, edge[d,1], A)
	      	lambdai[[c]]<- rowSums(sweep(allxmat[[c]], 2, bmat[[c]][, n3], '*'))
	      	}
	      				      	 
	       	const<-vapply(1:C, function(c){log(gamma[c])+lambdai[[c]][edge[d,2]]-log(sum(exp(lambdai[[c]])))-length(textlist[[d]])*log(sum(zc[[c]])-K+alpha)+sum(vapply(currentZ[[d]][c,], function(k){log(zc[[c]][k]-ifelse(zc[[c]][k]>0, 1, 0)+alpha*mvec[k])}, c(1)))}, c(1))
	        currentC[d]<-which(rmultinom(1, 1, exp(const))==TRUE)
	       }
	      #entropy[o,i]<- entropy.empirical(tabulate(currentC)) 
	      }
	      currentCmat[o,]<-currentC
	      
	      corpusCnew<-list()
	      corpusW<-list()
	      tableW <- list()	
	      zcnew <- list()       
	      cat("inner iteration 2", "\n")
	      for (i in 1:n2){	#inner iteration 2
	      for (c in 1:C){
		  corpusCnew[[c]]<-currentZ2(c, which(currentC==c), currentZ)
	      zcnew[[c]]<-tabulate(unlist(corpusCnew[[c]]), nbins=K) }
	      finalZ <- finalZ2(currentC, currentZ)
	      for (k in 1:K){corpusW[[k]] <- unlist(textlist)[which(unlist(finalZ)==k)]
		  tableW[[k]]<-tabulate(corpusW[[k]], nbins=W)}	
	      for (d in 1:nrow(edge)){
	      for (w in 1:length(textlist[[d]])){	     		
	      const2<- vapply(1:K, function(k){log(zcnew[[currentC[d]]][k]-ifelse(zcnew[[currentC[d]]][k]>0, 1, 0)+alpha*mvec[k])+log(tableW[[k]][textlist[[d]][w]]-ifelse(tableW[[k]][textlist[[d]][w]]>0, 1, 0)+delta*nvec[w])-log(length(corpusW[[k]])-length(unique(corpusW[[k]]))+delta)}, c(1)) 						
	      currentZ[[d]][currentC[d],w] <- which(rmultinom(1, 1, exp(const2))==TRUE)										
	     }
         }
	     }
	  logPLmat <-c(logPLmat, logPLC(finalZ2(currentC, currentZ), zcnew, currentC, tableW, corpusW, textlist, alpha, mvec, delta, nvec))
    
	        bold <- list()
	        allxmat<-list()
	      	for(c in 1:C){
	      	bold[[c]] <- bmat[[c]][, n3]
	      	edgeC[[c]]<- edge[which(currentC==c),]
	      	allxmat[[c]]<- list()
			for (d in 1:nrow(edgeC[[c]])){
			 pre <- x1(edgeC[[c]], edgeC[[c]][d,3])
	      	 allxmat[[c]][[d]]<-allxmatC(edgeC[[c]], pre, edge[d,1], A)
	      	    	}	      	
	      	}
	      	lambdaiold <-list()
	      	lambdainew <-list()	
	      	bnew <- list()		
	        cat("inner iteration 3", "\n")
	      	for (i in 1:n3){ #inner iteration 3
	      	for (c in 1:C){
	      	bnew[[c]]<-rmvnorm(1, bold[[c]], (delta_B)^2*diag(P))
	      	lambdaiold[[c]]<-t(vapply(1:nrow(edgeC[[c]]), function(d){rowSums(sweep(allxmat[[c]][[d]], 2, bold[[c]], '*'))}, rep(1, length(A))))
	      	lambdainew[[c]]<-t(vapply(1:nrow(edgeC[[c]]), function(d){rowSums(sweep(allxmat[[c]][[d]], 2, bnew[[c]], '*'))}, rep(1, length(A))))
	      	}
	      	      		      			
		  	loglikediff<- vapply(1:C, function(c){log(dmvnorm(bnew[[c]],rep(0,P), sigma^2*diag(P)))-log(dmvnorm(bold[[c]],rep(0,P), sigma^2*diag(P)))+sum(vapply(1:nrow(edgeC[[c]]), function(d){lambdainew[[c]][d, edgeC[[c]][d,2]]-log(sum(exp(lambdainew[[c]][d,])))-lambdaiold[[c]][d, edgeC[[c]][d,2]]+log(sum(exp(lambdaiold[[c]][d,])))}, c(1)))}, c(1))          				
		  	u<-runif(5, 0, 1)
		  	u=log(u);	
		  	for (c in 1:C){if (u[c] < loglikediff[c]){bmat[[c]][,i]<-bnew[[c]]; bold[[c]] <-bnew[[c]] } else {bmat[[c]][,i]<-bold[[c]]}
		  		}
		  		}
		  		}
	return(list(C=currentCmat, Z=sapply(1:nrow(edge), function(d){currentZ[[d]][currentC[d],]}), B=bmat, A=alphamat, M=mvecmat, L=logPLmat))
		  		}
unix.time(DareMCMC<- MCMC(100, 3, 3, 1500, 0.05))    #smaller number of iterations due to computation time

