library(mvtnorm)
library(gtools)
library(MCMCpack)
library(MCMCglmm)
library(reshape2)
setwd("/Users/bomin8319/Desktop/SU16/Rcode")		#set directory
load("Temporal_Email_Data.Rdata")
Vance<-Temporal_Email_Data$Vance					 
Vancetext <-load("Vance_Internal.Rdata")

#pre-processing document_edge_matrix -> duplicate multicast edges and 
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
node<-as.data.frame(Vance_node)					#node info (in case if we want to use this e.g. gender)	
A<- 1:18											#total actors over the corpus


edge<-matrix(as.numeric(as.matrix(edge)), nrow=nrow(edge), ncol=ncol(edge)) 		

#Duplicating the text information for multicast emails
Vance_text<-matrix(NA, nrow=nrow(edge), ncol=ncol(document_word_matrix))
iter=1
	for (i in 1:nrow(document_edge_matrix)){
	for (j in 2:ncol(document_edge_matrix)){
		if (document_edge_matrix[i,j]==1){Vance_text[iter,]=document_word_matrix[i,]; iter=iter+1}
	}
}
text<- Vance_text									
#vocabulary 	 				



#delete two documents with zero word : 1 and 81 -> total 269 documents between 18 actors, 620 vocabulary
edge <- edge[c(-1, -81),]
text <- text[c(-1, -81),]

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
	
# we have four observed data structures 
# edge=(i,j,t), text=Z, textlist=W, vocabulary
options(digits=8)
edge[,3]<-edge[,3]/(3600)

#define a function for dynamic network effects
x<- function(edge, i,j,t, weight=0.05){
intercept <- 1;
send<-sum(dexp(t-(edge[edge[,1]==i & edge[,2]==j & edge[,3] < t, 3]), rate=weight));
receive <-sum(dexp(t-(edge[edge[,1]==j & edge[,2]==i & edge[,3] < t, 3]), rate=weight));
twosend <- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==i & edge[,2]==h & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==h & edge[,2]==j & edge[,3] < t, 3]), rate=weight))}))
tworeceive <- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==h & edge[,2]==i & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==j & edge[,2]==h & edge[,3] < t, 3]), rate=weight))}))
sibling<- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==h & edge[,2]==i & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==h & edge[,2]==j & edge[,3] < t, 3]), rate=weight))}))
cosibling<- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==i & edge[,2]==h & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==j& edge[,2]==h & edge[,3] < t, 3]), rate=weight))}))
return(c(intercept, 20*send, 20*receive, 20*sum(twosend, tworeceive, sibling, cosibling)))}


#initialize C, B, and Z 
set.seed(100)
K=20; 						#assume 20 number of topics
C=3							#assume 5 number of interaction patterns
		
alpha<-5						#Dirichlet concentration prior for theta
mvec <- rep(1, K)			#Dirichlet base prior for theta - assign equal probability to every topic
	
P=4						#dimension of x(i,j,t) - we use 7 dynamic network statistics
	
sigma=1						#variance of Normal prior for beta

delta=5						#Dirichlet concentration prior for phi
W = length(vocabulary)		#total number of words 
nvec <- rep(1, W)			#Dirichlet base prior for phi - assign equal probability to every word

eta = 5						#Dirichlet concentration prior for gamma
lvec <- rep(1, C)			#Dirichlet base prior for gamma - assign equal probability to every interaction pattern
gamma <- rdirichlet(1,  eta*lvec)	#multinomial prior for gamma

#initialize beta and theta 
betainit<-list()
theta <- list()
for (c in 1:C){
	betainit[[c]]<-rmvnorm(1, rep(0, P), sigma^2*diag(P))			#beta follows multivariate normal
	theta[[c]]<- rdirichlet(1,alpha*mvec)						#theta follows dirichlet (alpha, m)
}


#initialize C 
cinit<- sample(1:C, nrow(edge), replace=TRUE)

#initialize Z - based on the initial assignment of interactin pattern C
zinit<-list()
for (d in 1:nrow(text)){
 	Md <- sum(text[d,])				#assign topic for each word in the document
 	zinit[[d]]<-matrix(0, nrow=C, ncol=Md)
 	for(c in 1:C){
 	zinit[[d]][c,]<-sapply(1:Md, function(m){which(rmultinom(1, 1, theta[[c]])==TRUE)})		#topic should be assigned jointly with c
}
}

set <- A
#Inference for C, B, and Z
MCMC =function(outer, n1, n2, n3, delta_B){
 
 	  #define the output matrix and assign initial values
	  currentC <- cinit	
	  currentZ <- zinit
	  bmat <- list()																			
	  for (c in 1:C){bmat[[c]] <- matrix(NA, nrow=P, ncol=n3); bmat[[c]][,]<-betainit[[c]]}

	  for(o in 1:outer){	#outer iteration
	  	 print(o)							
	  	
	  	 for(i in 1:n1){ 
	  	  if(i%%10==0) print(i)	#inter iteration 1
	      #C update given Z and B - within each document d 	
	      zc<-list()
		  allxmat<-list()
	       for (d in 1:nrow(edge)){
	       	currentC2<-currentC[-d]
	       	edgeC<- lapply(1:C, function(c){edge[which(currentC2==c),]})	
	       	corpusC<- lapply(1:C, function(c){sapply(which(currentC2==c), function(d){currentZ[[d]][c,]})})	
	       	for (c in 1:C){
	       	if(length(corpusC[[c]])>0){zc[[c]]<-tabulate(unlist(corpusC[[c]]), nbins=K)} else {zc[[c]] <- rep(0, K)}
	      	edgeC[[c]]<- rbind(edgeC[[c]], edge[d,]); edgeC[[c]]<- edgeC[[c]][order(edgeC[[c]][,3]),]
	      	allxmat[[c]]<-t(sapply(set, function(j){x(edgeC[[c]], edge[d,1] ,j, edge[d,3])}))
	      	for (z in currentZ[[d]][c,]){zc[[c]][z]<-zc[[c]][z]+1} 
	      	}
	       	lambdai<- lapply(1:C, function(c){rowSums(sweep(allxmat[[c]], 2, bmat[[c]][, n3], '*'))})	      
	       	const<-sapply(1:C, function(c){log(gamma[c])+lambdai[[c]][edge[d,2]]-log(sum(exp(lambdai[[c]])))-length(currentZ[[d]][c,])*log(sum(zc[[c]])-1+alpha)+sum(sapply(currentZ[[d]][c,], function(k){log(zc[[c]][k]-1+alpha*mvec[k])}))})											 
	       currentC[d]<-which(rmultinom(1, 1, exp(const))==TRUE)
	       }
	       }
	      
	      for (i in 1:n2){	#inner iteration 2
	      if(i%%10==0) print(i)
	      #Z update given C and B- within each word m in the document d		
	      corpusCnew<-list()
	      zcnew<-list()
	      for (c in 1:C){ corpusCnew[[c]] <- sapply(which(currentC==c), function(d){currentZ[[d]][c,]})			
	      if(length(corpusCnew[[c]])>0){zcnew[[c]]<-tabulate(unlist(corpusCnew[[c]]), nbins=K)} else {zcnew[[c]] <- rep(0, K) }	
	       }													     											
	      corpusW <- list()
	      for (k in 1:K){corpusW[[k]] <- unlist(sapply(1:nrow(edge), function(d){textlist[[d]][which(currentZ[[d]][c,]==k)]}))}	
	      for (d in 1:nrow(edge)){
	      for (w in 1:length(currentZ[[d]][currentC[d],])){	     																									
	      	const2<-sapply(1:K, function(k){(zcnew[[currentC[d]]][k]-1+alpha*mvec[k])*(tabulate(corpusW[[k]], nbins=W)[textlist[[d]][w]]-1+delta*nvec[w])/sum(length(corpusW[[k]])-1+delta)}) 																											
	     	currentZ[[d]][currentC[d],w] <- which(rmultinom(1, 1, const2)==TRUE)															
	     }
         }
	     }
	     
	      #B update given C and Z - within each interaction pattern C
	        bold<-lapply(1:C, function(c){bmat[[c]][, n3]})					#old beta from previous iteration
	      	edgeC<- lapply(1:C, function(c){edge[which(currentC==c),]})	
			allxmat<-list()
	      	for(c in 1:C){
	      	allxmat[[c]]<- list()
			for (d in 1:nrow(edgeC[[c]])){
	      	 allxmat[[c]][[d]]<-t(sapply(set, function(j){x(edgeC[[c]], edgeC[[c]][d,1] ,j, edgeC[[c]][d,3])}))
	      	}	      	
	      	}			
			for (i in 1:n3){
	      	if(i%%1000==0) print(i)	
	      	bnew<-lapply(1:C, function(c){rmvnorm(1, bold[[c]], (delta_B)^2*diag(P))})	   		      	      		      								
	      	lambdaiold <- lapply(1:C, function(c){t(sapply(1:nrow(edgeC[[c]]), function(d){rowSums(sweep(allxmat[[c]][[d]], 2, bold[[c]], '*'))}))}) 	
         	lambdainew<- lapply(1:C, function(c){t(sapply(1:nrow(edgeC[[c]]), function(d){rowSums(sweep(allxmat[[c]][[d]], 2, bnew[[c]], '*'))}))}) 	
		  	loglikediff<- sapply(1:C, function(c){log(dmvnorm(bnew[[c]],rep(0,P), sigma^2*diag(P)))-log(dmvnorm(bold[[c]],rep(0,P), sigma^2*diag(P)))+	
		  				  sum(sapply(1:nrow(edgeC[[c]]), function(d){lambdainew[[c]][d, edgeC[[c]][d,2]]-log(sum(exp(lambdainew[[c]][d,])))-lambdaiold[[c]][d, edgeC[[c]][d,2]]+log(sum(exp(lambdaiold[[c]][d,])))}))})          				
		  	u=log(runif(5)) 																												
		   for (c in 1:C){if (u[c] < loglikediff[c]){bmat[[c]][,i]<-bnew[[c]]; bold[[c]] <-bnew[[c]] } else {bmat[[c]][,i]<-bold[[c]]} }												#M-H update
	      }																													
		  }	
	   return(list(C=currentC, Z=sapply(1:nrow(edge), function(d){currentZ[[d]][currentC[d],]}), B=bmat))
       }


unix.time(VanceMCMC<- MCMC(100, 10, 10, 3500, 0.85)) #about 250 seconds for 1 outer iteration
save(VanceMCMC, file="VanceMCMC.RData")
sapply(1:3, function(c){nrow(unique(t(VanceMCMC$B[[c]])))})/(3500)		#acceptance rate for beta
CofD<- VanceMCMC$C			#table of IP for each document
B<-lapply(1:3, function(c){VanceMCMC$B[[c]][,501:3500]})
Bnew<-list()
for (c in 1:3){
	Bnew[[c]]<-matrix(NA, nrow=4, ncol=1000)
	for(i in 1:1000){Bnew[[c]][,i]<-B[[c]][,3*i]}
}

sapply(1:3, function(c){rowMeans(Bnew[[c]])})
exp(sapply(1:5, function(c){rowMeans(Bnew[[c]])}))

mcmc <- mcmc(t(Bnew[[3]]))
summary(mcmc)
geweke.diag(mcmc)

par(mfrow=c(2,4))
plot(Bnew[[3]][1,], type="l", main="Trace of intercept", xlab="Iterations", ylab="")
plot(Bnew[[3]][2,], type="l", main="Trace of send", xlab="Iterations", ylab="")
plot(Bnew[[3]][3,], type="l", main="Trace of receive", xlab="Iterations", ylab="")
plot(Bnew[[3]][4,], type="l", main="Trace of triangles", xlab="Iterations", ylab="")
hist(Bnew[[3]][1,],freq=FALSE, breaks=15, main="Density of intercept", xlab="beta_0")
hist(Bnew[[3]][2,],freq=FALSE, breaks=15, main="Density of send", xlab="beta_1")
hist(Bnew[[3]][3,],freq=FALSE, breaks=15, main="Density of receive", xlab="beta_2")
hist(Bnew[[3]][4,],freq=FALSE, breaks=15, main="Density of triangles", xlab="beta_3")



data<-list()
for(b in 1:4){
	data[[b]] <- sapply(1:3, function(c){cbind(Bnew[[c]][b,])})}
forbox <-melt(data)
par(xpd=T, mar=par()$mar+c(0,0,0,4))
boxplot<-boxplot(forbox$value~forbox$Var2+forbox$L1, at=c(1,2,3,5,6,7,9,10,11,13,14,15),
col=gray.colors(5), axes=FALSE)
title('Comparison of beta coefficients for IP1-IP3')
axis(2, labels=TRUE)
#axis(1, labels=FALSE)
box()
axis(side=1, at=c(2,6,10,14), labels=c('intercept','send','receive','triangles'), line=0.5, lwd=0 )
legend(locator(1), c("IP1","IP2","IP3"), col=gray.colors(5), pch=15)


#Summary of Z
Zsummary<-VanceMCMC$Z
for (d in 1:269){
	names(Zsummary[[d]])<-wordlist[[d]]
}
Ztable<-tabulate(unlist(Zsummary), nbins=20)														#table of word-topic assignment for each document
names(Ztable)<-c(1:20)
Ztable<-rbind(Ztable[order(Ztable, decreasing=TRUE)], Ztable[order(Ztable, decreasing=TRUE)]/sum(Ztable))

Zsummary1<-list()
iter=1
for (d in which(CofD==1)){
	Zsummary1[[iter]]<-{Zsummary[[d]]}; iter=iter+1
}
#Vance_edge[which(CofD==1), 1]
topicdist1<-tabulate(unlist(Zsummary1), nbins=20)
topicdist1<-topicdist1/sum(topicdist1)
names(topicdist1)<-c(1:20)

Zsummary2<-list()
iter=1
for (d in which(CofD==2)){
	Zsummary2[[iter]]<-{Zsummary[[d]]}; iter=iter+1
}
#Vance_edge[which(CofD==2), 1]
topicdist2<-tabulate(unlist(Zsummary2), nbins=20)
topicdist2<-topicdist2/sum(topicdist2)
names(topicdist2)<-c(1:20)

Zsummary3<-list()
iter=1
for (d in which(CofD==3)){
	Zsummary3[[iter]]<-{Zsummary[[d]]}; iter=iter+1
}
#Vance_edge[which(CofD==3), 1]
topicdist3<-tabulate(unlist(Zsummary3), nbins=20)
topicdist3<-topicdist3/sum(topicdist3)
names(topicdist3)<-c(1:20)

summary <- rbind(topicdist1,topicdist2, topicdist3)
barplot(summary, beside=TRUE)
title('Topic Distritubitions given IP1-IP3')
legend(locator(1), c("IP1","IP2","IP3"), col=gray.colors(5), pch=15)


#word summary
allwords<-unlist(Zsummary)
Wsummary<-sapply(1:20, function(k){allwords[allwords==k]})
Wtable<-sapply(1:20, function(k){table(names(Wsummary[[k]]))})
W6<-Wtable[[6]]/sum(Wtable[[6]])
W6[order(W6, decreasing=TRUE)][1:10]
W9<-Wtable[[9]]/sum(Wtable[[9]])
W9[order(W9, decreasing=TRUE)][1:10]
W12<-Wtable[[12]]/sum(Wtable[[12]])
W12[order(W12, decreasing=TRUE)][1:10]
W13<-Wtable[[13]]/sum(Wtable[[13]])
W13[order(W13, decreasing=TRUE)][1:10]
W1<-Wtable[[1]]/sum(Wtable[[1]])
W1[order(W1, decreasing=TRUE)][1:10]
W14<-Wtable[[14]]/sum(Wtable[[14]])
W14[order(W14, decreasing=TRUE)][1:10]
W10<-Wtable[[10]]/sum(Wtable[[10]])
W10[order(W10, decreasing=TRUE)][1:10]
W5<-Wtable[[5]]/sum(Wtable[[5]])
W5[order(W5, decreasing=TRUE)][1:10]
W19<-Wtable[[19]]/sum(Wtable[[19]])
W19[order(W19, decreasing=TRUE)][1:10]
W20<-Wtable[[20]]/sum(Wtable[[20]])
W20[order(W20, decreasing=TRUE)][1:10]

