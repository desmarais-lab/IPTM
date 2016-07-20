library(mvtnorm)
library(gtools)
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
vocabulary 	 				



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



#define a function for dynamic network effects
x<- function(i,j,t, weight=1){
intercept <- 1;
send<-sum(dexp(t-(edge[edge[,1]==i & edge[,2]==j & edge[,3] < t, 3]), rate=weight));
receive <-sum(dexp(t-(edge[edge[,1]==j & edge[,2]==i & edge[,3] < t, 3]), rate=weight));
twosend <- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==i & edge[,2]==h & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==h & edge[,2]==j & edge[,3] < t, 3]), rate=weight))}))
tworeceive <- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==h & edge[,2]==i & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==j & edge[,2]==h & edge[,3] < t, 3]), rate=weight))}))
sibling<- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==h & edge[,2]==i & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==h & edge[,2]==j & edge[,3] < t, 3]), rate=weight))}))
cosibling<- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==i & edge[,2]==h & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==j& edge[,2]==h & edge[,3] < t, 3]), rate=weight))}))
return(c(intercept, send, receive, twosend, tworeceive, sibling, cosibling))}


#preprocessing for x function in order to save time
#allxmat <- list()
#for(i in set){
#	allxmat[[i]]<- array(0, dim=c(length(edge[,3]),7,length(set)))
#	for (j in set){
#		if(!i==j){
#	allxmat[[i]][,,j]<-t(sapply(edge[,3], function(t){x(i,j,t)}))}
#	else {allxmat[[i]][,,j]<-t(sapply(edge[,3], function(t){c(1, rep(0,6))}))} 	
#	}	
#}

#save(allxmat, file="allxmat.RData")
load("allxmat.RData")						#allxmat[[sender]][d,,receiver] gives us x(sender, receiver, time) of dth document

#rawlambda function provides matrix of stochastic intensities (before exponentiation = beta%*%x(i,j,t)) at time of dth document 
#note that edgetime=d (row number of the document)
rawlambda <- function(beta, A, edgetime){
	lambdaij<-matrix(0, nrow=length(A),ncol=length(A))
	xmat <- list()
	 for(i in 1:length(A)){
	 		xmat[[i]] <- t(sapply(1:length(A), function(j){allxmat[[i]][edgetime,,j]}))
	 	for(j in 1:length(A)){
	 		lambdaij[i,] <- (rowSums(sweep(xmat[[i]], 2, beta, '*')))
	diag(lambdaij)<-0
		}
		}
	return(lambdaij)
}


#initialize C, B, and Z 
set.seed(100)
K=20; 						#assume 20 number of topics
C=5							#assume 5 number of interaction patterns
		
alpha<-5						#Dirichlet concentration prior for theta
mvec <- rep(1, K)			#Dirichlet base prior for theta - assign equal probability to every topic
	
P=7							#dimension of x(i,j,t) - we use 7 dynamic network statistics
	
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


#initialize C - 4th column of edge is the initial interaction pattern for the document C
edge<-cbind(edge, sample(1:C, nrow(edge), replace=TRUE))


#initialize Z - based on the initial assignment of interactin pattern C
zinit<-list()
for (d in 1:nrow(text)){
 	Md <- sum(text[d,])				#assign topic for each word in the document
 	zinit[[d]]<-sapply(1:Md, function(m){which(rmultinom(1, 1, theta[[edge[d,4]]])==TRUE)})		
}


#Inference for C, B, and Z
MCMC =function(n, delta_B){
 
 	  #define the output matrix and assign initial values
	  cmat <- matrix(NA, nrow=nrow(edge), ncol=n+1) 	
	  zmat <- list()
	  bmat <- list()
	  cmat[,1] <-edge[,4]																						
	  for (c in 1:C){bmat[[c]] <- matrix(NA, nrow=P, ncol=n+1); bmat[[c]][,1]<-betainit[[c]]}
	  for(d in 1:nrow(edge)){zmat[[d]]<-matrix(NA, nrow=n+1, ncol=length(zinit[[d]]));zmat[[d]][1,]<-zinit[[d]]}
	  
	  for(i in 1:n){
	      if (i %% 10 == 0) print(i)		
	      
	      #C update - within each document d
	      for (d in 1:nrow(edge)){ 
	      const<-c()																														#for each c, calculate the constant used for multinomial 
	      currentZ<-zmat[[d]][i,]																											#use previous Z
	      for (c in 1:C){							
	      	corpusC <- sapply(which(cmat[,i]==c), function(d){zmat[[d]][i,]})																#extract subset of topic assignments with interaction pattern c
	      if(length(corpusC)>0){zc<-tabulate(unlist(corpusC), nbins=K)} else {zc <- rep(0, K)}												#extract table for topics appear in the document
	        lambdai <- rowSums(sweep(t(sapply(1:length(A), function(j){allxmat[[edge[d,1]]][d,,j]})), 2, bmat[[c]][,i], '*'))				#beta%*%x(i,j,t) given previous beta
	      	const[c]<-log(gamma[c])+lambdai[edge[d,2]]-log(sum(exp(lambdai)))-length(currentZ)*log(sum(zc)-1+alpha)+sum(sapply(currentZ, function(k){log(zc[k]-1+alpha*mvec[k])})) #Eq. (11) 
	      	}
 			cmat[d,i+1]<-which(rmultinom(1, 1, exp(const))==TRUE)																			#Gibbs update for c 
	      }	 
	      
	      #Z update - within each word m in the document d  			
	      corpusCnew<-list()
	      zcnew<-list()
	      for (c in 1:C){ corpusCnew[[c]] <- sapply(which(cmat[,i+1]==c), function(d){zmat[[d]][i,]})										#extract subset of documents with updated c
	      if(length(corpusCnew[[c]])>0){zcnew[[c]]<-tabulate(unlist(corpusCnew[[c]]), nbins=K)} else {zcnew[[c]] <- rep(0, K) }				#extract table for topics appear in the documents (CK)
		  }													     											
	      corpusW <- list()
	      for (k in 1:K){corpusW[[k]] <- unlist(sapply(1:nrow(edge), function(d){textlist[[d]][which(zmat[[d]][i,]==k)]}))}					#extract observed words for each topic k (WK)
	      for (d in 1:nrow(edge)){
	      currentC <- cmat[d, i+1]																											#use the updated c for the document 
	      for (w in 1:length(zmat[[d]][i,])){	     																									
	      	const2<-sapply(1:K, function(k){(zcnew[[currentC]][k]-1+alpha*mvec[k])*(tabulate(corpusW[[k]], nbins=W)[textlist[[d]][w]]-1+delta*nvec[w])/sum(length(corpusW[[k]])-1+delta)}) 																																					#for each topic k, calculate the constants used for multinomial 
	     	zmat[[d]][i+1,w]<- (which(rmultinom(1, 1, const2)==TRUE))																		#Gibbs update for z
         }
         }
	    
	      #B update - within each interaction pattern C (but using all documents in corpus 1:D)
	      for(c in 1:C){	
	      	edgeC <- edge[which(cmat[,i+1]==c),]																								#extract subset of documents with interaction pattern c
	      	bold <- bmat[[c]][,i]																											#old beta from previous iteration
		  	bnew <- rmvnorm(1, bold, (delta_B)^2*diag(P))																					#new beta using rmvnorm
		  	lambdaiold<- t(sapply(1:nrow(edgeC), function(d){rowSums(sweep(t(sapply(1:length(A), function(j){allxmat[[edge[d,1]]][d,,j]})), 2, bold, '*'))}))			#sum of beta%*%x(i,j,t) given bold
         	lambdainew<- t(sapply(1:nrow(edgeC), function(d){rowSums(sweep(t(sapply(1:length(A), function(j){allxmat[[edge[d,1]]][d,,j]})), 2, bnew, '*'))}))			#sum of beta%*%x(i,j,t) given bnew
		  	loglikediff <- log(dmvnorm(bnew,rep(0,P), sigma^2*diag(P)))-log(dmvnorm(bold,rep(0,P), sigma^2*diag(P)))+	
		  				   sum(sapply(1:nrow(edgeC), function(d){lambdainew[d, edgeC[d,2]]}))-sum(sapply(1:nrow(edgeC), function(d){log(sum(exp(lambdainew)))}))-
		  				   sum(sapply(1:nrow(edgeC), function(d){lambdaiold[d, edgeC[d,2]]}))+sum(sapply(1:nrow(edgeC), function(d){log(sum(exp(lambdaiold)))}))			#difference in loglikelihood
	            u=runif(1) 																																				#compare with unifrom(0,1)
		  if (log(u) < loglikediff){bmat[[c]][,i+1]<-bnew } else {bmat[[c]][,i+1]<-bold}														#M-H update
	      }																													
   		  	  }
		  		
	   return(list(C=cmat, Z=zmat, B=bmat))
       }

unix.time(try2 <- MCMC(10^2, 1)) #about 200 seconds for 100 samples

unix.time(try <- MCMC(2*10^4, 1)) #about 10.36 hours for 20000 samples
save(try, file="try.RData")
load("try.RData")

#5000 samples dropped as burn-in
Cnew <- try$C[,-(1:5001)]
Znew<-list()
for(d in 1:269){Znew[[d]] <- as.matrix(try$Z[[d]][-(1:5001),])}
Bnew<-list()
for(c in 1:5){Bnew[[c]] <- try$B[[c]][,-(1:5001)]}

#thinning -> every 15th sample taken so we end up with 1000 sample
iter=1; Cnew2<-matrix(NA, nrow=269, ncol=1000)
for (i in 1:1000){Cnew2[,iter]<- Cnew[,15*i]; iter=iter+1}
Znew2<-list()
for (d in 1:269){Znew2[[d]]<-matrix(NA, nrow=1000, ncol=ncol(Znew[[d]]))
	iter=1;
	for (i in 1:1000){Znew2[[d]][iter,]<- Znew[[d]][15*i,]; iter=iter+1}}
Bnew2<-list()
for (c in 1:5){Bnew2[[c]]<-matrix(NA, nrow=7, ncol=1000)
	iter=1;
	for (i in 1:1000){Bnew2[[c]][,iter]<- Bnew[[c]][,15*i]; iter=iter+1}}

#Summary of result
sapply(1:5, function(c){nrow(unique(t(try$B[[c]])))})/(2*10^4)		#acceptance rate for beta
CofD<-sapply(1:269, function(d){Cnew2[d,1000]})				#table of IP for each document
data<-list()
for(b in 1:7){
	data[[b]] <- sapply(1:5, function(c){cbind(Bnew2[[c]][b,])})}
forbox <-melt(data)
boxplot<-boxplot(forbox$value~forbox$Var2+forbox$L1, at=c(1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,19,20,21,22,23,25,26,27,28,29,31,32,33,34,35,37,38,39,40,41),
col=gray.colors(5), axes=FALSE)
title('Comparison of beta coefficients for IP1-IP5')
axis(2, labels=TRUE)
#axis(1, labels=FALSE)
box()
axis(side=1, at=c(3,9,15,21,27,33,39), labels=c('intercept','send','receive','2-send','2-recieve','sibling', 'co-sibling'), line=0.5, lwd=0 )

#posterior mean for beta 
exp(sapply(1:5, function(c){colMeans(t(Bnew2[[c]]))}))				#Take exponential if we want to interpret this

Zsummary<-list()
for (d in 1:269){
	Zsummary[[d]]<-Znew2[[d]][1000,]									#summarize using the last iteration
	names(Zsummary[[d]])<-c(wordlist[[d]])
}
Zsummary																#table of word-topic assignment for each document

Zsummary1<-list()
iter=1
for (d in which(CofD==1)){
	Zsummary1[[iter]]<-{Zsummary[[d]]}; iter=iter+1
}
#Vance_edge[which(CofD==1), 1]
topicdist1<-tabulate(unlist(Zsummary1), nbins=20)
names(topicdist1)<-c(1:20)
topicdist1<-topicdist1[order(topicdist1, decreasing=TRUE)]/sum(topicdist1)
Zsum<- unlist(Zsummary1)
Zsum1 <- matrix(0, nrow=20, ncol=length(unique(names(Zsum)))); colnames(Zsum1)<-unique(names(Zsum))
for (i in 1:ncol(Zsum1)){
	listofi<-tabulate(Zsum[which(names(Zsum)==colnames(Zsum1)[i])], nbins=20)
	Zsum1[,i]<-listofi}
wordgiventopic15<-Zsum1[15,]/sum(Zsum1[15,])
wordgiventopic15<-wordgiventopic15[order(wordgiventopic15, decreasing=TRUE)]
wordgiventopic4<-Zsum1[4,]/sum(Zsum1[4,])
wordgiventopic4<-wordgiventopic4[order(wordgiventopic4, decreasing=TRUE)]
wordgiventopic8<-Zsum1[8,]/sum(Zsum1[8,])
wordgiventopic8<-wordgiventopic8[order(wordgiventopic8, decreasing=TRUE)]
wordgiventopic10<-Zsum1[10,]/sum(Zsum1[10,])
wordgiventopic10<-wordgiventopic10[order(wordgiventopic10, decreasing=TRUE)]
wordgiventopic15[1:10];wordgiventopic4[1:10];wordgiventopic8[1:10];wordgiventopic10[1:10]


Zsummary2<-list()
iter=1
for (d in which(CofD==2)){
	Zsummary2[[iter]]<-{Zsummary[[d]]}; iter=iter+1
}
#Vance_edge[which(CofD==2), 1]
topicdist2<-tabulate(unlist(Zsummary2), nbins=20)
names(topicdist2)<-c(1:20)
topicdist2<-topicdist2[order(topicdist2, decreasing=TRUE)]/sum(topicdist2)
Zsum<- unlist(Zsummary2)
Zsum2 <- matrix(0, nrow=20, ncol=length(unique(names(Zsum)))); colnames(Zsum2)<-unique(names(Zsum))
for (i in 1:ncol(Zsum2)){
	listofi<-tabulate(Zsum[which(names(Zsum)==colnames(Zsum2)[i])], nbins=20)
	Zsum2[,i]<-listofi}
wordgiventopic18<-Zsum2[18,]/sum(Zsum2[18,])
wordgiventopic18<-wordgiventopic18[order(wordgiventopic18, decreasing=TRUE)]
wordgiventopic19<-Zsum2[19,]/sum(Zsum2[19,])
wordgiventopic19<-wordgiventopic19[order(wordgiventopic19, decreasing=TRUE)]
wordgiventopic3<-Zsum2[3,]/sum(Zsum2[3,])
wordgiventopic3<-wordgiventopic3[order(wordgiventopic3, decreasing=TRUE)]
wordgiventopic13<-Zsum2[13,]/sum(Zsum2[13,])
wordgiventopic13<-wordgiventopic13[order(wordgiventopic13, decreasing=TRUE)]
wordgiventopic18[1:10];wordgiventopic19[1:10];wordgiventopic3[1:10];wordgiventopic13[1:10]


Zsummary3<-list()
iter=1
for (d in which(CofD==3)){
	Zsummary3[[iter]]<-{Zsummary[[d]]}; iter=iter+1
}
#Vance_edge[which(CofD==3), 1]
topicdist3<-tabulate(unlist(Zsummary3), nbins=20)
names(topicdist3)<-c(1:20)
topicdist3<-topicdist3[order(topicdist3, decreasing=TRUE)]/sum(topicdist3)
Zsum<- unlist(Zsummary3)
Zsum3 <- matrix(0, nrow=20, ncol=length(unique(names(Zsum)))); colnames(Zsum3)<-unique(names(Zsum))
for (i in 1:ncol(Zsum3)){
	listofi<-tabulate(Zsum[which(names(Zsum)==colnames(Zsum3)[i])], nbins=20)
	Zsum3[,i]<-listofi}
wordgiventopic9<-Zsum3[9,]/sum(Zsum3[9,])
wordgiventopic9<-wordgiventopic9[order(wordgiventopic9, decreasing=TRUE)]
wordgiventopic7<-Zsum3[7,]/sum(Zsum3[7,])
wordgiventopic7<-wordgiventopic7[order(wordgiventopic7, decreasing=TRUE)]
wordgiventopic11<-Zsum3[11,]/sum(Zsum3[11,])
wordgiventopic11<-wordgiventopic11[order(wordgiventopic3, decreasing=TRUE)]
wordgiventopic5<-Zsum3[5,]/sum(Zsum3[5,])
wordgiventopic5<-wordgiventopic5[order(wordgiventopic5, decreasing=TRUE)]
wordgiventopic9[1:10];wordgiventopic7[1:10];wordgiventopic11[1:10];wordgiventopic5[1:10]


Zsummary4<-list()
iter=1
for (d in which(CofD==4)){
	Zsummary4[[iter]]<-{Zsummary[[d]]}; iter=iter+1
}
#Vance_edge[which(CofD==4), 1]
topicdist4<-tabulate(unlist(Zsummary4), nbins=20)
names(topicdist4)<-c(1:20)
topicdist4<-topicdist4[order(topicdist4, decreasing=TRUE)]/sum(topicdist4)
Zsum<- unlist(Zsummary4)
Zsum4 <- matrix(0, nrow=20, ncol=length(unique(names(Zsum)))); colnames(Zsum4)<-unique(names(Zsum))
for (i in 1:ncol(Zsum4)){
	listofi<-tabulate(Zsum[which(names(Zsum)==colnames(Zsum4)[i])], nbins=20)
	Zsum4[,i]<-listofi}
wordgiventopic14<-Zsum4[14,]/sum(Zsum4[14,])
wordgiventopic14<-wordgiventopic14[order(wordgiventopic14, decreasing=TRUE)]
wordgiventopic12<-Zsum4[12,]/sum(Zsum4[12,])
wordgiventopic12<-wordgiventopic12[order(wordgiventopic12, decreasing=TRUE)]
wordgiventopic17<-Zsum4[17,]/sum(Zsum4[17,])
wordgiventopic17<-wordgiventopic17[order(wordgiventopic17, decreasing=TRUE)]
wordgiventopic1<-Zsum4[1,]/sum(Zsum4[1,])
wordgiventopic1<-wordgiventopic1[order(wordgiventopic1, decreasing=TRUE)]
wordgiventopic14[1:10];wordgiventopic12[1:10];wordgiventopic17[1:10];wordgiventopic1[1:10]


Zsummary5<-list()
iter=1
for (d in which(CofD==5)){
	Zsummary5[[iter]]<-{Zsummary[[d]]}; iter=iter+1
}
#Vance_edge[which(CofD==5), 1]
topicdist5<-tabulate(unlist(Zsummary5), nbins=20)
names(topicdist5)<-c(1:20)
topicdist5<-topicdist5[order(topicdist5, decreasing=TRUE)]/sum(topicdist5)
Zsum<- unlist(Zsummary5)
Zsum5 <- matrix(0, nrow=20, ncol=length(unique(names(Zsum)))); colnames(Zsum5)<-unique(names(Zsum))
for (i in 1:ncol(Zsum5)){
	listofi<-tabulate(Zsum[which(names(Zsum)==colnames(Zsum5)[i])], nbins=20)
	Zsum5[,i]<-listofi}
wordgiventopic1<-Zsum5[1,]/sum(Zsum5[1,])
wordgiventopic1<-wordgiventopic1[order(wordgiventopic1, decreasing=TRUE)]
wordgiventopic16<-Zsum5[16,]/sum(Zsum5[16,])
wordgiventopic16<-wordgiventopic16[order(wordgiventopic16, decreasing=TRUE)]
wordgiventopic2<-Zsum5[2,]/sum(Zsum5[2,])
wordgiventopic2<-wordgiventopic2[order(wordgiventopic2, decreasing=TRUE)]
wordgiventopic20<-Zsum5[20,]/sum(Zsum5[20,])
wordgiventopic20<-wordgiventopic20[order(wordgiventopic20, decreasing=TRUE)]
wordgiventopic1[1:10];wordgiventopic16[1:10];wordgiventopic2[1:10];wordgiventopic20[1:10]




###################################################################
Dare<-Temporal_Email_Data$Dare					 
Daretext <-load("Dare_Internal.Rdata")

#pre-processing document_edge_matrix -> duplicate multicast edges and 
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
A<- 1:27											#total actors over the corpus


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
vocabulary 	 				



#delete two documents with zero word -> total 4845 documents between 26 actors (actor 17 isolate), 2907 vocabulary
nwords<-sapply(1:nrow(text), function(d){sum(text[d,])})
edge <- edge[-which(nwords==0),]			#final document_word_matrix	
text <- text[-which(nwords==0),]			#final document_edge_matrix 

		

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



#define a function for dynamic network effects
x<- function(i,j,t, weight=1){
intercept <- 1;
send<-sum(dexp(t-(edge[edge[,1]==i & edge[,2]==j & edge[,3] < t, 3]), rate=weight));
receive <-sum(dexp(t-(edge[edge[,1]==j & edge[,2]==i & edge[,3] < t, 3]), rate=weight));
twosend <- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==i & edge[,2]==h & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==h & edge[,2]==j & edge[,3] < t, 3]), rate=weight))}))
tworeceive <- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==h & edge[,2]==i & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==j & edge[,2]==h & edge[,3] < t, 3]), rate=weight))}))
sibling<- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==h & edge[,2]==i & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==h & edge[,2]==j & edge[,3] < t, 3]), rate=weight))}))
cosibling<- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==i & edge[,2]==h & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==j& edge[,2]==h & edge[,3] < t, 3]), rate=weight))}))
return(c(intercept, send, receive, twosend, tworeceive, sibling, cosibling))}


#preprocessing for x function in order to save time
allxmat <- allxmat2
for(i in 14:27){
	print(i)
	allxmat[[i]]<- array(0, dim=c(length(edge[,3]),7,length(A)))
	for (j in A){ print(j)
		if(!i==j){
	allxmat[[i]][,,j]<-t(sapply(edge[,3], function(t){x(i,j,t)}))}
	else {allxmat[[i]][,,j]<-t(sapply(edge[,3], function(t){c(1, rep(0,6))}))} 	
	}	
}
allxmat2 <-allxmat
save(allxmat2, file="allxmat2.RData")
load("allxmat2.RData")						#allxmat[[sender]][d,,receiver] gives us x(sender, receiver, time) of dth document

#rawlambda function provides matrix of stochastic intensities (before exponentiation = beta%*%x(i,j,t)) at time of dth document 
#note that edgetime=d (row number of the document)
rawlambda <- function(beta, A, edgetime){
	lambdaij<-matrix(0, nrow=length(A),ncol=length(A))
	xmat <- list()
	 for(i in 1:length(A)){
	 		xmat[[i]] <- t(sapply(1:length(A), function(j){allxmat[[i]][edgetime,,j]}))
	 	for(j in 1:length(A)){
	 		lambdaij[i,] <- (rowSums(sweep(xmat[[i]], 2, beta, '*')))
	diag(lambdaij)<-0
		}
		}
	return(lambdaij)
}


#initialize C, B, and Z 
set.seed(100)
K=20; 						#assume 20 number of topics
C=5							#assume 5 number of interaction patterns
		
alpha<-5						#Dirichlet concentration prior for theta
mvec <- rep(1, K)			#Dirichlet base prior for theta - assign equal probability to every topic
	
P=7							#dimension of x(i,j,t) - we use 7 dynamic network statistics
	
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


#initialize C - 4th column of edge is the initial interaction pattern for the document C
edge<-cbind(edge, sample(1:C, nrow(edge), replace=TRUE))


#initialize Z - based on the initial assignment of interactin pattern C
zinit<-list()
for (d in 1:nrow(text)){
 	Md <- sum(text[d,])				#assign topic for each word in the document
 	zinit[[d]]<-sapply(1:Md, function(m){which(rmultinom(1, 1, theta[[edge[d,4]]])==TRUE)})		
}


#Inference for C, B, and Z
MCMC =function(n, delta_B){
 
 	  #define the output matrix and assign initial values
	  cmat <- matrix(NA, nrow=nrow(edge), ncol=n+1) 	
	  zmat <- list()
	  bmat <- list()
	  cmat[,1] <-edge[,4]																						
	  for (c in 1:C){bmat[[c]] <- matrix(NA, nrow=P, ncol=n+1); bmat[[c]][,1]<-betainit[[c]]}
	  for(d in 1:nrow(edge)){zmat[[d]]<-matrix(NA, nrow=n+1, ncol=length(zinit[[d]]));zmat[[d]][1,]<-zinit[[d]]}
	  
	  for(i in 1:n){
	      if (i %% 10 == 0) print(i)		
	      
	      #C update - within each document d
	      for (d in 1:nrow(edge)){ 
	      const<-c()																														#for each c, calculate the constant used for multinomial 
	      currentZ<-zmat[[d]][i,]																											#use previous Z
	      for (c in 1:C){							
	      	corpusC <- sapply(which(cmat[,i]==c), function(d){zmat[[d]][i,]})																#extract subset of topic assignments with interaction pattern c
	      if(length(corpusC)>0){zc<-tabulate(unlist(corpusC), nbins=K)} else {zc <- rep(0, K)}												#extract table for topics appear in the document
	        lambdai <- rowSums(sweep(t(sapply(1:length(A), function(j){allxmat[[edge[d,1]]][d,,j]})), 2, bmat[[c]][,i], '*'))				#beta%*%x(i,j,t) given previous beta
	      	const[c]<-log(gamma[c])+lambdai[edge[d,2]]-log(sum(exp(lambdai)))-length(currentZ)*log(sum(zc)-1+alpha)+sum(sapply(currentZ, function(k){log(zc[k]-1+alpha*mvec[k])})) #Eq. (11) 
	      	}
 			cmat[d,i+1]<-which(rmultinom(1, 1, exp(const))==TRUE)																			#Gibbs update for c 
	      }	 
	      
	      #Z update - within each word m in the document d  			
	      corpusCnew<-list()
	      zcnew<-list()
	      for (c in 1:C){ corpusCnew[[c]] <- sapply(which(cmat[,i+1]==c), function(d){zmat[[d]][i,]})										#extract subset of documents with updated c
	      if(length(corpusCnew[[c]])>0){zcnew[[c]]<-tabulate(unlist(corpusCnew[[c]]), nbins=K)} else {zcnew[[c]] <- rep(0, K) }				#extract table for topics appear in the documents (CK)
		  }													     											
	      corpusW <- list()
	      for (k in 1:K){corpusW[[k]] <- unlist(sapply(1:nrow(edge), function(d){textlist[[d]][which(zmat[[d]][i,]==k)]}))}					#extract observed words for each topic k (WK)
	      for (d in 1:nrow(edge)){
	      currentC <- cmat[d, i+1]																											#use the updated c for the document
	      iter <-1	 
	      for (w in zmat[[d]][i,]){	     																									
	      		const2<-sapply(1:K, function(k){(zcnew[[currentC]][k]-1+alpha*mvec[k])*(tabulate(corpusW[[k]], nbins=W)[textlist[[d]][iter]]-1+delta*nvec[iter])/sum(length(corpusW[[k]])-1+delta)}) 																																					#for each topic k, calculate the constants used for multinomial 
	     	zmat[[d]][i+1,iter]<- (which(rmultinom(1, 1, const2)==TRUE))																		#Gibbs update for z
	     	iter <- iter+1																													#iter is used to count mth word in document 
         }
         }
	    
	      #B update - within each interaction pattern C (but using all documents in corpus 1:D)
	      for(c in 1:C){	
	      	edgeC <- edge[which(cmat[,i+1]==c),]																								#extract subset of documents with interaction pattern c
	      	bold <- bmat[[c]][,i]																											#old beta from previous iteration
		  	bnew <- rmvnorm(1, bold, (delta_B)^2*diag(P))																					#new beta using rmvnorm
		  	lambdaiold<- t(sapply(1:nrow(edgeC), function(d){rowSums(sweep(t(sapply(1:length(A), function(j){allxmat[[edge[d,1]]][d,,j]})), 2, bold, '*'))}))			#sum of beta%*%x(i,j,t) given bold
         	lambdainew<- t(sapply(1:nrow(edgeC), function(d){rowSums(sweep(t(sapply(1:length(A), function(j){allxmat[[edge[d,1]]][d,,j]})), 2, bnew, '*'))}))			#sum of beta%*%x(i,j,t) given bnew
		  	loglikediff <- log(dmvnorm(bnew,rep(0,P), sigma^2*diag(P)))-log(dmvnorm(bold,rep(0,P), sigma^2*diag(P)))+	
		  				   sum(sapply(1:nrow(edgeC), function(d){lambdainew[d, edgeC[d,2]]}))-sum(sapply(1:nrow(edgeC), function(d){log(sum(exp(lambdainew)))}))-
		  				   sum(sapply(1:nrow(edgeC), function(d){lambdaiold[d, edgeC[d,2]]}))+sum(sapply(1:nrow(edgeC), function(d){log(sum(exp(lambdaiold)))}))			#difference in loglikelihood
	            u=runif(1) 																																				#compare with unifrom(0,1)
		  if (log(u) < loglikediff){bmat[[c]][,i+1]<-bnew } else {bmat[[c]][,i+1]<-bold}														#M-H update
	      }																													
   		  	  }
		  		
	   return(list(C=cmat, Z=zmat, B=bmat))
       }

unix.time(try2 <- MCMC(10^2, 1)) #about 200 seconds for 100 samples

unix.time(try <- MCMC(2*10^4, 1)) #about 10.36 hours for 20000 samples
save(try, file="try.RData")
sapply(1:269, function(d){tabulate(try$C[d,], nbin=5)})						#table of IP for each document
sapply(1:5, function(c){nrow(unique(t(try$B[[c]])))})/(2*10^4	)	#acceptance rate for beta
sapply(1:5, function(c){colMeans(t(try$B[[c]]))})					#posterior mean for beta 
Zsummary<-list()
for (d in 1:269){
	Zsummary[[d]]<-sapply(1:ncol(try$Z[[d]]), function(k){tabulate(try$Z[[d]][,k], nbins=20)})
	colnames(Zsummary[[d]])<-c(wordlist[[d]])
}
Zsummary																#table of word-topic assignment for each document

