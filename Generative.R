library(gtools)
library(mvtnorm)

setwd("/Users/bomin8319/Desktop/SU16/Rcode")
load("Temporal_Email_Data.Rdata")

#First, let's look at Dare county
Dare<-Temporal_Email_Data$Dare
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
node<-as.data.frame(Dare_node)
set<- 1:27

edge<-matrix(as.numeric(as.matrix(edge)), nrow=nrow(edge), ncol=ncol(edge))



#Assuming 100 topics and 10000 words
K=100
V=10000
phi<-list()
delta<-rep(10,V) #assume delta=10
for (k in 1:K){
	phi[[k]] <- rdirichlet(1, delta)
}


#function to obtain time-dependent covariates (7 dynamic statistics)
x<- function(i,j,t, weight=5*10^(-7)){
intercept <- 1;
send<-sum(dexp(t-(edge[edge[,1]==i & edge[,2]==j & edge[,3] < t, 3]), rate=weight));
receive <-sum(dexp(t-(edge[edge[,1]==j & edge[,2]==i & edge[,3] < t, 3]), rate=weight));
twosend <- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==i & edge[,2]==h & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==h & edge[,2]==j & edge[,3] < t, 3]), rate=weight))}))
tworeceive <- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==h & edge[,2]==i & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==j & edge[,2]==h & edge[,3] < t, 3]), rate=weight))}))
sibling<- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==h & edge[,2]==i & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==h & edge[,2]==j & edge[,3] < t, 3]), rate=weight))}))
cosibling<- sum(sapply(A[!A==i & !A==j], function(h){sum(dexp(t-(edge[edge[,1]==i & edge[,2]==h & edge[,3] < t, 3]), rate=weight))*sum(dexp(t-(edge[edge[,1]==j& edge[,2]==h & edge[,3] < t, 3]), rate=weight))}))
return(c(intercept, send, receive, twosend, tworeceive, sibling, cosibling))}



#Assuming 5 Interaction Patterns
#Use P=7 network statistics as defined in the article
C=5
P=7
sigma <- 1
alpha<-rep(5, K)
beta <- list()
theta <- list()
lambdamat<-list()
	
#after we fix time t	
#e.g. 
t=1354160524
	xmat <- list()
    for(i in 1:length(A)){
	xmat[[i]] <- matrix(0, nrow=length(A), ncol=7)
		for(j in 1:length(A)){
	 if(!i==j){xmat[[i]][j,]<- x(i, j, t)} else {xmat[[i]][j,]<- c(1, rep(0,6))}
		}
	}

lambda <- function(beta, A){
	lambdaij<-matrix(0, nrow=length(A),ncol=length(A))
	for(i in 1:length(A)){
		for(j in 1:length(A)){
			lambdaij[i,j]<-ifelse(!i==j, exp(beta%*%xmat[[i]][j,]), 0)
		}
	}
	return(lambdaij)
}

A <- set
for (c in 1:C){
	beta[[c]]<-rmvnorm(1, rep(0, P), diag(P))
	theta[[c]]<- rdirichlet(1,alpha)
	lambdamat[[c]]<-lambda(beta[[c]], A)
}

#Assuming D=100 documents
D=100 
eta<-rep(10,C) #assume delta=10
gamma <- rdirichlet(1, eta)
document <- list()
for (d in 1:D){
	cdocument<-which(rmultinom(1, 1,gamma)==TRUE)
	time <- rexp(1, sum(lambdamat[[cdocument]]))
	ijcomb <- which(rmultinom(1, 1, c(lambdamat[[cdocument]]/sum(lambdamat[[cdocument]])))==TRUE)
	j <-	 ijcomb %% length(A)	
	i <- (ijcomb-j)/length(A)+1
	document[[d]] <- list(timestamp=t+time, sender=i, receiver=j, IP=cdocument)

	length <- runif(1, 1, 100)	 #randomly choose the number of words in the document
	document[[d]]$words<- rep(0, length)
	for ( m in 1:length){ 
		z <- which(rmultinom(1, 1, theta[[cdocument]])==TRUE)
		w<- which(rmultinom(1, 1, phi[[z]])==TRUE)
		document[[d]]$words[m]<-w }
	
}


#problem: need to update xmat every time
#current version: at time t, based on the history of interaction until t, I want to write d=100 documents at the same time.
#(only epsilon difference and they do not change or update the history thus not affect each other)

