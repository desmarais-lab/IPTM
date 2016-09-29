library(entropy)
#Caclulate the VI between the assignments
sample1<-load("VanceMCMC3.RData")
sample2<-load("VanceMCMC2.RData")
C1 <- VanceMCMC$C
C2 <- VanceMCMC2$C


#Variation of Information for Interactin Patterns & Topics
VIforC <- function(C1, C2){
	C1table <- tabulate(C1)
	C2table <- tabulate(C2)
	sameC <- table(C1, C2)
	return(2*entropy.empirical(sameC)-entropy.empirical(C1table)-entropy.empirical(C2table) )
	}
	

Z1 <- unlist(VanceMCMC$Z)
Z2 <- unlist(VanceMCMC2$Z)

	
VIforC <- function(C1, C2){
	C1table <- tabulate(C1)
	C2table <- tabulate(C2)
	sameC <- table(C1, C2)
	return(2*entropy.empirical(sameC)-entropy.empirical(C1table)-entropy.empirical(C2table) )
	}	
	
	
	
#with hyperparameter optimization
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver21.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver22.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver23.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver24.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver25.RData")
C21<-unlist(VanceMCMC21$C)
C22<-unlist(VanceMCMC22$C)
C23<-unlist(VanceMCMC23$C)
C24<-unlist(VanceMCMC24$C)
C25<-unlist(VanceMCMC25$C)
mean(VIforC(C21, C22),VIforC(C21, C23),VIforC(C21, C24),VIforC(C21, C25),VIforC(C22, C23),VIforC(C22, C24),VIforC(C22, C25),
VIforC(C23, C24),VIforC(C23, C25),VIforC(C24, C25))

#without hyperparameter optimization(fixed)
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver11.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver12.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver13.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver14.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver15.RData")
C11<-unlist(VanceMCMC11$C)
C12<-unlist(VanceMCMC12$C)
C13<-unlist(VanceMCMC13$C)
C14<-unlist(VanceMCMC14$C)
C15<-unlist(VanceMCMC15$C)

mean(VIforC(C11, C12),VIforC(C11, C13),VIforC(C11, C14),VIforC(C11, C15),VIforC(C12, C13),VIforC(C12, C14),VIforC(C12, C15),
VIforC(C13, C14),VIforC(C13, C15),VIforC(C14, C15))

#across the models
without <- rbind(C21, C22, C23, C24, C25)
with <- rbind(C11, C12, C13, C14, C15)
VIs<-c()
for (i in 1:5){
	for (j in 1:5){
		VIs<- c(VIs, VIforC(without[i,], with[j,])) 
	}
}


	
#with hyperparameter optimization
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver21.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver22.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver23.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver24.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver25.RData")
C21<-unlist(VanceMCMC21$Z)
C22<-unlist(VanceMCMC22$Z)
C23<-unlist(VanceMCMC23$Z)
C24<-unlist(VanceMCMC24$Z)
C25<-unlist(VanceMCMC25$Z)
sd(c(VIforC(C21, C22),VIforC(C21, C23),VIforC(C21, C24),VIforC(C21, C25),VIforC(C22, C23),VIforC(C22, C24),VIforC(C22, C25),
VIforC(C23, C24),VIforC(C23, C25),VIforC(C24, C25)))

#without hyperparameter optimization
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver11.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver12.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver13.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver14.RData")
load("/Users/bomin8319/Desktop/SU16/Rcode/data/VanceMCMCserver15.RData")
C11<-unlist(VanceMCMC11$Z)
C12<-unlist(VanceMCMC12$Z)
C13<-unlist(VanceMCMC13$Z)
C14<-unlist(VanceMCMC14$Z)
C15<-unlist(VanceMCMC15$Z)

sd(c(VIforC(C11, C12),VIforC(C11, C13),VIforC(C11, C14),VIforC(C11, C15),VIforC(C12, C13),VIforC(C12, C14),VIforC(C12, C15),
VIforC(C13, C14),VIforC(C13, C15),VIforC(C14, C15)))

#across the models
without <- rbind(C21, C22, C23, C24, C25)
with <- rbind(C11, C12, C13, C14, C15)
VIs<-c()
for (i in 1:5){
	for (j in 1:5){
		if(!i==j){
		VIs<- c(VIs, VIforC(without[i,], with[j,])) }
	}
}


#plot
Lwith<- rbind(VanceMCMC21$L, VanceMCMC22$L, VanceMCMC23$L, VanceMCMC24$L, VanceMCMC25$L)
Lwithout <- rbind(VanceMCMC11$L, VanceMCMC12$L, VanceMCMC13$L, VanceMCMC14$L, VanceMCMC15$L)


plot(colMeans(Lwith), type='l', col='red', xlab="(Outer) Iterations", ylab="Log-Probability")
lines(colMeans(Lwithout), col='blue')
title('Comparison of Symmetric and Asymmetric Models')
legend(locator(1), c("Symmetric", "Asymmetric"), col=c('blue', 'red'), pch=15)


#hist
Awith <- rbind(VanceMCMC21$A, VanceMCMC22$A, VanceMCMC23$A, VanceMCMC24$A, VanceMCMC25$A)
alpha <- c(Awith[,-1])
hist(alpha, xlab="" )

#topic proportions
Mwith <- list()
Mwith[[1]]<-VanceMCMC21$M[200,]
Mwith[[2]]<-VanceMCMC22$M[200,]
Mwith[[3]]<-VanceMCMC23$M[200,]
Mwith[[4]]<-VanceMCMC24$M[200,]
Mwith[[5]]<-VanceMCMC25$M[200,]

#entropy
Ewith<-rbind(VanceMCMC21$E[,10],VanceMCMC22$E[,10],VanceMCMC23$E[,10],VanceMCMC24$E[,10],VanceMCMC25$E[,10])
Ewithout<-rbind(VanceMCMC11$E[,10],VanceMCMC12$E[,10],VanceMCMC13$E[,10],VanceMCMC14$E[,10],VanceMCMC15$E[,10])
plot(colMeans(Ewith), type='l')
lines(colMeans(Ewithout), col='blue')



#Now summarize the topic distributions based on the best result -> Second







par(mfrow=c(2,1))

#Summary of Z
CofD<- VanceMCMC12$C	
Zsummary<-VanceMCMC12$Z
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
title('Topic Distritubitions given IP1-IP3 (Symmetric)')
legend(locator(1), c("IP1","IP2","IP3"), col=gray.colors(3), pch=15)



CofD<- VanceMCMC22$C	
Zsummary<-VanceMCMC22$Z
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
title('Topic Distritubitions given IP1-IP3 (Asymmetric)')
legend(locator(1), c("IP1","IP2","IP3"), col=gray.colors(3), pch=15)



#words

CofD<- VanceMCMC12$C	
Zsummary<-VanceMCMC12$Z
for (d in 1:269){
	names(Zsummary[[d]])<-wordlist[[d]]
}
Ztable<-tabulate(unlist(Zsummary), nbins=20)														#table of word-topic assignment for each document
names(Ztable)<-c(1:20)
Ztable<-rbind(Ztable[order(Ztable, decreasing=TRUE)], Ztable[order(Ztable, decreasing=TRUE)]/sum(Ztable))

allwords<-unlist(Zsummary)
Wsummary<-sapply(1:20, function(k){allwords[allwords==k]})
Wtable<-sapply(1:20, function(k){table(names(Wsummary[[k]]))})

library(lda)
K=20
W=620
topicword <- matrix(0, nrow=K, ncol=W)
colnames(topicword)<- vocabulary
for (i in 1:length(allwords)){
	topicword[allwords[i], which(colnames(topicword)==names(allwords[i]))]<-topicword[allwords[i], which(colnames(topicword)==names(allwords[i]))]+1
}
top.topic.words(topicword, num.words=10, by.score=TRUE)



CofD<- VanceMCMC22$C	
Zsummary<-VanceMCMC22$Z
for (d in 1:269){
	names(Zsummary[[d]])<-wordlist[[d]]
}
Ztable<-tabulate(unlist(Zsummary), nbins=20)														#table of word-topic assignment for each document
names(Ztable)<-c(1:20)
Ztable<-rbind(Ztable[order(Ztable, decreasing=TRUE)], Ztable[order(Ztable, decreasing=TRUE)]/sum(Ztable))

hi <- rbind(Ztable, Ztable2)/sum(Ztable)
barplot(hi, beside=TRUE)
title('Topic Distritubitions combining IP1-IP3')
legend(locator(1), c("Symmetric","Asymmetric"), col=gray.colors(2), pch=15)
